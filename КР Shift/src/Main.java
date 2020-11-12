import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Scanner;
import java.util.StringTokenizer;

public class Main {
    public static void main(String[] args) {
        try {
            Scanner in1 = new Scanner(new File("input1.txt"));
            Scanner in2 = new Scanner(new File("input2.txt"));

            ArrayList<String> allFilesInFolder = new ArrayList<>();
            while (in1.hasNextLine()){
                allFilesInFolder.add(in1.nextLine());
            }
            System.out.println(allFilesInFolder);

            ArrayList<String> deletedFiles = new ArrayList<>();
            ArrayList<Integer> IndexesOfAllocation = new ArrayList<>();
            ArrayList<String> Comands = new ArrayList<>();
            IndexesOfAllocation.add(0);
            
            
            while(in2.hasNextLine()){
                String Str = in2.nextLine();
                if(Str.equals("delete")){
                    int lastAllocation = IndexesOfAllocation.size()-1;
                    for(int i = 0; i < IndexesOfAllocation.size(); i++){
                        if(IndexesOfAllocation.size() > 1){
                            if(IndexesOfAllocation.get(0) < IndexesOfAllocation.get(1)){
                                deletedFiles.add(allFilesInFolder.get(IndexesOfAllocation.get(0)));    //down
                                allFilesInFolder.remove(IndexesOfAllocation.get(0).intValue());
                            }
                            if(IndexesOfAllocation.get(0) > IndexesOfAllocation.get(1)){
                                deletedFiles.add(allFilesInFolder.get(IndexesOfAllocation.get(lastAllocation)));    //up
                                allFilesInFolder.remove(IndexesOfAllocation.get(lastAllocation).intValue());
                            }
                        }
                        else{
                            deletedFiles.add(allFilesInFolder.get(IndexesOfAllocation.get(0)));
                            allFilesInFolder.remove(IndexesOfAllocation.get(0).intValue());
                        }
                    }

                    int currentIndex = IndexesOfAllocation.get(0);
                    if(IndexesOfAllocation.size() > 1){
                        if(IndexesOfAllocation.get(0) < IndexesOfAllocation.get(1)){
                            currentIndex = IndexesOfAllocation.get(0);
                        }
                        if(IndexesOfAllocation.get(0) > IndexesOfAllocation.get(1)){
                            currentIndex = IndexesOfAllocation.get(IndexesOfAllocation.size() - 1);
                        }
                    }
                    IndexesOfAllocation.clear();
                    IndexesOfAllocation.add(currentIndex);
                    continue;
                }
                if(Str.equals("up")){
                    int Index;
                    Index =  IndexesOfAllocation.get(IndexesOfAllocation.size()-1) - 1;
                    IndexesOfAllocation.clear();
                    IndexesOfAllocation.add(Index);
                    continue;
                }
                if(Str.equals("down")){
                    int Index;
                    Index =  IndexesOfAllocation.get(IndexesOfAllocation.size()-1) + 1;
                    IndexesOfAllocation.clear();
                    IndexesOfAllocation.add(Index);
                    continue;
                }
                if(Str.substring(0, 6).equals("shift+")){
                    Comands.add(Str);
                    StringTokenizer st = new StringTokenizer(Str, "+");
                    st.nextToken();
                    int up = 0;
                    int down = 0;
                    while(st.hasMoreTokens()){
                        String tmp;
                        tmp = st.nextToken();
                        if(tmp.equals("up")){
                            up++;
                            continue;
                        }
                        if(tmp.equals("down")){
                            down++;
                        }
                    }
                    int res = down - up;
                    if(res > 0){
                        for(int i = 0; i < res; i++){
                            IndexesOfAllocation.add(IndexesOfAllocation.get(IndexesOfAllocation.size()-1) + 1);
                        }
                    }
                    if(res < 0){
                        res *= -1;
                        for(int i = 0; i < res; i++){
                            IndexesOfAllocation.add(IndexesOfAllocation.get(IndexesOfAllocation.size()-1) - 1);
                        }
                    }
                }
            }

            Comands.stream().sorted(Comparator.comparing(String::length)).forEach(s-> System.out.println(s));
            System.out.println(allFilesInFolder);
            System.out.println(deletedFiles);
        }
        catch (FileNotFoundException e){
            System.out.println("File is not found.");
        }
    }
}