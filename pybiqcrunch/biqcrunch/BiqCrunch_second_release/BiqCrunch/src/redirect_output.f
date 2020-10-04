      subroutine redirect_output()
                close(0)
                close(6)
                open(0,FILE = '/dev/null')
                open(6,FILE = '/dev/null')
        write(6,*) ii

        return
        end
