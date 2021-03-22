module file_date_mod

    use namelist_mod

    public date_advance

    contains

    subroutine date_advance(year , month , day , hour , minute)

        integer, intent(inout) :: year
        integer, intent(inout) :: month
        integer, intent(inout) :: day
        integer, intent(inout) :: hour
        integer, intent(inout) :: minute

        integer month_day(12) 
        
        data month_day /31,28,31,30,31,30,31,31,30,31,30,31/

        if (mod(year , 400) == 0) then 
          month_day(2) = 29
        else if (mod(year , 100) /= 0 .and. mod(year , 4) ==0 ) then
          month_day(2) = 29
        end if

        minute = minute + initial_interval 
        if (minute > 59) then
          minute = minute - 60
          hour = hour + 1
          if (hour > 23) then
            hour = hour - 24
            day = day + 1
            if (day > month_day(month)) then
              day = 1
              month = month + 1
              if (month > 12) then
                month = 1
                year = year + 1
              end if
            end if
          end if
        end if
    end subroutine date_advance

end module file_date_mod