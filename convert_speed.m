function out = convert_speed(val, direction)
% convert_speed - converts between m/s and knots
% direction: 'ms2kt' or 'kt2ms'

switch direction
    case 'ms2kt'
        out = val * 1.94384;
    case 'kt2ms'
        out = val / 1.94384;
    otherwise
        error('direction must be ms2kt or kt2ms')
end
end