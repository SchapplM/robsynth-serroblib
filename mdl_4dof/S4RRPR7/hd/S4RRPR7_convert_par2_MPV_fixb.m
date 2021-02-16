% Return the minimum parameter vector for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [21x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S4RRPR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_convert_par2_MPV_fixb: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR7_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR7_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t10 = pkin(6) ^ 2;
t11 = pkin(3) ^ 2;
t17 = Ifges(4,1) + Ifges(5,2);
t12 = (Ifges(4,2) - t17 + (-t10 + t11) * m(5));
t18 = (pkin(6) * mrSges(5,3));
t1 = -2 * t18 + t12;
t3 = pkin(6) * m(5) + mrSges(5,3);
t13 = t3 * pkin(3) + Ifges(4,4);
t8 = sin(pkin(7));
t9 = cos(pkin(7));
t19 = t8 * t9;
t16 = t13 * t19;
t6 = 2 * t18;
t5 = t9 ^ 2;
t15 = 0.2e1 * t5;
t2 = [t1 * t5 + 0.2e1 * t16 + m(5) * t10 + t6 + ((pkin(1) ^ 2 + pkin(5) ^ 2) * m(3)) + (2 * pkin(5) * mrSges(3,3)) + Ifges(3,2) + Ifges(2,3) + t17; pkin(1) * m(3) + mrSges(2,1); -pkin(5) * m(3) + mrSges(2,2) - mrSges(3,3); (-t12 + t6) * t15 - 0.4e1 * t16 + Ifges(3,1) - Ifges(3,2) + t1; -t1 * t19 + t13 * t15 + Ifges(3,4) - t13; t9 * Ifges(4,5) - t8 * Ifges(4,6) + Ifges(3,5); t8 * Ifges(4,5) + t9 * Ifges(4,6) + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + Ifges(5,2) + t6 + (t10 + t11) * m(5); mrSges(3,1); mrSges(3,2); pkin(3) * m(5) + mrSges(4,1); mrSges(4,2) - t3; mrSges(4,3); m(4) + m(5); Ifges(5,1) - Ifges(5,2); Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3); mrSges(5,1); mrSges(5,2);];
MPV = t2;
