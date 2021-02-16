% Return the minimum parameter vector for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [28x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPRR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR2_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t7 = (m(5) + m(6));
t2 = t7 * pkin(3) ^ 2;
t12 = t2 - Ifges(4,1) + Ifges(4,2);
t5 = cos(pkin(9));
t3 = t5 ^ 2;
t19 = t12 * t3;
t18 = 2 * pkin(7);
t4 = sin(pkin(9));
t17 = t4 * t5;
t9 = (pkin(4) ^ 2);
t16 = t9 * m(6);
t15 = (mrSges(5,3) + mrSges(6,3));
t13 = Ifges(4,4) * t17;
t11 = -pkin(7) * m(5) - (pkin(7) + pkin(8)) * m(6) - t15;
t8 = pkin(7) ^ 2;
t1 = t11 * pkin(3) + Ifges(4,5);
t6 = [t19 + 0.2e1 * t13 + ((t8 + t9 + (t18 + pkin(8)) * pkin(8)) * m(6)) + (m(5) * t8) + (t15 * t18) + (2 * mrSges(6,3) * pkin(8)) + Ifges(4,1) + ((pkin(1) ^ 2 + pkin(6) ^ 2) * m(3)) + (2 * pkin(6) * mrSges(3,3)) + Ifges(3,2) + Ifges(5,2) + Ifges(6,2) + Ifges(2,3); pkin(1) * m(3) + mrSges(2,1); -pkin(6) * m(3) + mrSges(2,2) - mrSges(3,3); Ifges(3,1) - Ifges(3,2) + t12 - 0.4e1 * t13 - 0.2e1 * t19; 0.2e1 * t3 * Ifges(4,4) - t12 * t17 + Ifges(3,4) - Ifges(4,4); -t4 * Ifges(4,6) + t1 * t5 + Ifges(3,5); t5 * Ifges(4,6) + t1 * t4 + Ifges(3,6); t2 + Ifges(3,3) + Ifges(4,3); mrSges(3,1); mrSges(3,2); t7 * pkin(3) + mrSges(4,1); mrSges(4,2); mrSges(4,3) - t11; m(4) + t7; Ifges(5,1) - Ifges(5,2) - t16; Ifges(5,4); (-m(6) * pkin(8) - mrSges(6,3)) * pkin(4) + Ifges(5,5); Ifges(5,6); Ifges(5,3) + t16; pkin(4) * m(6) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t6;
