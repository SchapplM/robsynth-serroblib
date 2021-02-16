% Return the minimum parameter vector for
% S5RRPRR8
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
% Datum: 2021-01-15 21:36
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPRR8_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR8_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR8_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t7 = m(5) + m(6);
t2 = t7 * pkin(3) ^ 2;
t14 = t2 - Ifges(4,1) + Ifges(4,2);
t6 = cos(pkin(9));
t3 = t6 ^ 2;
t19 = t14 * t3;
t5 = sin(pkin(9));
t18 = t5 * t6;
t17 = 2 * pkin(8) * mrSges(6,3) + Ifges(6,2);
t15 = Ifges(4,4) * t18;
t13 = pkin(8) * m(6) + mrSges(6,3);
t12 = -t7 * pkin(7) - mrSges(5,3);
t10 = (pkin(4) ^ 2);
t9 = (pkin(7) ^ 2);
t8 = pkin(8) ^ 2;
t1 = t12 * pkin(3) + Ifges(4,5);
t4 = [t19 + 0.2e1 * t15 + ((t9 + t10) * m(6)) + (m(5) * t9) + (2 * pkin(7) * mrSges(5,3)) + Ifges(4,1) + ((pkin(1) ^ 2 + pkin(6) ^ 2) * m(3)) + (2 * pkin(6) * mrSges(3,3)) + Ifges(3,2) + Ifges(5,2) + Ifges(2,3); pkin(1) * m(3) + mrSges(2,1); -pkin(6) * m(3) + mrSges(2,2) - mrSges(3,3); Ifges(3,1) - Ifges(3,2) + t14 - 0.4e1 * t15 - 0.2e1 * t19; 0.2e1 * t3 * Ifges(4,4) - t14 * t18 + Ifges(3,4) - Ifges(4,4); -t5 * Ifges(4,6) + t1 * t6 + Ifges(3,5); t6 * Ifges(4,6) + t1 * t5 + Ifges(3,6); t2 + Ifges(3,3) + Ifges(4,3); mrSges(3,1); mrSges(3,2); t7 * pkin(3) + mrSges(4,1); mrSges(4,2); mrSges(4,3) - t12; m(4) + t7; Ifges(5,1) - Ifges(5,2) + ((t8 - t10) * m(6)) + t17; t13 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t8 + t10) * m(6) + t17; pkin(4) * m(6) + mrSges(5,1); mrSges(5,2) - t13; Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t4;
