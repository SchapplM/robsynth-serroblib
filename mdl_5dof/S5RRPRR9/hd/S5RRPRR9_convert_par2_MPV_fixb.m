% Return the minimum parameter vector for
% S5RRPRR9
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
% Datum: 2021-01-15 21:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPRR9_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR9_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t12 = (mrSges(5,3) + mrSges(6,3));
t25 = (-(-pkin(8) - pkin(7)) * m(6) + pkin(7) * m(5));
t38 = t25 + t12;
t42 = t38 * pkin(3) + Ifges(4,4);
t17 = (pkin(7) ^ 2);
t29 = (-Ifges(5,2) - Ifges(6,2));
t36 = 2 * pkin(7);
t41 = m(5) * t17 + t12 * t36 - t29;
t11 = cos(pkin(9));
t6 = t11 ^ 2;
t40 = 0.2e1 * t6;
t39 = Ifges(4,1) - Ifges(4,2) + t41;
t19 = (pkin(3) ^ 2);
t15 = m(5) * t19;
t18 = (pkin(4) ^ 2);
t28 = (pkin(8) ^ 2 + t18);
t23 = -2 * pkin(7) * pkin(8) - t17 - t28;
t37 = -(t19 + t23) * m(6) - t15 + t39;
t35 = m(6) * pkin(8);
t14 = m(5) + m(6);
t33 = mrSges(6,3) * pkin(8);
t32 = t18 * m(6);
t31 = t12 * pkin(3) + Ifges(4,4);
t10 = sin(pkin(9));
t30 = t10 * t11;
t7 = 2 * t33;
t27 = -mrSges(6,3) - t35;
t22 = -t14 * t17 - Ifges(4,1) + t29;
t8 = -2 * t33;
t1 = [((t19 - t28) * m(6) + t15 + t8 + Ifges(4,2) + t22) * t6 + 0.2e1 * ((pkin(7) * t14 + t35) * pkin(3) + t31) * t30 + (t28 * m(6)) + t7 + ((pkin(1) ^ 2 + pkin(6) ^ 2) * m(3)) + (2 * pkin(6) * mrSges(3,3)) + Ifges(3,2) + Ifges(2,3) + (-t6 + 0.1e1) * (mrSges(5,3) - t27) * t36 - t22; pkin(1) * m(3) + mrSges(2,1); -pkin(6) * m(3) + mrSges(2,2) - mrSges(3,3); (t7 + t37) * t40 - 0.4e1 * (t25 * pkin(3) + t31) * t30 + t8 - Ifges(3,2) + Ifges(3,1) - t37; t42 * t40 - (t23 * m(6) + t14 * t19 - t39 + t8) * t30 + Ifges(3,4) - t42; t11 * Ifges(4,5) - t10 * Ifges(4,6) + Ifges(3,5); t10 * Ifges(4,5) + t11 * Ifges(4,6) + Ifges(3,6); (t19 - t23) * m(6) + t15 + t7 + Ifges(3,3) + Ifges(4,3) + t41; mrSges(3,1); mrSges(3,2); t14 * pkin(3) + mrSges(4,1); mrSges(4,2) - t38; mrSges(4,3); m(4) + t14; Ifges(5,1) - Ifges(5,2) - t32; Ifges(5,4); t27 * pkin(4) + Ifges(5,5); Ifges(5,6); Ifges(5,3) + t32; pkin(4) * m(6) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
