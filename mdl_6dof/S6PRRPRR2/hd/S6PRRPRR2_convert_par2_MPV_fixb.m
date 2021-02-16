% Return the minimum parameter vector for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [29x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:49
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRPRR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR2_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t12 = (mrSges(6,3) + mrSges(7,3));
t25 = (-(-pkin(10) - pkin(9)) * m(7) + pkin(9) * m(6));
t38 = t25 + t12;
t42 = t38 * pkin(4) + Ifges(5,4);
t17 = (pkin(9) ^ 2);
t29 = (-Ifges(6,2) - Ifges(7,2));
t36 = 2 * pkin(9);
t41 = m(6) * t17 + t12 * t36 - t29;
t11 = cos(pkin(12));
t6 = t11 ^ 2;
t40 = 0.2e1 * t6;
t39 = Ifges(5,1) - Ifges(5,2) + t41;
t19 = (pkin(4) ^ 2);
t15 = m(6) * t19;
t18 = (pkin(5) ^ 2);
t28 = (pkin(10) ^ 2 + t18);
t23 = -2 * pkin(9) * pkin(10) - t17 - t28;
t37 = -(t19 + t23) * m(7) - t15 + t39;
t35 = m(7) * pkin(10);
t14 = m(6) + m(7);
t33 = mrSges(7,3) * pkin(10);
t32 = t18 * m(7);
t31 = t12 * pkin(4) + Ifges(5,4);
t10 = sin(pkin(12));
t30 = t10 * t11;
t7 = 2 * t33;
t27 = -mrSges(7,3) - t35;
t22 = -t14 * t17 - Ifges(5,1) + t29;
t8 = -2 * t33;
t1 = [m(2) + m(3) + m(4); ((t19 - t28) * m(7) + t15 + t8 + Ifges(5,2) + t22) * t6 + 0.2e1 * ((pkin(9) * t14 + t35) * pkin(4) + t31) * t30 + (t28 * m(7)) + t7 + ((pkin(2) ^ 2 + pkin(8) ^ 2) * m(4)) + (2 * pkin(8) * mrSges(4,3)) + Ifges(4,2) + Ifges(3,3) + (-t6 + 0.1e1) * (mrSges(6,3) - t27) * t36 - t22; pkin(2) * m(4) + mrSges(3,1); -pkin(8) * m(4) + mrSges(3,2) - mrSges(4,3); (t7 + t37) * t40 - 0.4e1 * (t25 * pkin(4) + t31) * t30 + t8 - Ifges(4,2) + Ifges(4,1) - t37; t42 * t40 - (t23 * m(7) + t14 * t19 - t39 + t8) * t30 + Ifges(4,4) - t42; t11 * Ifges(5,5) - t10 * Ifges(5,6) + Ifges(4,5); t10 * Ifges(5,5) + t11 * Ifges(5,6) + Ifges(4,6); (t19 - t23) * m(7) + t15 + t7 + Ifges(4,3) + Ifges(5,3) + t41; mrSges(4,1); mrSges(4,2); t14 * pkin(4) + mrSges(5,1); mrSges(5,2) - t38; mrSges(5,3); m(5) + t14; Ifges(6,1) - Ifges(6,2) - t32; Ifges(6,4); t27 * pkin(5) + Ifges(6,5); Ifges(6,6); Ifges(6,3) + t32; pkin(5) * m(7) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV = t1;
