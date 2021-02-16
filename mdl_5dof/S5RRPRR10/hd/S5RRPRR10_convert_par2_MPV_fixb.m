% Return the minimum parameter vector for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2021-01-15 22:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPRR10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t14 = m(5) + m(6);
t17 = pkin(8) ^ 2;
t42 = (-t14 * t17 - Ifges(4,1) - Ifges(5,2));
t35 = t14 * pkin(8) + mrSges(5,3);
t34 = t35 * pkin(3) + Ifges(4,4);
t12 = cos(pkin(10));
t7 = t12 ^ 2;
t41 = 0.2e1 * t7;
t40 = (m(3) * pkin(7));
t39 = (Ifges(4,2) + t42);
t20 = pkin(3) ^ 2;
t15 = m(5) * t20;
t19 = pkin(4) ^ 2;
t38 = (-t15 - (-t19 + t20) * m(6));
t32 = (m(6) * t19);
t36 = -Ifges(3,2) - t32 - (2 * mrSges(3,3) + t40) * pkin(7) + t42;
t31 = (pkin(8) * mrSges(5,3));
t30 = 2 * pkin(9) * mrSges(6,3) + Ifges(6,2);
t11 = sin(pkin(10));
t28 = t11 * t12;
t9 = 2 * t31;
t26 = t34 * t28;
t24 = pkin(9) * m(6) + mrSges(6,3);
t10 = -2 * t31;
t23 = t10 + t39;
t16 = pkin(9) ^ 2;
t13 = cos(pkin(5));
t1 = [(t10 + t36) * t13 ^ 2 + t9 + m(3) * pkin(1) ^ 2 + Ifges(2,3) + 0.2e1 * (((m(5) / 0.2e1 + m(6) / 0.2e1) * t17 + t31 + (-t20 / 0.2e1 + t19 / 0.2e1) * m(6) - t15 / 0.2e1 + Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t7 - t26) * (t13 - 0.1e1) * (t13 + 0.1e1) - t36; pkin(1) * m(3) + mrSges(2,1); (-mrSges(3,3) - t40) * sin(pkin(5)) + mrSges(2,2); -0.4e1 * t26 - Ifges(3,2) + Ifges(3,1) + t23 + (t9 + t38 - t39) * t41 - t38; t34 * t41 - (t14 * t20 + t23 - t32) * t28 + Ifges(3,4) - t34; Ifges(4,5) * t12 - Ifges(4,6) * t11 + Ifges(3,5); Ifges(4,5) * t11 + Ifges(4,6) * t12 + Ifges(3,6); (t17 + t19 + t20) * m(6) + t15 + m(5) * t17 + t9 + Ifges(5,2) + Ifges(3,3) + Ifges(4,3); mrSges(3,1); mrSges(3,2); pkin(3) * t14 + mrSges(4,1); mrSges(4,2) - t35; mrSges(4,3); m(4) + t14; Ifges(5,1) - Ifges(5,2) + (t16 - t19) * m(6) + t30; pkin(4) * t24 + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t16 + t19) * m(6) + t30; pkin(4) * m(6) + mrSges(5,1); mrSges(5,2) - t24; Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
