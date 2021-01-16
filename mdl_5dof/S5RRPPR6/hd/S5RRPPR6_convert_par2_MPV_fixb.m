% Return the minimum parameter vector for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% MPV [25x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPPR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR6_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t39 = 2 * Ifges(5,4);
t16 = cos(pkin(8));
t10 = t16 ^ 2;
t38 = 0.2e1 * t10;
t13 = sin(pkin(9));
t15 = cos(pkin(9));
t8 = pkin(7) * m(6) + mrSges(6,3);
t6 = t8 * pkin(4) - Ifges(5,5);
t37 = t13 * Ifges(5,6) + t6 * t15 + Ifges(4,4);
t26 = t13 * t15;
t18 = pkin(4) ^ 2;
t7 = m(6) * t18 - Ifges(5,1) + Ifges(5,2);
t9 = t15 ^ 2;
t30 = t26 * t39 + t7 * t9;
t36 = Ifges(4,1) + Ifges(5,2) - t30;
t33 = -Ifges(4,2) - Ifges(5,3);
t28 = (pkin(7) * mrSges(6,3));
t11 = 2 * t28;
t32 = t11 + Ifges(6,2) + t36;
t17 = pkin(7) ^ 2;
t29 = m(6) * t17;
t14 = sin(pkin(8));
t25 = t14 * t16;
t24 = t37 * t25;
t22 = Ifges(6,2) + t29;
t20 = t15 * Ifges(5,6) - t6 * t13 - Ifges(4,6);
t1 = -(2 * t28) - t22 - t33 - t36;
t2 = -t7 * t26 + t9 * t39 - Ifges(5,4) + Ifges(4,5);
t3 = [t1 * t10 + 0.2e1 * t24 + (t17 + t18) * m(6) + ((pkin(1) ^ 2 + pkin(6) ^ 2) * m(3)) + (2 * pkin(6) * mrSges(3,3)) + Ifges(3,2) + Ifges(2,3) + t32; pkin(1) * m(3) + mrSges(2,1); -pkin(6) * m(3) + mrSges(2,2) - mrSges(3,3); (t29 + t32 + t33) * t38 - 0.4e1 * t24 - Ifges(3,2) + Ifges(3,1) + t1; -t1 * t25 + t37 * t38 + Ifges(3,4) - t37; t20 * t14 + t2 * t16 + Ifges(3,5); t2 * t14 - t20 * t16 + Ifges(3,6); Ifges(5,1) + Ifges(3,3) + Ifges(4,3) + t11 + t22 + t30; mrSges(3,1); mrSges(3,2); mrSges(4,1); mrSges(4,2); mrSges(4,3); m(4); pkin(4) * m(6) + mrSges(5,1); mrSges(5,2); mrSges(5,3) + t8; m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t3;
