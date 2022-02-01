% Return the minimum parameter vector for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% MPV [23x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPPRR4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t9 = (mrSges(5,3) + mrSges(6,3));
t22 = (pkin(6) + pkin(7)) * m(6) + m(5) * pkin(6) + t9;
t21 = -2 * pkin(7);
t6 = sin(pkin(8));
t20 = (mrSges(4,3) + t22) * t6;
t19 = m(5) + m(6);
t13 = (pkin(4) ^ 2);
t18 = t13 * m(6);
t12 = (pkin(6) ^ 2);
t16 = pkin(6) * t21 - pkin(7) ^ 2 - t12 - t13;
t14 = pkin(3) ^ 2;
t5 = sin(pkin(9));
t7 = cos(pkin(9));
t15 = -(m(5) * t12) + (mrSges(6,3) * t21) - (2 * t9 * pkin(6)) - Ifges(3,1) - Ifges(4,2) - Ifges(5,2) - Ifges(6,2) + (0.2e1 * t5 * Ifges(4,4) + (t14 * t19 - Ifges(4,1) + Ifges(4,2)) * t7) * t7;
t8 = cos(pkin(8));
t4 = 0.1e1 / t8;
t1 = [(t14 - t16) * m(6) + m(5) * t14 + Ifges(2,3) - t15 + (0.2e1 * (t5 * Ifges(4,6) + Ifges(3,4) + (pkin(3) * t22 - Ifges(4,5)) * t7) * t6 + (t16 * m(6) + Ifges(3,2) + Ifges(4,3) + t15) * t8) * t8; mrSges(2,1); mrSges(2,2); (mrSges(3,1) * t8 - mrSges(3,2) * t6) * t4; mrSges(3,3); m(3); (t7 * t20 + t8 * (pkin(3) * t19 + mrSges(4,1))) * t4; (mrSges(4,2) * t8 - t20 * t5) * t4; m(4) + t19; Ifges(5,1) - Ifges(5,2) - t18; Ifges(5,4); (-pkin(7) * m(6) - mrSges(6,3)) * pkin(4) + Ifges(5,5); Ifges(5,6); Ifges(5,3) + t18; pkin(4) * m(6) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
