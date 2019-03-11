% Return the minimum parameter vector for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% m_mdh [7x1]
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
% MPV [33x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t132 = (m(6) + m(7));
t152 = -pkin(9) * t132 - mrSges(6,3);
t127 = (m(5) + t132);
t136 = (pkin(8) ^ 2);
t139 = (pkin(3) ^ 2);
t121 = (mrSges(5,3) - t152);
t137 = (pkin(5) ^ 2);
t147 = (t137 * m(7) + Ifges(6,2));
t135 = (pkin(9) ^ 2);
t138 = (pkin(4) ^ 2);
t148 = (Ifges(5,2) + (t135 + t138) * t132);
t140 = 2 * pkin(9) * mrSges(6,3) + 2 * pkin(8) * t121 + t147 + t148;
t116 = Ifges(4,2) + (t136 + t139) * t127 + t140;
t117 = t136 * t127 + Ifges(4,1) + t140;
t151 = -t116 + t117;
t150 = -pkin(8) * t127 - t121;
t146 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t130 = sin(pkin(11));
t131 = cos(pkin(11));
t145 = t130 * t131;
t125 = t130 ^ 2;
t126 = t131 ^ 2;
t144 = t126 - t125;
t143 = Ifges(4,4) * t145;
t142 = pkin(10) * m(7) + mrSges(7,3);
t122 = pkin(3) * t127 + mrSges(4,1);
t141 = -t130 * mrSges(4,2) + t131 * t122;
t134 = pkin(10) ^ 2;
t118 = t150 * pkin(3) + Ifges(4,5);
t1 = [Ifges(2,3) + Ifges(3,2) + t125 * t117 + 0.2e1 * t143 + t126 * t116 + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); t151 * t144 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t143; t144 * Ifges(4,4) + t151 * t145 + Ifges(3,4); -t130 * Ifges(4,6) + t131 * t118 + Ifges(3,5); t131 * Ifges(4,6) + t130 * t118 + Ifges(3,6); 0.2e1 * pkin(2) * t141 + (t139 * t127) + Ifges(3,3) + Ifges(4,3); mrSges(3,1) + t141; t131 * mrSges(4,2) + t130 * t122 + mrSges(3,2); mrSges(4,3) - t150; m(4) + t127; t135 * t132 + Ifges(5,1) - t148; Ifges(5,4); t152 * pkin(4) + Ifges(5,5); Ifges(5,6); t138 * t132 + Ifges(5,3); pkin(4) * t132 + mrSges(5,1); mrSges(5,2); m(7) * t134 + Ifges(6,1) + t146 - t147; t142 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t134 + t137) * m(7) + t146; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t142; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
