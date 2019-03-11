% Return the minimum parameter vector for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% MPV [30x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t141 = (pkin(8) ^ 2);
t143 = (pkin(3) ^ 2);
t142 = (pkin(5) ^ 2);
t150 = (t142 * m(7) + Ifges(5,2) + Ifges(6,3));
t145 = 2 * pkin(8) * mrSges(5,3) + t150;
t122 = Ifges(4,2) + (t141 + t143) * m(5) + t145;
t123 = t141 * m(5) + Ifges(4,1) + t145;
t155 = -t122 + t123;
t154 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t135 = sin(pkin(11));
t137 = cos(pkin(11));
t153 = t135 * t137;
t136 = sin(pkin(10));
t138 = cos(pkin(10));
t152 = t136 * t138;
t130 = t136 ^ 2;
t132 = t138 ^ 2;
t151 = t132 - t130;
t149 = Ifges(6,4) * t153;
t148 = Ifges(4,4) * t152;
t147 = -pkin(8) * m(5) - mrSges(5,3);
t146 = -pkin(9) * m(7) - mrSges(7,3);
t128 = m(5) * pkin(3) + mrSges(4,1);
t144 = -t136 * mrSges(4,2) + t138 * t128;
t140 = pkin(9) ^ 2;
t131 = t137 ^ 2;
t129 = t135 ^ 2;
t127 = t147 * pkin(3) + Ifges(4,5);
t126 = t146 * pkin(5) + Ifges(6,5);
t125 = m(7) * t140 + Ifges(6,1) + t154;
t124 = Ifges(6,2) + (t140 + t142) * m(7) + t154;
t1 = [Ifges(2,3) + Ifges(3,2) + t130 * t123 + 0.2e1 * t148 + t132 * t122 + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); t155 * t151 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t148; t151 * Ifges(4,4) + t155 * t152 + Ifges(3,4); -t136 * Ifges(4,6) + t138 * t127 + Ifges(3,5); t138 * Ifges(4,6) + t136 * t127 + Ifges(3,6); (t143 * m(5)) + 0.2e1 * pkin(2) * t144 + Ifges(3,3) + Ifges(4,3); mrSges(3,1) + t144; t138 * mrSges(4,2) + t136 * t128 + mrSges(3,2); mrSges(4,3) - t147; m(4) + m(5); t129 * t124 + t131 * t125 + Ifges(5,1) - 0.2e1 * t149 - t150; t135 * Ifges(6,6) - t137 * t126 + Ifges(5,4); Ifges(5,5) + (t131 - t129) * Ifges(6,4) + (-t124 + t125) * t153; -t137 * Ifges(6,6) - t135 * t126 + Ifges(5,6); t131 * t124 + t129 * t125 + Ifges(5,3) + 0.2e1 * t149; mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t146; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
