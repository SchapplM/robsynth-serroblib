% Return the minimum parameter vector for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t141 = (pkin(9) ^ 2);
t143 = (pkin(5) ^ 2);
t154 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t126 = Ifges(6,2) + (t141 + t143) * m(7) + t154;
t127 = m(7) * t141 + Ifges(6,1) + t154;
t136 = sin(pkin(11));
t131 = t136 ^ 2;
t138 = cos(pkin(11));
t133 = t138 ^ 2;
t153 = t136 * t138;
t150 = Ifges(6,4) * t153;
t123 = t131 * t126 + t133 * t127 + Ifges(5,1) - 0.2e1 * t150;
t129 = t143 * m(7) + Ifges(5,2) + Ifges(6,3);
t155 = t123 - t129;
t137 = sin(pkin(10));
t139 = cos(pkin(10));
t152 = t137 * t139;
t132 = t137 ^ 2;
t134 = t139 ^ 2;
t151 = t134 - t132;
t149 = -pkin(8) * m(4) - mrSges(4,3);
t148 = -pkin(9) * m(7) - mrSges(7,3);
t128 = t148 * pkin(5) + Ifges(6,5);
t125 = t136 * Ifges(6,6) - t138 * t128 + Ifges(5,4);
t147 = t125 * t152;
t146 = (mrSges(3,3) - t149);
t145 = t139 * mrSges(5,1) - t137 * mrSges(5,2);
t144 = pkin(2) ^ 2;
t142 = pkin(8) ^ 2;
t140 = (m(3) + m(4));
t130 = t142 + t144;
t124 = -t138 * Ifges(6,6) - t136 * t128 + Ifges(5,6);
t122 = Ifges(5,5) + (t133 - t131) * Ifges(6,4) + (-t126 + t127) * t153;
t1 = [Ifges(2,3) + Ifges(3,2) + Ifges(4,2) + t132 * t123 + 0.2e1 * t147 + t134 * t129 + (2 * pkin(8) * mrSges(4,3)) + t130 * m(4) + (2 * pkin(7) * t146) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * t140); pkin(1) * t140 + mrSges(2,1); -pkin(7) * t140 + mrSges(2,2) - t146; Ifges(3,1) - Ifges(3,2) + (-t130 + t142) * m(4); Ifges(3,4); t149 * pkin(2) + Ifges(3,5); Ifges(3,6); t144 * m(4) + Ifges(3,3); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2); t155 * t151 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t147; t151 * t125 + t155 * t152 + Ifges(4,4); t139 * t122 - t137 * t124 + Ifges(4,5); t137 * t122 + t139 * t124 + Ifges(4,6); 0.2e1 * pkin(3) * t145 + t133 * t126 + t131 * t127 + Ifges(4,3) + Ifges(5,3) + 0.2e1 * t150; mrSges(4,1) + t145; t137 * mrSges(5,1) + t139 * mrSges(5,2) + mrSges(4,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t148; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
