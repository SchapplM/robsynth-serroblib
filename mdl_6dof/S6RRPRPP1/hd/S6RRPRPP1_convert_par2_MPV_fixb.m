% Return the minimum parameter vector for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% MPV [25x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRPP1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t143 = pkin(8) ^ 2;
t141 = Ifges(6,4) - Ifges(7,5);
t134 = sin(pkin(10));
t136 = cos(pkin(10));
t154 = t134 * t136;
t149 = t141 * t154;
t129 = t134 ^ 2;
t131 = t136 ^ 2;
t139 = Ifges(6,2) + Ifges(7,3);
t142 = Ifges(6,1) + Ifges(7,1);
t151 = t129 * t142 + t131 * t139 + Ifges(5,2);
t147 = (2 * pkin(8) * mrSges(5,3)) + 0.2e1 * t149 + t151;
t121 = t143 * m(5) + Ifges(4,1) + t147;
t144 = pkin(3) ^ 2;
t127 = t144 * m(5) + Ifges(4,2);
t155 = t121 - t127;
t135 = sin(pkin(9));
t137 = cos(pkin(9));
t153 = t135 * t137;
t130 = t135 ^ 2;
t132 = t137 ^ 2;
t152 = t132 - t130;
t150 = pkin(8) * m(5) + mrSges(5,3);
t123 = pkin(3) * t150 + Ifges(4,4);
t148 = t123 * t153;
t146 = t136 * mrSges(6,1) - t134 * mrSges(6,2);
t126 = mrSges(4,2) - t150;
t128 = m(5) * pkin(3) + mrSges(4,1);
t145 = -t135 * t126 + t137 * t128;
t140 = Ifges(6,5) + Ifges(7,4);
t138 = Ifges(6,6) - Ifges(7,6);
t1 = [Ifges(2,3) + Ifges(3,2) + t130 * t121 + 0.2e1 * t148 + t132 * t127 + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); t155 * t152 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t148; t152 * t123 + t155 * t153 + Ifges(3,4); t137 * Ifges(4,5) - t135 * Ifges(4,6) + Ifges(3,5); t135 * Ifges(4,5) + t137 * Ifges(4,6) + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + (t143 + t144) * m(5) + 0.2e1 * t145 * pkin(2) + t147; mrSges(3,1) + t145; t137 * t126 + t135 * t128 + mrSges(3,2); mrSges(4,3); m(4) + m(5); t129 * t139 + t131 * t142 + Ifges(5,1) - 0.4e1 * t149 - t151; Ifges(5,4) + (t131 - t129) * t141 + (-t139 + t142) * t154; -t134 * t138 + t136 * t140 + Ifges(5,5); t134 * t140 + t136 * t138 + Ifges(5,6); 0.2e1 * pkin(4) * t146 + Ifges(7,2) + Ifges(5,3) + Ifges(6,3); mrSges(5,1) + t146; t134 * mrSges(6,1) + t136 * mrSges(6,2) + mrSges(5,2); mrSges(6,3); m(6); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
