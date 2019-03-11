% Return the minimum parameter vector for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% MPV [38x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRRR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t153 = -pkin(11) * m(7) - mrSges(7,3);
t127 = (mrSges(6,3) - t153);
t135 = (m(6) + m(7));
t152 = -pkin(10) * t135 - t127;
t138 = (pkin(10) ^ 2);
t142 = (pkin(4) ^ 2);
t150 = (Ifges(5,2) + (t138 + t142) * t135);
t137 = (pkin(11) ^ 2);
t141 = (pkin(5) ^ 2);
t149 = (Ifges(6,2) + (t137 + t141) * m(7));
t133 = (m(5) + t135);
t129 = (m(4) + t133);
t148 = -pkin(8) * t129 - mrSges(4,3);
t123 = (mrSges(5,3) - t152);
t147 = pkin(9) * t133 + t123;
t146 = (mrSges(3,3) - t148);
t145 = 2 * pkin(11) * mrSges(7,3) + 2 * pkin(9) * t123 + 2 * pkin(10) * t127 + Ifges(7,2) + t149 + t150;
t144 = (pkin(2) ^ 2);
t143 = (pkin(3) ^ 2);
t140 = (pkin(8) ^ 2);
t139 = pkin(9) ^ 2;
t132 = (t140 + t144);
t128 = (m(3) + t129);
t1 = [Ifges(2,3) + Ifges(3,2) + Ifges(4,2) + t143 * t133 + 2 * pkin(8) * mrSges(4,3) + t132 * t129 + 2 * pkin(7) * t146 + (pkin(1) ^ 2 + pkin(7) ^ 2) * t128; pkin(1) * t128 + mrSges(2,1); -pkin(7) * t128 + mrSges(2,2) - t146; Ifges(3,1) - Ifges(3,2) + (-t132 + t140) * t129; Ifges(3,4); t148 * pkin(2) + Ifges(3,5); Ifges(3,6); t144 * t129 + Ifges(3,3); pkin(2) * t129 + mrSges(3,1); mrSges(3,2); Ifges(4,1) - Ifges(4,2) + (t139 - t143) * t133 + t145; t147 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t139 + t143) * t133 + t145; pkin(3) * t133 + mrSges(4,1); mrSges(4,2) - t147; t138 * t135 + Ifges(5,1) - t150; Ifges(5,4); t152 * pkin(4) + Ifges(5,5); Ifges(5,6); t142 * t135 + Ifges(5,3); pkin(4) * t135 + mrSges(5,1); mrSges(5,2); m(7) * t137 + Ifges(6,1) - t149; Ifges(6,4); t153 * pkin(5) + Ifges(6,5); Ifges(6,6); t141 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
