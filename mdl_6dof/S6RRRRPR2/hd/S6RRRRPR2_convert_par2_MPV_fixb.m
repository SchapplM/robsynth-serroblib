% Return the minimum parameter vector for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% MPV [35x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRPR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR2_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR2_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t150 = -pkin(9) * m(5) - mrSges(5,3);
t125 = (mrSges(4,3) - t150);
t134 = (m(4) + m(5));
t149 = -pkin(8) * t134 - t125;
t147 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t132 = sin(pkin(11));
t133 = cos(pkin(11));
t146 = t132 * t133;
t145 = Ifges(6,4) * t146;
t144 = -pkin(10) * m(7) - mrSges(7,3);
t139 = (pkin(5) ^ 2);
t143 = t139 * m(7) + Ifges(5,2) + Ifges(6,3);
t142 = (mrSges(3,3) - t149);
t141 = (pkin(2) ^ 2);
t140 = (pkin(3) ^ 2);
t138 = (pkin(8) ^ 2);
t137 = (pkin(9) ^ 2);
t136 = pkin(10) ^ 2;
t130 = (m(3) + t134);
t129 = t133 ^ 2;
t128 = t132 ^ 2;
t127 = (t138 + t141);
t126 = (t137 + t140);
t124 = t144 * pkin(5) + Ifges(6,5);
t123 = m(7) * t136 + Ifges(6,1) + t147;
t122 = Ifges(6,2) + (t136 + t139) * m(7) + t147;
t1 = [Ifges(2,3) + Ifges(3,2) + Ifges(4,2) + 2 * pkin(9) * mrSges(5,3) + t126 * m(5) + 2 * pkin(8) * t125 + t127 * t134 + 2 * pkin(7) * t142 + (pkin(1) ^ 2 + pkin(7) ^ 2) * t130 + t143; pkin(1) * t130 + mrSges(2,1); -pkin(7) * t130 + mrSges(2,2) - t142; Ifges(3,1) - Ifges(3,2) + (-t127 + t138) * t134; Ifges(3,4); t149 * pkin(2) + Ifges(3,5); Ifges(3,6); t141 * t134 + Ifges(3,3); pkin(2) * t134 + mrSges(3,1); mrSges(3,2); Ifges(4,1) - Ifges(4,2) + (-t126 + t137) * m(5); Ifges(4,4); t150 * pkin(3) + Ifges(4,5); Ifges(4,6); t140 * m(5) + Ifges(4,3); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2); t128 * t122 + t129 * t123 + Ifges(5,1) - t143 - 0.2e1 * t145; t132 * Ifges(6,6) - t133 * t124 + Ifges(5,4); Ifges(5,5) + (t129 - t128) * Ifges(6,4) + (-t122 + t123) * t146; -t133 * Ifges(6,6) - t132 * t124 + Ifges(5,6); t129 * t122 + t128 * t123 + Ifges(5,3) + 0.2e1 * t145; mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t144; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
