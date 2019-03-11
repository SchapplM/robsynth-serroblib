% Return the minimum parameter vector for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRPR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR6_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR6_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t137 = (m(4) + m(5));
t143 = (pkin(3) ^ 2);
t153 = (t143 * m(5) + Ifges(4,2));
t152 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t133 = sin(pkin(11));
t135 = cos(pkin(11));
t151 = t133 * t135;
t150 = Ifges(6,4) * t151;
t149 = 2 * pkin(7) * mrSges(4,3) + t153;
t139 = pkin(9) ^ 2;
t142 = (pkin(5) ^ 2);
t123 = Ifges(6,2) + (t139 + t142) * m(7) + t152;
t125 = m(7) * t139 + Ifges(6,1) + t152;
t128 = t133 ^ 2;
t129 = t135 ^ 2;
t148 = t129 * t123 + t128 * t125 + Ifges(5,2);
t147 = pkin(8) * m(5) + mrSges(5,3);
t146 = -pkin(9) * m(7) - mrSges(7,3);
t145 = (2 * pkin(8) * mrSges(5,3)) + t148 + 0.2e1 * t150;
t127 = m(7) * pkin(5) + mrSges(6,1);
t144 = -t133 * mrSges(6,2) + t135 * t127;
t141 = pkin(7) ^ 2;
t140 = pkin(8) ^ 2;
t136 = cos(pkin(10));
t134 = sin(pkin(10));
t126 = t146 * pkin(5) + Ifges(6,5);
t1 = [Ifges(2,3) + t136 ^ 2 * (Ifges(3,2) + (pkin(2) ^ 2 + t141) * t137 + t149) + (0.2e1 * t136 * Ifges(3,4) + (t141 * t137 + Ifges(3,1) + t149) * t134) * t134; mrSges(2,1); mrSges(2,2); pkin(2) * t137 + mrSges(3,1); mrSges(3,2); pkin(7) * t137 + mrSges(3,3) + mrSges(4,3); m(3) + t137; (t140 * m(5)) + Ifges(4,1) + t145 - t153; t147 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + ((t140 + t143) * m(5)) + t145; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t147; t128 * t123 + t129 * t125 + Ifges(5,1) - t148 - 0.4e1 * t150; Ifges(5,4) + (t129 - t128) * Ifges(6,4) + (-t123 + t125) * t151; -t133 * Ifges(6,6) + t135 * t126 + Ifges(5,5); t135 * Ifges(6,6) + t133 * t126 + Ifges(5,6); (t142 * m(7)) + 0.2e1 * pkin(4) * t144 + Ifges(5,3) + Ifges(6,3); mrSges(5,1) + t144; t135 * mrSges(6,2) + t133 * t127 + mrSges(5,2); mrSges(6,3) - t146; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
