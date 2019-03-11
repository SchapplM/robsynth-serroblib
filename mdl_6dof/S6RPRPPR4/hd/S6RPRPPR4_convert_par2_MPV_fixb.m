% Return the minimum parameter vector for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% MPV [29x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPPR4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR4_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t149 = 2 * pkin(8) * mrSges(7,3) + Ifges(7,2);
t134 = sin(pkin(10));
t136 = cos(pkin(10));
t148 = t134 * t136;
t147 = pkin(8) * m(7) + mrSges(7,3);
t139 = Ifges(5,4) - Ifges(6,5);
t146 = t139 * t148;
t143 = (pkin(5) ^ 2);
t145 = t143 * m(7) + Ifges(4,2) + Ifges(6,2) + Ifges(5,3);
t144 = 2 * pkin(7) * mrSges(4,3) + t145;
t142 = pkin(7) ^ 2;
t141 = pkin(8) ^ 2;
t138 = Ifges(5,6) - Ifges(6,6);
t137 = cos(pkin(9));
t135 = sin(pkin(9));
t131 = t136 ^ 2;
t130 = t134 ^ 2;
t129 = t147 * pkin(5) + Ifges(6,4) + Ifges(5,5);
t128 = m(7) * t141 + Ifges(5,1) + Ifges(6,1) + t149;
t127 = Ifges(5,2) + Ifges(6,3) + (t141 + t143) * m(7) + t149;
t1 = [Ifges(2,3) + t137 ^ 2 * (Ifges(3,2) + (pkin(2) ^ 2 + t142) * m(4) + t144) + (0.2e1 * t137 * Ifges(3,4) + (t142 * m(4) + Ifges(3,1) + t144) * t135) * t135; mrSges(2,1); mrSges(2,2); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2); pkin(7) * m(4) + mrSges(3,3) + mrSges(4,3); m(3) + m(4); t130 * t127 + t131 * t128 + Ifges(4,1) - t145 - 0.2e1 * t146; -t136 * t129 + t134 * t138 + Ifges(4,4); Ifges(4,5) + (t131 - t130) * t139 + (-t127 + t128) * t148; -t134 * t129 - t136 * t138 + Ifges(4,6); t131 * t127 + t130 * t128 + Ifges(4,3) + 0.2e1 * t146; mrSges(4,1); mrSges(4,2); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t147; mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
