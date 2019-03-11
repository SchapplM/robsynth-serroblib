% Return the minimum parameter vector for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPRR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t140 = -pkin(9) * m(7) - mrSges(7,3);
t115 = (mrSges(6,3) - t140);
t126 = (m(6) + m(7));
t139 = -pkin(8) * t126 - t115;
t128 = (pkin(9) ^ 2);
t130 = (pkin(5) ^ 2);
t137 = (Ifges(6,2) + (t128 + t130) * m(7));
t122 = sin(pkin(11));
t124 = cos(pkin(11));
t136 = t122 * t124;
t135 = Ifges(5,4) * t136;
t131 = (pkin(4) ^ 2);
t134 = t131 * t126 + Ifges(4,2) + Ifges(5,3);
t133 = 2 * pkin(9) * mrSges(7,3) + 2 * pkin(8) * t115 + Ifges(7,2) + t137;
t116 = -pkin(7) * m(4) + mrSges(3,2) - mrSges(4,3);
t117 = m(4) * pkin(2) + mrSges(3,1);
t123 = sin(pkin(10));
t125 = cos(pkin(10));
t132 = -t123 * t116 + t125 * t117;
t129 = pkin(8) ^ 2;
t120 = t124 ^ 2;
t119 = t122 ^ 2;
t112 = t139 * pkin(4) + Ifges(5,5);
t111 = t129 * t126 + Ifges(5,1) + t133;
t110 = Ifges(5,2) + (t129 + t131) * t126 + t133;
t1 = [Ifges(2,3) + Ifges(3,3) + (2 * pkin(7) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(7) ^ 2) * m(4)) + 0.2e1 * t132 * pkin(1) + t134; mrSges(2,1) + t132; t125 * t116 + t123 * t117 + mrSges(2,2); m(3) + m(4); t119 * t110 + t120 * t111 + Ifges(4,1) - t134 - 0.2e1 * t135; t122 * Ifges(5,6) - t124 * t112 + Ifges(4,4); Ifges(4,5) + (t120 - t119) * Ifges(5,4) + (-t110 + t111) * t136; -t124 * Ifges(5,6) - t122 * t112 + Ifges(4,6); t120 * t110 + t119 * t111 + Ifges(4,3) + 0.2e1 * t135; mrSges(4,1); mrSges(4,2); pkin(4) * t126 + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t139; m(5) + t126; m(7) * t128 + Ifges(6,1) - t137; Ifges(6,4); t140 * pkin(5) + Ifges(6,5); Ifges(6,6); t130 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
