% Return the minimum parameter vector for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPPR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR2_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t129 = (pkin(9) ^ 2);
t131 = (pkin(5) ^ 2);
t140 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t117 = Ifges(5,2) + Ifges(6,3) + (t129 + t131) * m(7) + t140;
t119 = t131 * m(7) + Ifges(5,1) + Ifges(6,2);
t141 = -t117 + t119;
t124 = sin(pkin(10));
t125 = cos(pkin(10));
t139 = t124 * t125;
t121 = t124 ^ 2;
t122 = t125 ^ 2;
t138 = t122 - t121;
t137 = -pkin(8) * m(4) - mrSges(4,3);
t136 = -pkin(9) * m(7) - mrSges(7,3);
t127 = Ifges(5,4) + Ifges(6,6);
t135 = t127 * t139;
t134 = (mrSges(3,3) - t137);
t133 = t125 * mrSges(5,1) - t124 * mrSges(5,2);
t132 = pkin(2) ^ 2;
t130 = pkin(8) ^ 2;
t128 = (m(3) + m(4));
t126 = Ifges(5,6) - Ifges(6,5);
t120 = t130 + t132;
t118 = t136 * pkin(5) - Ifges(6,4) + Ifges(5,5);
t1 = [Ifges(2,3) + Ifges(3,2) + Ifges(4,2) + t121 * t119 + 0.2e1 * t135 + t122 * t117 + (2 * pkin(8) * mrSges(4,3)) + t120 * m(4) + (2 * pkin(7) * t134) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * t128); pkin(1) * t128 + mrSges(2,1); -pkin(7) * t128 + mrSges(2,2) - t134; Ifges(3,1) - Ifges(3,2) + (-t120 + t130) * m(4); Ifges(3,4); t137 * pkin(2) + Ifges(3,5); Ifges(3,6); t132 * m(4) + Ifges(3,3); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2); t141 * t138 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t135; t138 * t127 + t141 * t139 + Ifges(4,4); t125 * t118 - t124 * t126 + Ifges(4,5); t124 * t118 + t125 * t126 + Ifges(4,6); (m(7) * t129) + 0.2e1 * pkin(3) * t133 + Ifges(6,1) + Ifges(4,3) + Ifges(5,3) + t140; mrSges(4,1) + t133; t124 * mrSges(5,1) + t125 * mrSges(5,2) + mrSges(4,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) + t136; mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
