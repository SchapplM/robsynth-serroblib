% Return the minimum parameter vector for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% MPV [24x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRPP1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t127 = (m(4) + m(5));
t118 = sin(pkin(10));
t120 = cos(pkin(10));
t136 = t118 * t120;
t115 = t118 ^ 2;
t116 = t120 ^ 2;
t123 = Ifges(6,2) + Ifges(7,3);
t126 = Ifges(6,1) + Ifges(7,1);
t135 = t115 * t126 + t116 * t123 + Ifges(5,2);
t134 = pkin(8) * m(5) + mrSges(5,3);
t125 = Ifges(6,4) - Ifges(7,5);
t133 = t125 * t136;
t132 = (2 * pkin(8) * mrSges(5,3)) + 0.2e1 * t133 + t135;
t131 = t120 * mrSges(6,1) - t118 * mrSges(6,2);
t111 = -pkin(7) * t127 + mrSges(3,2) - mrSges(4,3);
t114 = pkin(2) * t127 + mrSges(3,1);
t119 = sin(pkin(9));
t121 = cos(pkin(9));
t130 = -t119 * t111 + t121 * t114;
t129 = pkin(3) ^ 2;
t128 = pkin(8) ^ 2;
t124 = Ifges(6,5) + Ifges(7,4);
t122 = Ifges(6,6) - Ifges(7,6);
t1 = [Ifges(2,3) + Ifges(3,3) + t129 * m(5) + Ifges(4,2) + (2 * pkin(7) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(7) ^ 2) * t127) + 0.2e1 * t130 * pkin(1); mrSges(2,1) + t130; t121 * t111 + t119 * t114 + mrSges(2,2); m(3) + t127; Ifges(4,1) - Ifges(4,2) + (t128 - t129) * m(5) + t132; t134 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t128 + t129) * m(5) + t132; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t134; t115 * t123 + t116 * t126 + Ifges(5,1) - 0.4e1 * t133 - t135; Ifges(5,4) + (t116 - t115) * t125 + (-t123 + t126) * t136; -t118 * t122 + t120 * t124 + Ifges(5,5); t118 * t124 + t120 * t122 + Ifges(5,6); 0.2e1 * pkin(4) * t131 + Ifges(7,2) + Ifges(5,3) + Ifges(6,3); mrSges(5,1) + t131; t118 * mrSges(6,1) + t120 * mrSges(6,2) + mrSges(5,2); mrSges(6,3); m(6); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
