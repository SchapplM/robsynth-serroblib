% Return the minimum parameter vector for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% MPV [26x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRPP6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP6_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t120 = m(4) + m(5);
t129 = (pkin(7) * t120);
t113 = sin(pkin(9));
t114 = cos(pkin(9));
t128 = t113 * t114;
t110 = t113 ^ 2;
t111 = t114 ^ 2;
t116 = Ifges(6,2) + Ifges(7,3);
t119 = Ifges(6,1) + Ifges(7,1);
t127 = t110 * t119 + t111 * t116 + Ifges(5,2);
t126 = pkin(8) * m(5) + mrSges(5,3);
t118 = Ifges(6,4) - Ifges(7,5);
t125 = t118 * t128;
t124 = (2 * pkin(8) * mrSges(5,3)) + 0.2e1 * t125 + t127;
t123 = t114 * mrSges(6,1) - t113 * mrSges(6,2);
t122 = (pkin(3) ^ 2);
t121 = pkin(8) ^ 2;
t117 = Ifges(6,5) + Ifges(7,4);
t115 = Ifges(6,6) - Ifges(7,6);
t1 = [t122 * m(5) + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + (2 * mrSges(4,3) + t129) * pkin(7); mrSges(2,1); mrSges(2,2); mrSges(3,2) - mrSges(4,3) - t129; mrSges(3,3); m(3) + t120; Ifges(4,1) - Ifges(4,2) + ((t121 - t122) * m(5)) + t124; t126 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + ((t121 + t122) * m(5)) + t124; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t126; t110 * t116 + t111 * t119 + Ifges(5,1) - 0.4e1 * t125 - t127; Ifges(5,4) + (t111 - t110) * t118 + (-t116 + t119) * t128; -t113 * t115 + t114 * t117 + Ifges(5,5); t113 * t117 + t114 * t115 + Ifges(5,6); 0.2e1 * pkin(4) * t123 + Ifges(7,2) + Ifges(5,3) + Ifges(6,3); mrSges(5,1) + t123; t113 * mrSges(6,1) + t114 * mrSges(6,2) + mrSges(5,2); mrSges(6,3); m(6); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
