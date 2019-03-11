% Return the minimum parameter vector for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPPRPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t129 = 2 * pkin(8) * mrSges(7,3) + Ifges(7,2);
t113 = sin(pkin(11));
t116 = cos(pkin(11));
t128 = t113 * t116;
t122 = (pkin(5) ^ 2);
t127 = (t122 * m(7) + Ifges(5,2) + Ifges(6,3));
t126 = Ifges(6,4) * t128;
t125 = -pkin(8) * m(7) - mrSges(7,3);
t124 = 2 * pkin(7) * mrSges(5,3) + t127;
t115 = sin(pkin(9));
t118 = cos(pkin(9));
t123 = t118 * mrSges(3,1) - t115 * mrSges(3,2);
t121 = pkin(7) ^ 2;
t120 = pkin(8) ^ 2;
t117 = cos(pkin(10));
t114 = sin(pkin(10));
t110 = t116 ^ 2;
t109 = t113 ^ 2;
t108 = t125 * pkin(5) + Ifges(6,5);
t107 = m(7) * t120 + Ifges(6,1) + t129;
t106 = Ifges(6,2) + (t120 + t122) * m(7) + t129;
t1 = [Ifges(2,3) + Ifges(3,3) + t117 ^ 2 * (Ifges(4,2) + (pkin(3) ^ 2 + t121) * m(5) + t124) + (0.2e1 * t117 * Ifges(4,4) + (t121 * m(5) + Ifges(4,1) + t124) * t114) * t114 + 0.2e1 * t123 * pkin(1); mrSges(2,1) + t123; t115 * mrSges(3,1) + t118 * mrSges(3,2) + mrSges(2,2); m(3); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2); pkin(7) * m(5) + mrSges(4,3) + mrSges(5,3); m(4) + m(5); t109 * t106 + t110 * t107 + Ifges(5,1) - 0.2e1 * t126 - t127; t113 * Ifges(6,6) - t116 * t108 + Ifges(5,4); Ifges(5,5) + (t110 - t109) * Ifges(6,4) + (-t106 + t107) * t128; -t116 * Ifges(6,6) - t113 * t108 + Ifges(5,6); t110 * t106 + t109 * t107 + Ifges(5,3) + 0.2e1 * t126; mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t125; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
