% Return the minimum parameter vector for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRPR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR7_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR7_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t110 = m(4) + m(5);
t124 = (pkin(7) * t110);
t113 = (pkin(5) ^ 2);
t102 = (t113 * m(7) + Ifges(6,2));
t111 = (pkin(9) ^ 2);
t122 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t99 = m(7) * t111 + Ifges(6,1) + t122;
t123 = -t102 + t99;
t108 = sin(pkin(10));
t109 = cos(pkin(10));
t121 = t108 * t109;
t105 = t108 ^ 2;
t106 = t109 ^ 2;
t120 = t106 - t105;
t119 = -pkin(8) * m(5) - mrSges(5,3);
t118 = pkin(9) * m(7) + mrSges(7,3);
t100 = t118 * pkin(5) + Ifges(6,4);
t117 = t100 * t121;
t116 = (mrSges(4,3) - t119);
t101 = mrSges(6,2) - t118;
t103 = m(7) * pkin(5) + mrSges(6,1);
t115 = -t108 * t101 + t109 * t103;
t114 = pkin(3) ^ 2;
t112 = pkin(8) ^ 2;
t104 = t112 + t114;
t1 = [t104 * m(5) + (2 * pkin(8) * mrSges(5,3)) + t106 * t102 + t105 * t99 + Ifges(3,1) + Ifges(4,2) + Ifges(5,2) + Ifges(2,3) + 0.2e1 * t117 + ((2 * t116 + t124) * pkin(7)); mrSges(2,1); mrSges(2,2); mrSges(3,2) - t116 - t124; mrSges(3,3); m(3) + t110; Ifges(4,1) - Ifges(4,2) + (-t104 + t112) * m(5); Ifges(4,4); t119 * pkin(3) + Ifges(4,5); Ifges(4,6); t114 * m(5) + Ifges(4,3); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2); t123 * t120 + Ifges(5,1) - Ifges(5,2) - 0.4e1 * t117; t120 * t100 + t123 * t121 + Ifges(5,4); t109 * Ifges(6,5) - t108 * Ifges(6,6) + Ifges(5,5); t108 * Ifges(6,5) + t109 * Ifges(6,6) + Ifges(5,6); Ifges(5,3) + Ifges(6,3) + ((t111 + t113) * m(7)) + 0.2e1 * t115 * pkin(4) + t122; mrSges(5,1) + t115; t109 * t101 + t108 * t103 + mrSges(5,2); mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
