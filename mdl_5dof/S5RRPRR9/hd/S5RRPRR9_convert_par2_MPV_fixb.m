% Return the minimum parameter vector for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [26x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPRR9_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR9_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t113 = (m(5) + m(6));
t117 = (pkin(3) ^ 2);
t104 = (t117 * t113 + Ifges(4,2));
t115 = (pkin(7) ^ 2);
t122 = -pkin(8) * m(6) - mrSges(6,3);
t106 = (mrSges(5,3) - t122);
t114 = (pkin(8) ^ 2);
t116 = (pkin(4) ^ 2);
t125 = (Ifges(5,2) + (t114 + t116) * m(6));
t119 = 2 * pkin(8) * mrSges(6,3) + 2 * pkin(7) * t106 + Ifges(6,2) + t125;
t99 = t115 * t113 + Ifges(4,1) + t119;
t126 = -t104 + t99;
t111 = sin(pkin(9));
t112 = cos(pkin(9));
t124 = t111 * t112;
t108 = t111 ^ 2;
t109 = t112 ^ 2;
t123 = t109 - t108;
t120 = pkin(7) * t113 + t106;
t100 = t120 * pkin(3) + Ifges(4,4);
t121 = t100 * t124;
t101 = mrSges(4,2) - t120;
t103 = pkin(3) * t113 + mrSges(4,1);
t118 = -t111 * t101 + t112 * t103;
t1 = [Ifges(2,3) + Ifges(3,2) + t108 * t99 + 0.2e1 * t121 + t109 * t104 + (2 * pkin(6) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(6) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(6) * m(3) + mrSges(2,2) - mrSges(3,3); t126 * t123 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t121; t123 * t100 + t126 * t124 + Ifges(3,4); t112 * Ifges(4,5) - t111 * Ifges(4,6) + Ifges(3,5); t111 * Ifges(4,5) + t112 * Ifges(4,6) + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + ((t115 + t117) * t113) + 0.2e1 * t118 * pkin(2) + t119; mrSges(3,1) + t118; t112 * t101 + t111 * t103 + mrSges(3,2); mrSges(4,3); m(4) + t113; m(6) * t114 + Ifges(5,1) - t125; Ifges(5,4); t122 * pkin(4) + Ifges(5,5); Ifges(5,6); t116 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
