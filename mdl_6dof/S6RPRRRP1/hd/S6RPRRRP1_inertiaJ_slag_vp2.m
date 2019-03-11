% Calculate joint inertia matrix for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:33
% EndTime: 2019-03-09 05:55:34
% DurationCPUTime: 1.03s
% Computational Cost: add. (1273->217), mult. (2377->297), div. (0->0), fcn. (2253->8), ass. (0->95)
t174 = Ifges(6,1) + Ifges(7,1);
t173 = Ifges(7,4) + Ifges(6,5);
t100 = sin(qJ(5));
t103 = cos(qJ(5));
t137 = t100 ^ 2 + t103 ^ 2;
t172 = mrSges(6,3) + mrSges(7,2);
t171 = Ifges(7,2) + Ifges(6,3);
t133 = Ifges(7,5) * t100;
t135 = Ifges(6,4) * t100;
t101 = sin(qJ(4));
t102 = sin(qJ(3));
t104 = cos(qJ(4));
t105 = cos(qJ(3));
t65 = t101 * t102 - t104 * t105;
t67 = t101 * t105 + t102 * t104;
t170 = (t174 * t103 + t133 - t135) * t67 + t173 * t65;
t132 = Ifges(7,5) * t103;
t134 = Ifges(6,4) * t103;
t169 = t174 * t100 - t132 + t134;
t168 = (Ifges(6,6) - Ifges(7,6)) * t103 + t173 * t100;
t85 = pkin(3) * t101 + pkin(9);
t167 = t137 * t85;
t166 = t137 * t67;
t157 = pkin(4) * t65;
t99 = cos(pkin(10));
t83 = -pkin(1) * t99 - pkin(2);
t68 = -pkin(3) * t105 + t83;
t22 = -pkin(9) * t67 + t157 + t68;
t98 = sin(pkin(10));
t82 = pkin(1) * t98 + pkin(7);
t155 = pkin(8) + t82;
t127 = t155 * t102;
t55 = t155 * t105;
t26 = -t101 * t127 + t104 * t55;
t6 = -t100 * t26 + t103 * t22;
t7 = t100 * t22 + t103 * t26;
t122 = -t100 * t6 + t103 * t7;
t165 = t105 ^ 2;
t164 = t172 * t137;
t24 = t101 * t55 + t104 * t127;
t163 = t24 ^ 2;
t162 = t65 ^ 2;
t161 = 0.2e1 * t24;
t160 = 0.2e1 * t68;
t72 = -mrSges(7,1) * t103 - mrSges(7,3) * t100;
t159 = 0.2e1 * t72;
t154 = m(7) * t100;
t153 = Ifges(6,6) * t65;
t152 = pkin(3) * t104;
t3 = qJ(6) * t65 + t7;
t150 = t103 * t3;
t148 = t24 * t65;
t143 = t167 * t67;
t142 = t166 * pkin(9);
t141 = t167 * pkin(9);
t140 = t137 * t85 ^ 2;
t138 = t137 * pkin(9) ^ 2;
t136 = t102 ^ 2 + t165;
t88 = t100 * mrSges(7,2);
t131 = t100 * t67;
t130 = t103 * t67;
t129 = qJ(6) * t103;
t31 = -t65 * mrSges(7,1) + mrSges(7,2) * t130;
t124 = Ifges(7,6) * t131 + t173 * t130 + t171 * t65;
t4 = -pkin(5) * t65 - t6;
t123 = t100 * t4 + t150;
t121 = -t105 * mrSges(4,1) + t102 * mrSges(4,2);
t120 = t100 * mrSges(6,1) + t103 * mrSges(6,2);
t119 = t100 * mrSges(7,1) - t103 * mrSges(7,3);
t118 = pkin(5) * t100 - t129;
t69 = -pkin(5) * t103 - qJ(6) * t100 - pkin(4);
t74 = -Ifges(7,3) * t103 + t133;
t75 = Ifges(6,2) * t103 + t135;
t117 = Ifges(5,3) + (-t74 + t75) * t103 + t169 * t100;
t116 = (mrSges(5,1) * t104 - mrSges(5,2) * t101) * pkin(3);
t29 = -mrSges(6,2) * t65 - mrSges(6,3) * t131;
t30 = mrSges(6,1) * t65 - mrSges(6,3) * t130;
t32 = -mrSges(7,2) * t131 + mrSges(7,3) * t65;
t114 = (t29 + t32) * t103 + (-t30 + t31) * t100;
t56 = t65 * mrSges(5,1);
t73 = -mrSges(6,1) * t103 + mrSges(6,2) * t100;
t113 = -t67 * mrSges(5,2) - t56 + (t72 + t73) * t65 + t172 * t166;
t112 = mrSges(7,2) * t129 - pkin(5) * t88 + t168;
t111 = 0.2e1 * t164;
t110 = -m(7) * t118 - t119 - t120;
t14 = Ifges(7,6) * t65 + (Ifges(7,3) * t100 + t132) * t67;
t15 = t153 + (-Ifges(6,2) * t100 + t134) * t67;
t9 = t118 * t67 + t24;
t109 = -t26 * mrSges(5,2) + mrSges(7,2) * t150 + Ifges(5,5) * t67 + t4 * t88 + t9 * t72 + (t73 - mrSges(5,1)) * t24 + t170 * t100 / 0.2e1 + (t74 / 0.2e1 - t75 / 0.2e1) * t131 + t169 * t130 / 0.2e1 + (-t14 / 0.2e1 + t15 / 0.2e1) * t103 + t122 * mrSges(6,3) + (-Ifges(5,6) + t168 / 0.2e1) * t65;
t86 = -pkin(4) - t152;
t61 = t67 ^ 2;
t54 = t69 - t152;
t28 = t120 * t67;
t27 = t119 * t67;
t1 = [0.2e1 * t83 * t121 + Ifges(4,2) * t165 + t56 * t160 + 0.2e1 * t9 * t27 + t28 * t161 + 0.2e1 * t7 * t29 + 0.2e1 * t6 * t30 + 0.2e1 * t4 * t31 + 0.2e1 * t3 * t32 + Ifges(2,3) + Ifges(3,3) + (-0.2e1 * mrSges(5,3) * t26 + Ifges(5,2) * t65 + t124) * t65 + (mrSges(5,2) * t160 + mrSges(5,3) * t161 + Ifges(5,1) * t67 - 0.2e1 * Ifges(5,4) * t65 + t170 * t103 + (t14 - t15 - t153) * t100) * t67 + m(4) * (t136 * t82 ^ 2 + t83 ^ 2) + m(5) * (t26 ^ 2 + t68 ^ 2 + t163) + m(6) * (t6 ^ 2 + t7 ^ 2 + t163) + m(7) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) + m(3) * (t98 ^ 2 + t99 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t102 + 0.2e1 * Ifges(4,4) * t105) * t102 + 0.2e1 * (mrSges(3,1) * t99 - mrSges(3,2) * t98) * pkin(1) + 0.2e1 * t136 * t82 * mrSges(4,3); (t27 + t28) * t65 + t114 * t67 + m(6) * (t122 * t67 + t148) + m(7) * (t123 * t67 + t65 * t9) + m(5) * (t26 * t67 + t148); m(3) + m(5) * (t61 + t162) + m(4) * t136 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t137 * t61 + t162); t114 * t85 + (-mrSges(4,2) * t82 + Ifges(4,6)) * t105 + (-mrSges(4,1) * t82 + Ifges(4,5)) * t102 + t109 + t86 * t28 + t54 * t27 + m(6) * (t122 * t85 + t24 * t86) + m(7) * (t123 * t85 + t54 * t9) + (m(5) * (t101 * t26 - t104 * t24) + (-t101 * t65 - t104 * t67) * mrSges(5,3)) * pkin(3); m(6) * (t65 * t86 + t143) + m(7) * (t54 * t65 + t143) + m(5) * (t101 * t67 - t104 * t65) * pkin(3) + t113 - t121; t54 * t159 + 0.2e1 * t86 * t73 + Ifges(4,3) + 0.2e1 * t116 + m(7) * (t54 ^ 2 + t140) + m(6) * (t86 ^ 2 + t140) + m(5) * (t101 ^ 2 + t104 ^ 2) * pkin(3) ^ 2 + t111 * t85 + t117; t114 * pkin(9) + t109 + t69 * t27 - pkin(4) * t28 + m(6) * (-pkin(4) * t24 + pkin(9) * t122) + m(7) * (pkin(9) * t123 + t69 * t9); m(6) * (t142 - t157) + m(7) * (t65 * t69 + t142) + t113; (t86 - pkin(4)) * t73 + (t54 + t69) * t72 + t116 + m(7) * (t54 * t69 + t141) + m(6) * (-pkin(4) * t86 + t141) + t117 + (pkin(9) + t85) * t164; -0.2e1 * pkin(4) * t73 + t69 * t159 + m(7) * (t69 ^ 2 + t138) + m(6) * (pkin(4) ^ 2 + t138) + t111 * pkin(9) + t117; -Ifges(6,6) * t131 + m(7) * (-pkin(5) * t4 + qJ(6) * t3) - t4 * mrSges(7,1) - t7 * mrSges(6,2) - pkin(5) * t31 + t3 * mrSges(7,3) + qJ(6) * t32 + t6 * mrSges(6,1) + t124; t110 * t67; t110 * t85 + t112; pkin(9) * t110 + t112; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t171; m(7) * t4 + t31; m(7) * t131; t154 * t85 + t88; pkin(9) * t154 + t88; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
