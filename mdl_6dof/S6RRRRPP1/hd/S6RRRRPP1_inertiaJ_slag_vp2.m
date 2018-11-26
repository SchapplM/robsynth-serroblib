% Calculate joint inertia matrix for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:04
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:04:18
% EndTime: 2018-11-23 18:04:19
% DurationCPUTime: 1.50s
% Computational Cost: add. (2492->301), mult. (4757->413), div. (0->0), fcn. (5101->8), ass. (0->113)
t204 = Ifges(6,1) + Ifges(7,1);
t203 = -Ifges(6,4) + Ifges(7,5);
t202 = Ifges(7,4) + Ifges(6,5);
t201 = -Ifges(6,6) + Ifges(7,6);
t195 = (mrSges(6,3) + mrSges(7,2));
t200 = 2 * t195;
t138 = sin(qJ(3));
t139 = sin(qJ(2));
t141 = cos(qJ(3));
t142 = cos(qJ(2));
t110 = t138 * t139 - t141 * t142;
t135 = sin(pkin(10));
t136 = cos(pkin(10));
t137 = sin(qJ(4));
t140 = cos(qJ(4));
t106 = t135 * t140 + t136 * t137;
t111 = t138 * t142 + t139 * t141;
t55 = t106 * t111;
t105 = t135 * t137 - t136 * t140;
t56 = t105 * t111;
t199 = t202 * t110 + t203 * t55 - t204 * t56;
t198 = t203 * t105 + t204 * t106;
t127 = -pkin(2) * t142 - pkin(1);
t66 = pkin(3) * t110 - pkin(9) * t111 + t127;
t188 = -pkin(8) - pkin(7);
t118 = t188 * t142;
t164 = t188 * t139;
t85 = -t141 * t118 + t138 * t164;
t29 = -t137 * t85 + t140 * t66;
t30 = t137 * t66 + t140 * t85;
t158 = -t137 * t29 + t140 * t30;
t197 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t196 = Ifges(5,5) * t137 + Ifges(5,6) * t140 + t201 * t105 + t202 * t106;
t83 = -t118 * t138 - t141 * t164;
t194 = t83 ^ 2;
t69 = t105 * mrSges(7,1) - t106 * mrSges(7,3);
t193 = 0.2e1 * t69;
t70 = t105 * mrSges(6,1) + t106 * mrSges(6,2);
t192 = 0.2e1 * t70;
t191 = 0.2e1 * t83;
t190 = 0.2e1 * t127;
t185 = pkin(2) * t141;
t184 = pkin(4) * t135;
t128 = t140 * qJ(5);
t14 = pkin(4) * t110 - t111 * t128 + t29;
t172 = t111 * t137;
t20 = -qJ(5) * t172 + t30;
t6 = t135 * t14 + t136 * t20;
t32 = -mrSges(6,2) * t110 - mrSges(6,3) * t55;
t35 = -mrSges(7,2) * t55 + mrSges(7,3) * t110;
t183 = t32 + t35;
t33 = mrSges(6,1) * t110 + mrSges(6,3) * t56;
t34 = -t110 * mrSges(7,1) - t56 * mrSges(7,2);
t182 = t34 - t33;
t179 = mrSges(6,3) * t105;
t178 = Ifges(5,4) * t137;
t177 = Ifges(5,4) * t140;
t176 = t105 * mrSges(7,2);
t93 = t106 * mrSges(7,2);
t175 = t106 * mrSges(6,3);
t171 = t111 * t140;
t169 = t137 ^ 2 + t140 ^ 2;
t168 = t139 ^ 2 + t142 ^ 2;
t124 = pkin(2) * t138 + pkin(9);
t161 = (-qJ(5) - t124) * t137;
t99 = t124 * t140 + t128;
t60 = t135 * t99 - t136 * t161;
t62 = t135 * t161 + t136 * t99;
t167 = t60 ^ 2 + t62 ^ 2;
t115 = pkin(9) * t140 + t128;
t162 = (-qJ(5) - pkin(9)) * t137;
t78 = t115 * t135 - t136 * t162;
t80 = t136 * t115 + t135 * t162;
t166 = t78 ^ 2 + t80 ^ 2;
t165 = t136 * t175;
t126 = -pkin(4) * t140 - pkin(3);
t22 = t55 * mrSges(6,1) - t56 * mrSges(6,2);
t21 = t55 * mrSges(7,1) + t56 * mrSges(7,3);
t163 = t60 * t78 + t80 * t62;
t160 = t169 * t124;
t46 = pkin(4) * t172 + t83;
t159 = mrSges(5,1) * t137 + mrSges(5,2) * t140;
t5 = -t135 * t20 + t136 * t14;
t67 = -mrSges(5,2) * t110 - mrSges(5,3) * t172;
t68 = mrSges(5,1) * t110 - mrSges(5,3) * t171;
t157 = -t137 * t68 + t140 * t67;
t156 = 0.2e1 * mrSges(5,3) * t169;
t155 = (t141 * mrSges(4,1) - t138 * mrSges(4,2)) * pkin(2);
t116 = Ifges(5,2) * t140 + t178;
t117 = Ifges(5,1) * t137 + t177;
t71 = Ifges(7,5) * t106 + Ifges(7,3) * t105;
t72 = Ifges(6,4) * t106 - Ifges(6,2) * t105;
t154 = t140 * t116 + t137 * t117 + Ifges(4,3) + t198 * t106 + (t71 - t72) * t105;
t153 = t69 + t70;
t150 = Ifges(5,5) * t171 + t197 * t110 + t201 * t55 - t202 * t56;
t63 = pkin(5) * t105 - qJ(6) * t106 + t126;
t119 = qJ(6) + t184;
t122 = -pkin(4) * t136 - pkin(5);
t148 = -t119 * t176 + t122 * t93 - t179 * t184 + t196;
t114 = -mrSges(5,1) * t140 + mrSges(5,2) * t137;
t15 = -Ifges(7,5) * t56 + Ifges(7,6) * t110 + Ifges(7,3) * t55;
t16 = -Ifges(6,4) * t56 - Ifges(6,2) * t55 + Ifges(6,6) * t110;
t3 = qJ(6) * t110 + t6;
t38 = Ifges(5,6) * t110 + (-Ifges(5,2) * t137 + t177) * t111;
t39 = Ifges(5,5) * t110 + (Ifges(5,1) * t140 - t178) * t111;
t4 = -pkin(5) * t110 - t5;
t8 = pkin(5) * t55 + qJ(6) * t56 + t46;
t147 = t140 * t38 / 0.2e1 + t137 * t39 / 0.2e1 + Ifges(4,5) * t111 - t85 * mrSges(4,2) + t8 * t69 + t46 * t70 - t5 * t175 - t3 * t176 - t116 * t172 / 0.2e1 + t117 * t171 / 0.2e1 + t4 * t93 - t6 * t179 + (t114 - mrSges(4,1)) * t83 + (t71 / 0.2e1 - t72 / 0.2e1) * t55 - t198 * t56 / 0.2e1 + t199 * t106 / 0.2e1 + (t15 / 0.2e1 - t16 / 0.2e1) * t105 + t158 * mrSges(5,3) + (-Ifges(4,6) + t196 / 0.2e1) * t110;
t125 = -pkin(3) - t185;
t113 = t126 - t185;
t65 = t159 * t111;
t54 = t63 - t185;
t1 = [-0.2e1 * pkin(1) * (-mrSges(3,1) * t142 + mrSges(3,2) * t139) + t142 * (Ifges(3,4) * t139 + Ifges(3,2) * t142) + t139 * (Ifges(3,1) * t139 + Ifges(3,4) * t142) + t65 * t191 + 0.2e1 * t30 * t67 + 0.2e1 * t29 * t68 + 0.2e1 * t6 * t32 + 0.2e1 * t5 * t33 + 0.2e1 * t4 * t34 + 0.2e1 * t3 * t35 + 0.2e1 * t46 * t22 + 0.2e1 * t8 * t21 + Ifges(2,3) - t199 * t56 + (t15 - t16) * t55 + 0.2e1 * t168 * pkin(7) * mrSges(3,3) + (mrSges(4,2) * t190 + mrSges(4,3) * t191 + Ifges(4,1) * t111 - t137 * t38 + t140 * t39) * t111 + (mrSges(4,1) * t190 - 0.2e1 * t85 * mrSges(4,3) + Ifges(4,2) * t110 + (-Ifges(5,6) * t137 - (2 * Ifges(4,4))) * t111 + t150) * t110 + m(3) * (pkin(7) ^ 2 * t168 + pkin(1) ^ 2) + m(4) * (t127 ^ 2 + t85 ^ 2 + t194) + m(5) * (t29 ^ 2 + t30 ^ 2 + t194) + m(6) * (t46 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2); t183 * t62 + t182 * t60 + t157 * t124 + m(6) * (t113 * t46 - t5 * t60 + t6 * t62) + m(7) * (t3 * t62 + t4 * t60 + t54 * t8) + Ifges(3,5) * t139 + Ifges(3,6) * t142 + t125 * t65 + (-t139 * mrSges(3,1) - t142 * mrSges(3,2)) * pkin(7) + t113 * t22 + t54 * t21 + t147 + m(5) * (t124 * t158 + t125 * t83) + (m(4) * (t138 * t85 - t141 * t83) + (-t110 * t138 - t111 * t141) * mrSges(4,3)) * pkin(2); m(6) * (t113 ^ 2 + t167) + m(7) * (t54 ^ 2 + t167) + t124 * t156 + 0.2e1 * t125 * t114 + 0.2e1 * t155 + t113 * t192 + t54 * t193 + m(4) * (t138 ^ 2 + t141 ^ 2) * pkin(2) ^ 2 + m(5) * (t124 ^ 2 * t169 + t125 ^ 2) + t154 + Ifges(3,3) + (-t62 * t105 + t60 * t106) * t200; t126 * t22 + t183 * t80 + t182 * t78 + t157 * pkin(9) + m(6) * (t126 * t46 - t5 * t78 + t6 * t80) + m(7) * (t3 * t80 + t4 * t78 + t63 * t8) + t63 * t21 - pkin(3) * t65 + t147 + m(5) * (-pkin(3) * t83 + pkin(9) * t158); (pkin(9) * t169 + t160) * mrSges(5,3) + m(5) * (-pkin(3) * t125 + pkin(9) * t160) + t154 + t155 + (t125 - pkin(3)) * t114 + (t63 + t54) * t69 + (t126 + t113) * t70 + m(6) * (t113 * t126 + t163) + m(7) * (t54 * t63 + t163) + t195 * ((t60 + t78) * t106 + (-t62 - t80) * t105); -0.2e1 * pkin(3) * t114 + t126 * t192 + t63 * t193 + pkin(9) * t156 + m(6) * (t126 ^ 2 + t166) + m(7) * (t63 ^ 2 + t166) + m(5) * (pkin(9) ^ 2 * t169 + pkin(3) ^ 2) + t154 + (-t80 * t105 + t78 * t106) * t200; m(7) * (t119 * t3 + t122 * t4) + t119 * t35 + t122 * t34 - Ifges(5,6) * t172 + t29 * mrSges(5,1) - t30 * mrSges(5,2) + t3 * mrSges(7,3) - t4 * mrSges(7,1) + t5 * mrSges(6,1) - t6 * mrSges(6,2) + t150 + (m(6) * (t135 * t6 + t136 * t5) + t136 * t33 + t135 * t32) * pkin(4); m(7) * (t119 * t62 + t122 * t60) - t159 * t124 - t60 * mrSges(7,1) - t60 * mrSges(6,1) + t62 * mrSges(7,3) - t62 * mrSges(6,2) + t148 + (-t165 + m(6) * (t135 * t62 - t136 * t60)) * pkin(4); m(7) * (t119 * t80 + t122 * t78) - t159 * pkin(9) - t78 * mrSges(7,1) - t78 * mrSges(6,1) - t80 * mrSges(6,2) + t80 * mrSges(7,3) + t148 + (-t165 + m(6) * (t135 * t80 - t136 * t78)) * pkin(4); -0.2e1 * t122 * mrSges(7,1) + 0.2e1 * t119 * mrSges(7,3) + m(7) * (t119 ^ 2 + t122 ^ 2) + (0.2e1 * mrSges(6,1) * t136 - 0.2e1 * mrSges(6,2) * t135 + m(6) * (t135 ^ 2 + t136 ^ 2) * pkin(4)) * pkin(4) + t197; m(6) * t46 + m(7) * t8 + t21 + t22; m(6) * t113 + m(7) * t54 + t153; m(6) * t126 + m(7) * t63 + t153; 0; m(6) + m(7); m(7) * t4 + t34; m(7) * t60 + t93; m(7) * t78 + t93; m(7) * t122 - mrSges(7,1); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
