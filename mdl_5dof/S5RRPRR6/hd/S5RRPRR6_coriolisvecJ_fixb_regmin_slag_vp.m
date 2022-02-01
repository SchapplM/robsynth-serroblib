% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:25
% EndTime: 2022-01-20 11:17:29
% DurationCPUTime: 1.22s
% Computational Cost: add. (1562->181), mult. (2853->288), div. (0->0), fcn. (1831->8), ass. (0->153)
t93 = sin(pkin(9));
t88 = t93 ^ 2;
t94 = cos(pkin(9));
t171 = t94 ^ 2 + t88;
t100 = cos(qJ(2));
t168 = pkin(1) * qJD(2);
t139 = qJD(1) * t168;
t90 = qJD(1) + qJD(2);
t66 = t90 * qJD(3) + t100 * t139;
t199 = t171 * t66;
t96 = sin(qJ(4));
t164 = qJD(4) * t96;
t144 = t93 * t164;
t95 = sin(qJ(5));
t179 = t95 * t96;
t154 = t93 * t179;
t198 = -qJD(5) * t154 - t95 * t144;
t169 = pkin(1) * qJD(1);
t116 = -t100 * t169 + qJD(3);
t186 = t90 * t93;
t99 = cos(qJ(4));
t163 = qJD(4) * t99;
t165 = qJD(3) * t94;
t166 = t100 * t94;
t72 = -t94 * pkin(3) - t93 * pkin(7) - pkin(2);
t97 = sin(qJ(2));
t197 = t72 * t163 + t99 * t165 - (t99 * t166 + t96 * t97) * t169;
t160 = qJD(4) + qJD(5);
t180 = t94 * t99;
t146 = qJ(3) * t180;
t196 = (t96 * t72 + t146) * qJD(4) + t96 * t165 + (-t96 * t166 + t97 * t99) * t169;
t195 = pkin(4) * t93;
t194 = pkin(8) * t93;
t193 = t100 * pkin(1);
t185 = t90 * t96;
t70 = t90 * qJ(3) + t97 * t169;
t41 = (pkin(4) * t185 + t70) * t93;
t98 = cos(qJ(5));
t176 = t98 * t99;
t153 = t93 * t176;
t129 = t90 * t153;
t43 = t90 * t154 - t129;
t192 = t41 * t43;
t114 = t95 * t99 + t96 * t98;
t44 = t114 * t186;
t191 = t43 * t44;
t182 = t94 * t90;
t78 = -qJD(4) + t182;
t71 = -qJD(5) + t78;
t190 = t71 * t94;
t87 = t90 ^ 2;
t189 = t88 * t87;
t188 = t88 * t90;
t187 = t88 * t99;
t184 = t90 * t99;
t183 = t93 * t99;
t181 = t94 * t96;
t178 = t96 * t66;
t151 = t70 * t180;
t40 = t72 * t90 + t116;
t108 = -t96 * t40 - t151;
t159 = pkin(8) * t186;
t22 = -t96 * t159 - t108;
t177 = t98 * t22;
t143 = t70 * t164;
t126 = t97 * t139;
t149 = t96 * t126 + t40 * t163 + t66 * t180;
t175 = t66 * t187 + (-t94 * t143 + t149) * t94;
t172 = t198 * t90;
t170 = t96 ^ 2 - t99 ^ 2;
t167 = qJ(3) * t96;
t162 = qJD(4) + t78;
t161 = qJ(3) * qJD(4);
t158 = pkin(8) * t183;
t155 = t90 * t183;
t104 = -pkin(8) * t155 - t70 * t181;
t35 = t99 * t40;
t21 = t104 + t35;
t15 = -t78 * pkin(4) + t21;
t18 = qJD(5) * t95 * t22;
t101 = t160 * t114;
t28 = t101 * t93;
t145 = t90 * t163;
t38 = (pkin(4) * t145 + t66) * t93;
t113 = t176 - t179;
t55 = t113 * t93;
t8 = t104 * qJD(4) + t149;
t124 = t99 * t126 - t94 * t178;
t9 = (-t151 + (-t40 + t159) * t96) * qJD(4) + t124;
t157 = (t95 * t9 - t18 + (qJD(5) * t15 + t8) * t98) * t94 - t41 * t28 + t38 * t55;
t156 = t88 * t184;
t86 = t97 * pkin(1) + qJ(3);
t152 = t86 * t181;
t150 = t86 * t180;
t147 = t97 * t168;
t60 = t72 - t193;
t83 = t100 * t168 + qJD(3);
t148 = t96 * t147 + t60 * t163 + t83 * t180;
t140 = -t95 * t8 + t98 * t9;
t138 = t96 * t161;
t137 = pkin(4) * t71 - t15;
t136 = t60 - t194;
t135 = t72 - t194;
t133 = t171 * t83;
t132 = t171 * t100;
t131 = pkin(4) * t155;
t130 = t171 * qJD(3);
t115 = -t95 * t15 - t177;
t2 = t115 * qJD(5) + t140;
t109 = t160 * t153;
t29 = t109 + t198;
t54 = t114 * t93;
t125 = -t2 * t94 + t41 * t29 + t38 * t54;
t123 = t99 * t147 - t83 * t181;
t13 = t108 * qJD(4) + t124;
t122 = -t13 * t94 + (t163 * t70 + t178) * t88;
t121 = t162 * t186;
t120 = qJD(5) * (t135 * t99 + (-pkin(4) - t167) * t94) + (-t94 * t167 - t158) * qJD(4) + t197;
t79 = pkin(8) * t144;
t119 = qJD(5) * (t135 * t96 + t146) - t79 + t196;
t118 = (-qJD(2) + t90) * t169;
t117 = (-qJD(1) - t90) * t168;
t80 = t163 * t195;
t112 = -t116 * t93 - t80;
t111 = t97 * t117;
t110 = t97 * t118;
t107 = -t96 * t60 - t150;
t105 = qJD(4) * t93 * (t78 + t182);
t103 = -t78 ^ 2 - t189;
t102 = t18 + (t22 * t71 - t9) * t95 + t41 * t44;
t23 = t90 * t28;
t85 = t96 * t195;
t68 = t93 * qJ(3) + t85;
t67 = -t90 * pkin(2) + t116;
t61 = -0.2e1 * t88 * t96 * t145;
t56 = t93 * t86 + t85;
t46 = t93 * t83 + t80;
t42 = 0.2e1 * t170 * qJD(4) * t188;
t37 = t99 * t105;
t36 = t96 * t105;
t32 = t136 * t96 + t150;
t27 = t136 * t99 + (-t86 * t96 - pkin(4)) * t94;
t24 = t90 * t109 + t172;
t20 = t107 * qJD(4) + t123 + t79;
t19 = (-t152 - t158) * qJD(4) + t148;
t14 = t43 ^ 2 - t44 ^ 2;
t11 = -t160 * t129 + t43 * t71 - t172;
t10 = -t44 * t71 - t23;
t6 = t24 * t94 + t29 * t71;
t5 = t23 * t94 + t28 * t71;
t4 = -t23 * t55 + t43 * t28;
t3 = t23 * t54 - t55 * t24 + t28 * t44 + t43 * t29;
t1 = [0, 0, 0, 0, t111, t100 * t117, t94 * t111, t90 * t133 + t199, t70 * t133 + t86 * t199 + (t67 + (-pkin(2) - t193) * qJD(1)) * t147, t61, t42, t36, t37, 0, t88 * t83 * t185 - t123 * t78 + (-t107 * t78 + t86 * t156) * qJD(4) + t122, (-qJD(4) * t152 + t148) * t78 + (t83 * t184 + (-t86 * t90 - t70) * t164) * t88 + t175, t4, t3, t5, t6, 0, -(-t95 * t19 + t98 * t20 + (-t27 * t95 - t32 * t98) * qJD(5)) * t71 + t46 * t44 + t56 * t24 + t125, (t98 * t19 + t95 * t20 + (t27 * t98 - t32 * t95) * qJD(5)) * t71 - t46 * t43 - t56 * t23 + t157; 0, 0, 0, 0, t110, t100 * t118, t94 * t110, (-t132 * t169 + t130) * t90 + t199, t70 * t130 + qJ(3) * t199 + ((-pkin(2) * qJD(2) - t67) * t97 - t70 * t132) * t169, t61, t42, t36, t37, 0, t196 * t78 + (t116 * t96 + t99 * t161) * t188 + t122, (-t94 * t138 + t197) * t78 + (-t143 + (t116 * t99 - t138) * t90) * t88 + t175, t4, t3, t5, t6, 0, t68 * t24 + (t119 * t98 + t120 * t95) * t71 - t112 * t44 + t125, -t68 * t23 + (-t119 * t95 + t120 * t98) * t71 + t112 * t43 + t157; 0, 0, 0, 0, 0, 0, 0, -t171 * t87, -t171 * t90 * t70 + t126, 0, 0, 0, 0, 0, t103 * t96, t103 * t99, 0, 0, 0, 0, 0, t101 * t71 + (-t114 * t190 - t93 * t44) * t90, t43 * t186 + (t160 * t71 - t190 * t90) * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, t87 * t96 * t187, -t170 * t189, -t96 * t121, -t99 * t121, 0, t108 * t78 - t70 * t156 + t13, -t35 * t78 + (t162 * t94 + t188) * t96 * t70 - t149, -t191, t14, t10, t11, 0, (-t95 * t21 - t177) * t71 - t44 * t131 + t192 + (t137 * t95 - t177) * qJD(5) + t140, t43 * t131 + (qJD(5) * t137 - t21 * t71 - t8) * t98 + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, t14, t10, t11, 0, t115 * t71 + t192 + t2, (-t8 + (-qJD(5) - t71) * t15) * t98 + t102;];
tauc_reg = t1;
