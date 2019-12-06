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
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:36:10
% EndTime: 2019-12-05 18:36:15
% DurationCPUTime: 1.28s
% Computational Cost: add. (1566->183), mult. (2869->291), div. (0->0), fcn. (1839->8), ass. (0->155)
t94 = sin(pkin(9));
t89 = t94 ^ 2;
t95 = cos(pkin(9));
t177 = t95 ^ 2 + t89;
t101 = cos(qJ(2));
t174 = pkin(1) * qJD(2);
t142 = qJD(1) * t174;
t91 = qJD(1) + qJD(2);
t66 = t91 * qJD(3) + t101 * t142;
t201 = t177 * t66;
t97 = sin(qJ(4));
t165 = qJD(4) * t97;
t148 = t94 * t165;
t96 = sin(qJ(5));
t184 = t96 * t97;
t158 = t94 * t184;
t200 = -qJD(5) * t158 - t96 * t148;
t175 = pkin(1) * qJD(1);
t117 = -t101 * t175 + qJD(3);
t188 = t91 * t94;
t100 = cos(qJ(4));
t163 = qJD(4) * t100;
t166 = qJD(3) * t95;
t167 = t101 * t95;
t72 = -t95 * pkin(3) - t94 * pkin(7) - pkin(2);
t98 = sin(qJ(2));
t199 = t100 * t166 + t72 * t163 - (t100 * t167 + t97 * t98) * t175;
t162 = qJD(4) + qJD(5);
t169 = t100 * t95;
t146 = qJ(3) * t169;
t198 = (t97 * t72 + t146) * qJD(4) + t97 * t166 + (t100 * t98 - t97 * t167) * t175;
t197 = pkin(4) * t94;
t196 = pkin(8) * t94;
t195 = t98 * pkin(1);
t194 = t101 * pkin(1);
t187 = t91 * t97;
t70 = t91 * qJ(3) + t98 * t175;
t41 = (pkin(4) * t187 + t70) * t94;
t99 = cos(qJ(5));
t168 = t100 * t99;
t151 = t94 * t168;
t130 = t91 * t151;
t43 = t91 * t158 - t130;
t193 = t41 * t43;
t114 = t100 * t96 + t97 * t99;
t44 = t114 * t188;
t192 = t43 * t44;
t186 = t95 * t91;
t79 = -qJD(4) + t186;
t71 = -qJD(5) + t79;
t191 = t71 * t95;
t88 = t91 ^ 2;
t190 = t89 * t88;
t189 = t89 * t91;
t185 = t95 * t97;
t183 = t97 * t66;
t150 = t70 * t169;
t40 = t72 * t91 + t117;
t109 = -t97 * t40 - t150;
t161 = pkin(8) * t188;
t22 = -t97 * t161 - t109;
t182 = t99 * t22;
t147 = t70 * t165;
t127 = t98 * t142;
t156 = t97 * t127 + t40 * t163 + t66 * t169;
t172 = t100 * t89;
t181 = t66 * t172 + (-t95 * t147 + t156) * t95;
t178 = t200 * t91;
t176 = -t100 ^ 2 + t97 ^ 2;
t173 = qJ(3) * t97;
t171 = t100 * t91;
t170 = t100 * t94;
t164 = qJD(4) + t79;
t152 = t91 * t170;
t105 = -pkin(8) * t152 - t70 * t185;
t35 = t100 * t40;
t21 = t105 + t35;
t15 = -t79 * pkin(4) + t21;
t18 = qJD(5) * t96 * t22;
t102 = t162 * t114;
t28 = t102 * t94;
t140 = t91 * t163;
t38 = (pkin(4) * t140 + t66) * t94;
t115 = t168 - t184;
t55 = t115 * t94;
t8 = t105 * qJD(4) + t156;
t125 = t100 * t127 - t95 * t183;
t9 = (-t150 + (-t40 + t161) * t97) * qJD(4) + t125;
t160 = (t96 * t9 - t18 + (qJD(5) * t15 + t8) * t99) * t95 - t41 * t28 + t38 * t55;
t159 = pkin(8) * t170;
t87 = qJ(3) + t195;
t157 = t87 * t185;
t154 = t98 * t174;
t60 = t72 - t194;
t84 = t101 * t174 + qJD(3);
t155 = t97 * t154 + t60 * t163 + t84 * t169;
t153 = t89 * t171;
t149 = t87 * t169;
t143 = -t96 * t8 + t99 * t9;
t141 = qJ(3) * t165;
t139 = pkin(4) * t71 - t15;
t138 = t60 - t196;
t137 = t72 - t196;
t135 = t177 * t84;
t134 = t177 * t101;
t133 = t188 * t195;
t132 = t177 * qJD(3);
t131 = pkin(4) * t152;
t116 = -t96 * t15 - t182;
t2 = t116 * qJD(5) + t143;
t110 = t162 * t151;
t29 = t110 + t200;
t54 = t114 * t94;
t126 = -t2 * t95 + t41 * t29 + t38 * t54;
t124 = t100 * t154 - t84 * t185;
t13 = t109 * qJD(4) + t125;
t123 = -t13 * t95 + (t163 * t70 + t183) * t89;
t122 = t164 * t188;
t121 = qJD(5) * ((-pkin(4) - t173) * t95 + t137 * t100) + (-t95 * t173 - t159) * qJD(4) + t199;
t80 = pkin(8) * t148;
t120 = qJD(5) * (t137 * t97 + t146) - t80 + t198;
t119 = (-qJD(2) + t91) * t175;
t118 = (-qJD(1) - t91) * t174;
t81 = t163 * t197;
t113 = -t117 * t94 - t81;
t112 = t98 * t118;
t111 = t98 * t119;
t108 = -t97 * t60 - t149;
t107 = qJD(4) * t94 * (t79 + t186);
t104 = -t79 ^ 2 - t190;
t103 = t18 + (t22 * t71 - t9) * t96 + t41 * t44;
t23 = t91 * t28;
t86 = t97 * t197;
t76 = t94 * t127;
t68 = t94 * qJ(3) + t86;
t67 = -t91 * pkin(2) + t117;
t61 = -0.2e1 * t89 * t97 * t140;
t56 = t94 * t87 + t86;
t46 = t94 * t84 + t81;
t42 = 0.2e1 * t176 * qJD(4) * t189;
t37 = t100 * t107;
t36 = t97 * t107;
t32 = t138 * t97 + t149;
t27 = (-t87 * t97 - pkin(4)) * t95 + t138 * t100;
t24 = t91 * t110 + t178;
t20 = t108 * qJD(4) + t124 + t80;
t19 = (-t157 - t159) * qJD(4) + t155;
t14 = t43 ^ 2 - t44 ^ 2;
t11 = -t162 * t130 + t43 * t71 - t178;
t10 = -t44 * t71 - t23;
t6 = t24 * t95 + t29 * t71;
t5 = t23 * t95 + t28 * t71;
t4 = -t23 * t55 + t43 * t28;
t3 = t23 * t54 - t55 * t24 + t28 * t44 + t43 * t29;
t1 = [0, 0, 0, 0, t112, t101 * t118, t95 * t112, qJD(2) * t133 + t76, t91 * t135 + t201, t70 * t135 + t87 * t201 + (t67 + (-pkin(2) - t194) * qJD(1)) * t154, t61, t42, t36, t37, 0, t89 * t84 * t187 - t124 * t79 + (-t108 * t79 + t87 * t153) * qJD(4) + t123, (-qJD(4) * t157 + t155) * t79 + (t84 * t171 + (-t87 * t91 - t70) * t165) * t89 + t181, t4, t3, t5, t6, 0, -(-t96 * t19 + t99 * t20 + (-t27 * t96 - t32 * t99) * qJD(5)) * t71 + t46 * t44 + t56 * t24 + t126, (t99 * t19 + t96 * t20 + (t27 * t99 - t32 * t96) * qJD(5)) * t71 - t46 * t43 - t56 * t23 + t160; 0, 0, 0, 0, t111, t101 * t119, t95 * t111, -qJD(1) * t133 + t76, (-t134 * t175 + t132) * t91 + t201, t70 * t132 + qJ(3) * t201 + ((-pkin(2) * qJD(2) - t67) * t98 - t70 * t134) * t175, t61, t42, t36, t37, 0, t198 * t79 + (qJ(3) * t163 + t117 * t97) * t189 + t123, (-t95 * t141 + t199) * t79 + (-t147 + (t117 * t100 - t141) * t91) * t89 + t181, t4, t3, t5, t6, 0, t68 * t24 + (t120 * t99 + t121 * t96) * t71 - t113 * t44 + t126, -t68 * t23 + (-t120 * t96 + t121 * t99) * t71 + t113 * t43 + t160; 0, 0, 0, 0, 0, 0, 0, 0, -t177 * t88, -t177 * t91 * t70 + t127, 0, 0, 0, 0, 0, t104 * t97, t104 * t100, 0, 0, 0, 0, 0, t102 * t71 + (-t114 * t191 - t94 * t44) * t91, t43 * t188 + (t162 * t71 - t191 * t91) * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 * t97 * t172, -t176 * t190, -t97 * t122, -t100 * t122, 0, t109 * t79 - t70 * t153 + t13, -t35 * t79 + (t164 * t95 + t189) * t97 * t70 - t156, -t192, t14, t10, t11, 0, (-t96 * t21 - t182) * t71 - t44 * t131 + t193 + (t139 * t96 - t182) * qJD(5) + t143, t43 * t131 + (t139 * qJD(5) - t21 * t71 - t8) * t99 + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192, t14, t10, t11, 0, t116 * t71 + t193 + t2, (-t8 + (-qJD(5) - t71) * t15) * t99 + t103;];
tauc_reg = t1;
