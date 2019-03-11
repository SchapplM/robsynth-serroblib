% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRR11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:32:49
% EndTime: 2019-03-09 14:32:57
% DurationCPUTime: 2.43s
% Computational Cost: add. (3086->238), mult. (6887->427), div. (0->0), fcn. (6139->8), ass. (0->154)
t120 = cos(qJ(2));
t105 = qJD(2) * t120;
t116 = sin(qJ(2));
t113 = sin(qJ(6));
t117 = cos(qJ(6));
t118 = cos(qJ(5));
t119 = cos(qJ(4));
t180 = t118 * t119;
t114 = sin(qJ(5));
t115 = sin(qJ(4));
t183 = t114 * t115;
t130 = -t180 + t183;
t82 = t114 * t119 + t118 * t115;
t52 = -t113 * t82 - t117 * t130;
t200 = qJD(4) + qJD(5);
t53 = t200 * t82;
t170 = qJD(5) * t114;
t174 = qJD(4) * t115;
t54 = -t114 * t174 - t115 * t170 + t180 * t200;
t15 = -t52 * qJD(6) + t113 * t53 - t117 * t54;
t51 = -t113 * t130 + t117 * t82;
t201 = -t51 * t105 + t15 * t116;
t195 = pkin(3) + pkin(7);
t112 = t120 ^ 2;
t141 = qJD(2) * (t116 ^ 2 - t112);
t109 = t115 ^ 2;
t178 = -t119 ^ 2 + t109;
t140 = t178 * qJD(4);
t124 = -t51 * qJD(6) - t113 * t54 - t117 * t53;
t199 = t52 * t105 + t116 * t124;
t121 = -pkin(2) - pkin(8);
t181 = t116 * qJ(3);
t131 = -t120 * t121 + t181;
t166 = t120 * qJD(3);
t93 = t195 * t120;
t198 = qJD(2) * t131 - t93 * qJD(4) - t166;
t197 = t200 * t120;
t196 = 0.2e1 * qJD(3);
t194 = t116 * pkin(4);
t193 = t116 * pkin(5);
t192 = t118 * pkin(4);
t191 = pkin(9) - t121;
t75 = -pkin(1) - t131;
t148 = pkin(9) * t120 - t75;
t92 = t195 * t116;
t85 = t119 * t92;
t42 = t148 * t115 + t194 + t85;
t179 = t119 * t120;
t84 = t115 * t92;
t189 = t119 * t75 + t84;
t45 = -pkin(9) * t179 + t189;
t43 = t118 * t45;
t190 = t114 * t42 + t43;
t88 = t191 * t115;
t89 = t191 * t119;
t188 = -t114 * t89 - t118 * t88;
t60 = -t118 * t179 + t120 * t183;
t14 = t60 * pkin(10) + t190;
t186 = t117 * t14;
t185 = qJ(3) * t120;
t184 = t113 * t114;
t182 = t114 * t117;
t100 = t115 * pkin(4) + qJ(3);
t176 = qJD(2) * t116;
t175 = qJD(2) * t119;
t173 = qJD(4) * t119;
t172 = qJD(4) * t120;
t171 = qJD(4) * t121;
t169 = qJD(5) * t118;
t168 = qJD(6) * t113;
t167 = qJD(6) * t117;
t96 = pkin(4) * t173 + qJD(3);
t165 = -0.2e1 * pkin(1) * qJD(2);
t164 = pkin(5) * t105;
t163 = pkin(4) * t170;
t162 = pkin(4) * t169;
t161 = pkin(5) * t168;
t160 = pkin(5) * t167;
t159 = pkin(7) * t176;
t67 = pkin(4) * t179 + t93;
t32 = t130 * t197 + t82 * t176;
t139 = pkin(2) * t176 - t116 * qJD(3);
t57 = (pkin(8) * t116 - t185) * qJD(2) + t139;
t102 = pkin(7) * t105;
t87 = pkin(3) * t105 + t102;
t143 = -t115 * t57 + t119 * t87;
t21 = (-pkin(9) * t115 * t116 + pkin(4) * t120) * qJD(2) + (t148 * t119 - t84) * qJD(4) + t143;
t151 = t116 * t175;
t155 = t115 * t172;
t127 = t151 + t155;
t27 = -t115 * t87 - t119 * t57 - t92 * t173 + t75 * t174;
t25 = t127 * pkin(9) - t27;
t146 = -t114 * t25 + t118 * t21;
t9 = -t190 * qJD(5) + t146;
t4 = -t32 * pkin(10) + t164 + t9;
t153 = t115 * t176;
t33 = -t114 * t153 + t118 * t151 + t197 * t82;
t8 = -t114 * t21 - t118 * t25 - t42 * t169 + t45 * t170;
t5 = t33 * pkin(10) - t8;
t156 = -t113 * t5 + t117 * t4;
t154 = t119 * t172;
t152 = t116 * t105;
t150 = t115 * t173;
t145 = -t114 * t45 + t118 * t42;
t61 = t82 * t120;
t13 = t61 * pkin(10) + t145 + t193;
t149 = -t13 - t193;
t12 = t14 * t168;
t147 = -t113 * t4 + t12;
t144 = t114 * t88 - t118 * t89;
t101 = pkin(5) + t192;
t142 = qJD(6) * (-pkin(5) - t101);
t138 = t115 * t151;
t137 = -t120 * pkin(2) - t181;
t136 = t113 * t13 + t186;
t40 = pkin(10) * t130 + t144;
t41 = -t82 * pkin(10) + t188;
t135 = t113 * t41 - t117 * t40;
t134 = t113 * t40 + t117 * t41;
t38 = -t113 * t61 - t117 * t60;
t39 = t113 * t60 - t117 * t61;
t128 = -t82 * t105 - t54 * t116;
t80 = t191 * t174;
t81 = qJD(4) * t89;
t29 = -t114 * t80 + t118 * t81 + t89 * t169 - t88 * t170;
t86 = t195 * t176;
t126 = -t86 + (-t116 * t121 - t185) * qJD(4);
t125 = t137 * qJD(2) + t166;
t2 = -t136 * qJD(6) + t156;
t30 = -t188 * qJD(5) + t114 * t81 + t118 * t80;
t55 = -pkin(4) * t155 + (-pkin(4) * t119 - t195) * t176;
t123 = (t114 * t168 + (-t117 * t118 + t184) * qJD(5)) * pkin(4);
t122 = (-t114 * t167 + (-t113 * t118 - t182) * qJD(5)) * pkin(4);
t95 = 0.2e1 * t152;
t90 = -pkin(1) + t137;
t71 = -t115 * t105 - t116 * t173;
t70 = t119 * t105 - t116 * t174;
t63 = -qJ(3) * t105 + t139;
t59 = t82 * pkin(5) + t100;
t48 = -t60 * pkin(5) + t67;
t47 = -t101 * t168 + t122;
t46 = -t101 * t167 + t123;
t44 = t54 * pkin(5) + t96;
t34 = -t105 * t130 - t53 * t116;
t28 = -t189 * qJD(4) + t143;
t26 = -t33 * pkin(5) + t55;
t23 = t53 * pkin(10) + t30;
t22 = -t54 * pkin(10) - t29;
t11 = t39 * qJD(6) + t113 * t32 - t117 * t33;
t10 = -t38 * qJD(6) + t113 * t33 + t117 * t32;
t7 = -t134 * qJD(6) - t113 * t22 + t117 * t23;
t6 = t135 * qJD(6) - t113 * t23 - t117 * t22;
t1 = (-qJD(6) * t13 - t5) * t117 + t147;
t3 = [0, 0, 0, t95, -0.2e1 * t141, 0, 0, 0, t116 * t165, t120 * t165, 0, 0.2e1 * t63 * t120 - 0.2e1 * t90 * t176, -0.2e1 * t90 * t105 - 0.2e1 * t63 * t116, 0.2e1 * t90 * t63, -0.2e1 * t109 * t152 + 0.2e1 * t112 * t150, -0.2e1 * t112 * t140 - 0.4e1 * t120 * t138, 0.2e1 * t115 * t141 - 0.2e1 * t116 * t154, 0.2e1 * t116 * t155 + 0.2e1 * t119 * t141, t95, 0.2e1 * (-t93 * t175 + t28) * t116 + 0.2e1 * ((-t115 * t75 + t85) * qJD(2) - t86 * t119 - t93 * t174) * t120, 0.2e1 * (t93 * t115 * qJD(2) + t27) * t116 + 0.2e1 * (-qJD(2) * t189 + t86 * t115 - t173 * t93) * t120, -0.2e1 * t61 * t32, 0.2e1 * t32 * t60 - 0.2e1 * t61 * t33, -0.2e1 * t105 * t61 + 0.2e1 * t32 * t116, 0.2e1 * t105 * t60 + 0.2e1 * t33 * t116, t95, 0.2e1 * t105 * t145 + 0.2e1 * t9 * t116 - 0.2e1 * t67 * t33 - 0.2e1 * t55 * t60, -0.2e1 * t105 * t190 + 0.2e1 * t8 * t116 + 0.2e1 * t67 * t32 - 0.2e1 * t55 * t61, 0.2e1 * t39 * t10, -0.2e1 * t10 * t38 - 0.2e1 * t39 * t11, 0.2e1 * t10 * t116 + 0.2e1 * t105 * t39, -0.2e1 * t105 * t38 - 0.2e1 * t11 * t116, t95, 0.2e1 * t2 * t116 + 0.2e1 * (-t113 * t14 + t117 * t13) * t105 + 0.2e1 * t26 * t38 + 0.2e1 * t48 * t11, 0.2e1 * t1 * t116 + 0.2e1 * t48 * t10 - 0.2e1 * t105 * t136 + 0.2e1 * t26 * t39; 0, 0, 0, 0, 0, t105, -t176, 0, -t102, t159, t125, t102, -t159, t125 * pkin(7), t120 * t140 + t138, 0.4e1 * t120 * t150 - t178 * t176, t70, t71, 0, t126 * t115 - t119 * t198, t115 * t198 + t126 * t119, -t130 * t32 + t61 * t53, -t130 * t33 - t32 * t82 - t53 * t60 + t61 * t54, t34, t128, 0, -t100 * t33 + t105 * t144 + t30 * t116 + t67 * t54 + t55 * t82 - t96 * t60, t100 * t32 - t105 * t188 + t29 * t116 - t130 * t55 - t67 * t53 - t96 * t61, t10 * t52 + t124 * t39, -t10 * t51 - t52 * t11 - t124 * t38 + t15 * t39, t199, t201, 0, -t105 * t135 + t59 * t11 + t7 * t116 - t15 * t48 + t26 * t51 + t44 * t38, t59 * t10 - t105 * t134 + t6 * t116 + t124 * t48 + t26 * t52 + t44 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, qJ(3) * t196, -0.2e1 * t150, 0.2e1 * t140, 0, 0, 0, 0.2e1 * qJ(3) * t173 + 0.2e1 * qJD(3) * t115, -0.2e1 * qJ(3) * t174 + 0.2e1 * qJD(3) * t119, 0.2e1 * t130 * t53, 0.2e1 * t130 * t54 + 0.2e1 * t53 * t82, 0, 0, 0, 0.2e1 * t100 * t54 + 0.2e1 * t96 * t82, -0.2e1 * t100 * t53 - 0.2e1 * t130 * t96, 0.2e1 * t52 * t124, -0.2e1 * t124 * t51 + 0.2e1 * t15 * t52, 0, 0, 0, -0.2e1 * t15 * t59 + 0.2e1 * t44 * t51, 0.2e1 * t124 * t59 + 0.2e1 * t44 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, 0, 0, t102, 0, 0, 0, 0, 0, t70, t71, 0, 0, 0, 0, 0, t34, t128, 0, 0, 0, 0, 0, t199, t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153 - t154, t127, t105, t28, t27, 0, 0, t32, t33, t105, t105 * t192 + (-t43 + (-t42 - t194) * t114) * qJD(5) + t146 (-t105 * t114 - t116 * t169) * pkin(4) + t8, 0, 0, t10, -t11, t105, t47 * t116 + (-pkin(4) * t184 + t117 * t101) * t105 + t2, t46 * t116 - (pkin(4) * t182 + t113 * t101) * t105 - t117 * t5 - t13 * t167 + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174, -t173, 0, -t115 * t171, -t119 * t171, 0, 0, -t53, -t54, 0, t30, t29, 0, 0, t124, t15, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174, -t173, 0, 0, 0, 0, 0, -t53, -t54, 0, 0, 0, 0, 0, t124, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t163, -0.2e1 * t162, 0, 0, 0, 0, 0, 0.2e1 * t47, 0.2e1 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t33, t105, t9, t8, 0, 0, t10, -t11, t105, t117 * t164 + (t113 * t149 - t186) * qJD(6) + t156, t12 + (-t4 - t164) * t113 + (qJD(6) * t149 - t5) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t54, 0, t30, t29, 0, 0, t124, t15, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t54, 0, 0, 0, 0, 0, t124, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, -t162, 0, 0, 0, 0, 0, t113 * t142 + t122, t117 * t142 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t161, -0.2e1 * t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t105, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t15, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161, -t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
