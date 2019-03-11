% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:58
% EndTime: 2019-03-09 01:45:02
% DurationCPUTime: 1.89s
% Computational Cost: add. (3825->285), mult. (7817->382), div. (0->0), fcn. (5057->8), ass. (0->160)
t195 = 2 * qJD(3);
t88 = sin(pkin(9)) * pkin(1) + qJ(3);
t81 = qJD(1) * t88;
t106 = cos(qJ(6));
t101 = sin(pkin(10));
t105 = sin(qJ(4));
t107 = cos(qJ(4));
t173 = cos(pkin(10));
t116 = -t101 * t107 - t173 * t105;
t69 = t116 * qJD(1);
t201 = qJD(6) - t69;
t143 = t106 * t201;
t104 = sin(qJ(6));
t166 = qJD(1) * t105;
t150 = t101 * t166;
t145 = t173 * t107;
t134 = qJD(1) * t145;
t82 = qJD(4) * t134;
t60 = qJD(4) * t150 - t82;
t178 = t104 * t60;
t206 = t143 * t201 - t178;
t159 = t105 * qJD(5);
t163 = qJD(4) * t107;
t160 = t105 * qJD(2);
t86 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t76 = t86 * qJD(1) + qJD(3);
t62 = t76 * t163;
t49 = -qJD(4) * t160 + t62;
t37 = (-qJ(5) * t163 - t159) * qJD(1) + t49;
t149 = t173 * t37;
t156 = t107 * qJD(5);
t157 = t107 * qJD(2);
t199 = (qJ(5) * qJD(1) - t76) * t105 - t157;
t200 = -qJD(1) * t156 + t199 * qJD(4);
t12 = t200 * t101 + t149;
t74 = t116 * qJD(4);
t61 = qJD(1) * t74;
t155 = qJD(1) * qJD(4);
t148 = t107 * t155;
t98 = qJD(3) * qJD(1);
t79 = pkin(4) * t148 + t98;
t26 = -t60 * pkin(5) - t61 * pkin(8) + t79;
t42 = t173 * t199;
t165 = qJD(1) * t107;
t55 = t107 * t76 - t160;
t46 = -qJ(5) * t165 + t55;
t44 = qJD(4) * pkin(4) + t46;
t16 = t101 * t44 - t42;
t14 = qJD(4) * pkin(8) + t16;
t67 = pkin(4) * t166 + qJD(5) + t81;
t72 = t134 - t150;
t25 = -t69 * pkin(5) - t72 * pkin(8) + t67;
t6 = t104 * t25 + t106 * t14;
t2 = -qJD(6) * t6 - t104 * t12 + t106 * t26;
t205 = t201 * t6 + t2;
t122 = t104 * t14 - t106 * t25;
t1 = -t122 * qJD(6) + t104 * t26 + t106 * t12;
t204 = t122 * t201 + t1;
t158 = t106 * qJD(4);
t162 = qJD(6) * t104;
t27 = -qJD(6) * t158 - t106 * t61 + t72 * t162;
t53 = t104 * qJD(4) + t106 * t72;
t164 = qJD(4) * t105;
t71 = -qJD(4) * t145 + t101 * t164;
t183 = t116 * t27 - t53 * t71;
t144 = t104 * t201;
t202 = t53 * t144;
t78 = -t101 * t105 + t145;
t41 = t78 * t60;
t127 = -t201 * t74 + t41;
t153 = t78 * t162;
t198 = -t106 * t127 - t153 * t201;
t56 = t105 * t76 + t157;
t50 = t56 * qJD(4);
t197 = (t105 * t55 - t107 * t56) * qJD(4) - t49 * t105 + t50 * t107;
t196 = t72 ^ 2;
t11 = t101 * t37 - t173 * t200;
t175 = qJ(5) - t86;
t75 = t175 * t105;
t35 = -t101 * t75 + t175 * t145;
t194 = t11 * t35;
t193 = t11 * t116;
t192 = t11 * t78;
t51 = t104 * t72 - t158;
t191 = t51 * t69;
t190 = t53 * t51;
t189 = t53 * t72;
t188 = t72 * t51;
t187 = t72 * t69;
t186 = t74 * t51;
t185 = t74 * t53;
t184 = t116 * t60;
t161 = qJD(6) * t106;
t177 = t104 * t61;
t28 = t53 * qJD(6) + t177;
t182 = -t104 * t28 - t51 * t161;
t176 = t28 * t106;
t181 = -t106 * t186 - t78 * t176;
t180 = t74 * t69 + t41;
t179 = t101 * t199;
t57 = t106 * t60;
t174 = -t105 ^ 2 + t107 ^ 2;
t172 = qJD(6) * t51;
t108 = qJD(4) ^ 2;
t171 = t108 * t105;
t170 = t108 * t107;
t63 = t71 * qJD(4);
t169 = t81 * qJD(1);
t83 = pkin(4) * t163 + qJD(3);
t109 = qJD(1) ^ 2;
t167 = -t108 - t109;
t80 = t105 * pkin(4) + t88;
t154 = pkin(4) * t165;
t152 = t78 * t161;
t151 = t107 * t109 * t105;
t146 = t175 * t107;
t141 = -qJD(6) * t116 + qJD(1);
t139 = t53 * t152;
t138 = t105 * t148;
t90 = t101 * pkin(4) + pkin(8);
t137 = qJD(6) * t201 * t90 + t11;
t136 = -t1 * t116 - t6 * t71;
t135 = t116 * t2 - t122 * t71;
t133 = t104 * t6 - t106 * t122;
t132 = -t104 * t122 - t106 * t6;
t15 = t173 * t44 + t179;
t13 = -qJD(4) * pkin(5) - t15;
t131 = t13 * t74 + t192;
t130 = -t78 * t27 + t185;
t129 = t116 * t28 + t71 * t51;
t128 = t78 * t28 + t186;
t126 = t69 * t71 + t184;
t125 = -t116 * t61 - t72 * t71;
t124 = t78 * t61 + t72 * t74;
t34 = -pkin(5) * t116 - t78 * pkin(8) + t80;
t36 = -t101 * t146 - t173 * t75;
t9 = -t104 * t36 + t106 * t34;
t10 = t104 * t34 + t106 * t36;
t119 = t81 * t195;
t118 = -t57 + (t104 * t69 - t162) * t201;
t117 = t13 * t201 + t60 * t90;
t115 = t175 * t164 - t156;
t114 = -t116 * t12 + t15 * t74 - t16 * t71 - t192;
t112 = t127 * t104 - t152 * t201;
t111 = -t133 * qJD(6) + t1 * t106 - t2 * t104;
t91 = -t173 * pkin(4) - pkin(5);
t68 = t69 ^ 2;
t64 = t74 * qJD(4);
t54 = -qJD(4) * t146 - t159;
t32 = t72 * pkin(5) - t69 * pkin(8) + t154;
t30 = -t71 * pkin(5) - t74 * pkin(8) + t83;
t23 = t101 * t115 + t173 * t54;
t22 = t101 * t54 - t173 * t115;
t19 = t173 * t46 + t179;
t18 = t101 * t46 - t42;
t8 = t104 * t32 + t106 * t19;
t7 = -t104 * t19 + t106 * t32;
t4 = -t10 * qJD(6) - t104 * t23 + t106 * t30;
t3 = t9 * qJD(6) + t104 * t30 + t106 * t23;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t98, t119, -0.2e1 * t138, -0.2e1 * t174 * t155, -t171, 0.2e1 * t138, -t170, 0, t81 * t163 - t86 * t171 + (t105 * t195 + t88 * t163) * qJD(1), -t81 * t164 - t86 * t170 + (t107 * t195 - t88 * t164) * qJD(1), t197, -t197 * t86 + t119, t124, -t125 + t180, t64, t126, t63, 0, -t22 * qJD(4) - t116 * t79 - t80 * t60 - t67 * t71 - t83 * t69, -t23 * qJD(4) + t80 * t61 + t67 * t74 + t83 * t72 + t79 * t78, t22 * t72 + t23 * t69 + t35 * t61 + t36 * t60 - t114, t12 * t36 - t15 * t22 + t16 * t23 + t67 * t83 + t79 * t80 + t194, t130 * t106 - t53 * t153, -t139 + (-t185 + (t27 + t172) * t78) * t104 + t181, t183 + t198, t104 * t128 + t152 * t51, t112 + t129, -t201 * t71 + t184, t104 * t131 + t13 * t152 + t201 * t4 + t22 * t51 + t35 * t28 - t9 * t60 - t135, t10 * t60 + t106 * t131 - t13 * t153 - t201 * t3 + t22 * t53 - t35 * t27 - t136, -t10 * t28 + t9 * t27 - t3 * t51 - t4 * t53 - t133 * t74 + (qJD(6) * t132 - t1 * t104 - t2 * t106) * t78, t1 * t10 - t122 * t4 + t13 * t22 + t2 * t9 + t6 * t3 + t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, t171, 0, t50 * t105 + t49 * t107 + (-t105 * t56 - t107 * t55) * qJD(4), 0, 0, 0, 0, 0, 0, t63, -t64, t125 + t180, t12 * t78 + t15 * t71 + t16 * t74 - t193, 0, 0, 0, 0, 0, 0, t112 - t129, t183 - t198, t139 + (t185 + (-t27 + t172) * t78) * t104 + t181, t111 * t78 - t13 * t71 - t132 * t74 - t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, -t169, 0, 0, 0, 0, 0, 0, t167 * t105, t167 * t107, 0, -t197 - t169, 0, 0, 0, 0, 0, 0, qJD(1) * t69 + t64, -qJD(1) * t72 + t63, -t124 - t126, -t67 * qJD(1) + t114, 0, 0, 0, 0, 0, 0, -t116 * t178 + (t104 * t71 - t106 * t141) * t201 - t128, -t116 * t57 + (t104 * t141 + t106 * t71) * t201 - t130 (t141 * t53 + t129) * t106 + (t141 * t51 + t183) * t104 (t122 * t141 + t136) * t106 + (-t141 * t6 + t135) * t104 - t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, t174 * t109, 0, -t151, 0, 0, -t81 * t165, t81 * t166 - t62 + (t55 + t160) * qJD(4), 0, 0, -t187, -t68 + t196, 0, t187, -t82 + (t72 + t150) * qJD(4), 0, t18 * qJD(4) + t69 * t154 - t67 * t72 - t11, -t149 - t67 * t69 + (-pkin(4) * t72 + qJD(5) * t101) * t165 + (-t101 * (qJ(5) * t166 - t56) + t19) * qJD(4) (t16 - t18) * t72 - (-t15 + t19) * t69 + (t101 * t60 - t173 * t61) * pkin(4), t15 * t18 - t16 * t19 + (t101 * t12 - t173 * t11 - t67 * t165) * pkin(4), -t27 * t104 + t53 * t143 (-t27 + t191) * t106 - t202 + t182, -t189 + t206, t144 * t51 - t176, t118 + t188, -t201 * t72, t104 * t117 - t106 * t137 + t122 * t72 - t18 * t51 - t201 * t7 + t91 * t28, t104 * t137 + t106 * t117 - t18 * t53 + t201 * t8 - t91 * t27 + t6 * t72, t8 * t51 + t7 * t53 + (-t28 * t90 - t122 * t69 + t1 + (t53 * t90 + t122) * qJD(6)) * t106 + (-t27 * t90 + t6 * t69 - t2 + (t51 * t90 - t6) * qJD(6)) * t104, t11 * t91 + t111 * t90 + t122 * t7 - t13 * t18 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 + (t72 - t150) * qJD(4), 0.2e1 * t69 * qJD(4), -t68 - t196, t15 * t72 - t16 * t69 + t79, 0, 0, 0, 0, 0, 0, t118 - t188, -t189 - t206 (t27 + t191) * t106 + t202 + t182, t204 * t104 + t205 * t106 - t13 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, -t51 ^ 2 + t53 ^ 2, t201 * t51 - t27, -t190, -t177 + (-qJD(6) + t201) * t53, -t60, -t13 * t53 + t205, t13 * t51 - t204, 0, 0;];
tauc_reg  = t5;
