% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:25:33
% EndTime: 2019-12-05 18:25:42
% DurationCPUTime: 2.02s
% Computational Cost: add. (3045->289), mult. (7473->409), div. (0->0), fcn. (4783->6), ass. (0->175)
t115 = pkin(2) + pkin(1);
t162 = qJD(2) + qJD(4);
t110 = sin(qJ(5));
t113 = cos(qJ(5));
t111 = sin(qJ(4));
t112 = sin(qJ(2));
t193 = pkin(3) + qJ(3);
t85 = t193 * t112;
t79 = qJD(1) * t85;
t124 = t115 * qJD(2) - t79;
t207 = cos(qJ(4));
t114 = cos(qJ(2));
t86 = t193 * t114;
t80 = qJD(1) * t86;
t158 = t207 * t80;
t41 = t111 * t124 + t158;
t35 = t162 * pkin(4) + t41;
t78 = t111 * t114 + t207 * t112;
t71 = t78 * qJD(1);
t89 = t115 * t114;
t75 = -qJD(1) * t89 + qJD(3);
t48 = -pkin(4) * t71 + t75;
t130 = t110 * t35 - t113 * t48;
t121 = t207 * t124;
t166 = qJD(4) * t111;
t145 = qJD(2) * t193;
t167 = qJD(3) * t112;
t126 = -t114 * t145 - t167;
t60 = t126 * qJD(1);
t63 = qJD(3) * t114 - t112 * t145;
t61 = t63 * qJD(1);
t14 = qJD(4) * t121 + t111 * t60 - t80 * t166 + t207 * t61;
t177 = t111 * t112;
t133 = t162 * t177;
t153 = t207 * t114;
t96 = qJD(1) * t153;
t191 = t162 * t96;
t42 = qJD(1) * t133 - t191;
t163 = qJD(1) * qJD(2);
t148 = t112 * t163;
t100 = pkin(1) * t148;
t73 = pkin(2) * t148 + t100;
t32 = pkin(4) * t42 + t73;
t3 = -t130 * qJD(5) + t110 * t32 + t113 * t14;
t171 = qJD(1) * t112;
t151 = t111 * t171;
t69 = -t96 + t151;
t65 = qJD(5) + t69;
t215 = t130 * t65 + t3;
t17 = t110 * t48 + t113 * t35;
t4 = -qJD(5) * t17 - t110 * t14 + t113 * t32;
t214 = t17 * t65 + t4;
t127 = -t110 * t162 - t113 * t71;
t144 = t110 * t65;
t213 = t127 * t144;
t212 = t69 * t162;
t165 = qJD(5) * t110;
t50 = t162 * t78;
t43 = t50 * qJD(1);
t211 = -t113 * t43 + t65 * t165;
t210 = t110 * t130 + t113 * t17;
t22 = -t127 * qJD(5) - t110 * t42;
t209 = t71 ^ 2;
t109 = t114 ^ 2;
t208 = 0.2e1 * t109;
t2 = t3 * t113;
t182 = t111 * t80;
t40 = -t121 + t182;
t45 = -t111 * t79 + t158;
t204 = t40 * t45;
t203 = t40 * t71;
t77 = -t153 + t177;
t202 = t43 * t77;
t201 = t43 * t78;
t140 = t113 * t162;
t52 = t110 * t71 - t140;
t200 = t52 * t69;
t199 = t52 * t71;
t198 = t127 * t52;
t197 = t127 * t71;
t196 = t65 * t71;
t195 = t71 * t69;
t194 = t75 * t71;
t164 = qJD(5) * t113;
t192 = -t110 * t22 - t52 * t164;
t190 = qJD(2) * pkin(1);
t21 = -qJD(5) * t140 + t113 * t42 + t71 * t165;
t188 = t110 * t21;
t186 = t110 * t43;
t185 = t110 * t52;
t184 = t110 * t69;
t183 = t111 * t40;
t180 = t113 * t22;
t179 = t113 * t127;
t178 = t113 * t69;
t117 = qJD(1) ^ 2;
t176 = t112 * t117;
t175 = t114 * t117;
t116 = qJD(2) ^ 2;
t174 = t116 * t112;
t107 = t116 * t114;
t170 = qJD(1) * t114;
t160 = pkin(1) * t170;
t93 = qJD(3) - t160;
t173 = -qJD(3) + t93;
t81 = t115 * t171;
t169 = qJD(2) * t112;
t82 = t115 * t169;
t108 = t112 ^ 2;
t172 = t108 - t109;
t168 = qJD(2) * t114;
t161 = t130 * t178 - t17 * t184 + t2;
t57 = t207 * t60;
t146 = t111 * t61 - t57;
t15 = t41 * qJD(4) + t146;
t159 = t207 * t15;
t84 = -qJ(3) * t171 + t190;
t157 = t84 * t168;
t155 = t78 * t165;
t154 = t78 * t164;
t36 = t40 * t165;
t37 = t40 * t164;
t152 = qJ(3) * t168;
t150 = t130 * t71 + t36;
t149 = t207 * qJD(4);
t147 = t114 * t163;
t143 = t113 * t65;
t141 = (-qJD(3) - t93) * qJD(1);
t139 = t15 * t110 + t17 * t71 + t37;
t138 = t115 * t149;
t137 = t112 * t147;
t128 = -t111 * t86 - t207 * t85;
t56 = -t111 * t85 + t207 * t86;
t20 = t56 * qJD(4) + t111 * t63 - t207 * t126;
t136 = -t128 * t15 + t20 * t40;
t49 = -qJD(2) * t153 - t114 * t149 + t133;
t135 = t15 * t78 - t40 * t49;
t134 = -t49 * t65 + t201;
t131 = t110 * t17 - t113 * t130;
t58 = -pkin(4) * t78 - t89;
t28 = t110 * t58 + t113 * t56;
t27 = -t110 * t56 + t113 * t58;
t129 = -t65 * t184 - t211;
t125 = -t162 * t151 + t191;
t123 = -t131 * qJD(5) - t4 * t110;
t94 = t111 * t115 + pkin(4);
t122 = -t65 * t138 + t40 * t69 - t43 * t94;
t120 = t123 + t2;
t119 = t75 * t69 - t14;
t118 = qJ(3) ^ 2;
t95 = t112 * t175;
t88 = -0.2e1 * t137;
t87 = 0.2e1 * t137;
t83 = t172 * t117;
t74 = (-t152 - t167) * qJD(1);
t68 = -0.2e1 * t172 * t163;
t66 = t69 ^ 2;
t51 = pkin(4) * t69 + t81;
t46 = -t207 * t79 - t182;
t34 = pkin(4) * t49 + t82;
t33 = -t66 + t209;
t29 = t125 + t212;
t26 = pkin(4) * t184 - t113 * t40;
t25 = pkin(4) * t178 + t110 * t40;
t24 = t110 * t51 + t113 * t46;
t23 = -t110 * t46 + t113 * t51;
t19 = t128 * qJD(4) + t111 * t126 + t207 * t63;
t10 = t65 * t143 + t186 + t197;
t9 = t129 + t199;
t8 = t52 * t144 - t180;
t7 = -t127 * t143 - t188;
t6 = -t28 * qJD(5) - t110 * t19 + t113 * t34;
t5 = t27 * qJD(5) + t110 * t34 + t113 * t19;
t1 = (-t21 - t200) * t113 + t213 + t192;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t68, t107, t88, -t174, 0, 0, 0, 0, 0, t87, t68, t107, t88, -t174, 0, -qJ(3) * t107 + (-0.3e1 * t160 + t173) * t169, qJ(3) * t174 + (t173 * t114 + (0.2e1 * t108 - t109) * qJD(1) * pkin(1)) * qJD(2), -t157 - t112 * t74 + (-0.2e1 * t112 * t152 + (t108 + t208) * qJD(3)) * qJD(1), (qJD(1) * qJD(3) * t208 - t157) * qJ(3) + (-qJ(3) * t74 - qJD(3) * t84 + (-0.2e1 * t118 * t170 + (t93 - t160) * pkin(1)) * qJD(2)) * t112, -t42 * t78 - t49 * t71, t42 * t77 + t49 * t69 - t50 * t71 - t201, -t49 * t162, t50 * t69 + t202, -t50 * t162, 0, -t20 * t162 - t89 * t43 + t75 * t50 + t82 * t69 + t73 * t77, -t19 * t162 + t89 * t42 - t75 * t49 + t82 * t71 + t73 * t78, t128 * t42 - t14 * t77 - t19 * t69 + t20 * t71 - t41 * t50 - t43 * t56 + t135, t14 * t56 + t19 * t41 - t73 * t89 + t75 * t82 + t136, t127 * t155 + (t127 * t49 - t21 * t78) * t113, (-t110 * t127 + t113 * t52) * t49 + (t188 - t180 + (t179 + t185) * qJD(5)) * t78, t113 * t134 - t127 * t50 - t155 * t65 - t21 * t77, t52 * t154 + (t22 * t78 - t49 * t52) * t110, -t110 * t134 - t154 * t65 - t22 * t77 - t50 * t52, t50 * t65 + t202, t110 * t135 - t128 * t22 - t130 * t50 + t20 * t52 + t27 * t43 + t37 * t78 + t4 * t77 + t6 * t65, t113 * t135 - t127 * t20 + t128 * t21 - t17 * t50 - t28 * t43 - t3 * t77 - t36 * t78 - t5 * t65, t21 * t27 - t22 * t28 - t5 * t52 + t127 * t6 + t131 * t49 + (-t210 * qJD(5) - t110 * t3 - t113 * t4) * t78, -t130 * t6 + t17 * t5 + t27 * t4 + t28 * t3 + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, t83, 0, t95, 0, 0, 0, 0, 0, 0, -t95, t83, 0, t95, 0, 0, (pkin(1) * t175 + t141) * t112, -pkin(1) * t108 * t117 + t114 * t141, (qJ(3) * t176 + (t84 - t190) * qJD(1)) * t114, (qJ(3) * qJD(1) * t84 + t118 * t176) * t114 + (-t93 * t171 + t74) * pkin(1), t195, t33, t29, -t195, 0, 0, -t80 * t149 + t57 - t81 * t69 - t194 + t45 * t162 + (-t61 + (-t115 * t162 - t124) * qJD(4)) * t111, -t81 * t71 + t119 + (-t138 + t46) * t162, (t41 - t45) * t71 + (t40 + t46) * t69 + (t207 * t42 - t111 * t43 + (t111 * t71 - t207 * t69) * qJD(4)) * t115, -t204 - t41 * t46 - t75 * t81 + (-t159 + t111 * t14 + (t207 * t41 + t183) * qJD(4)) * t115, t7, t1, t10, t8, t9, -t196, -t23 * t65 - t45 * t52 + (t52 * t166 - t207 * t22) * t115 + (-qJD(5) * t65 * t94 - t15) * t113 + t122 * t110 + t150, t45 * t127 + (t165 * t94 + t24) * t65 + (-t127 * t166 + t207 * t21) * t115 + t122 * t113 + t139, -t23 * t127 + t24 * t52 + (-t52 * t138 - t22 * t94 + (-t127 * t94 + t130) * qJD(5)) * t113 + (-t127 * t138 - t21 * t94 - t4 + (t52 * t94 - t17) * qJD(5)) * t110 + t161, t130 * t23 - t17 * t24 - t204 + t120 * t94 + (-t159 + (t210 * t207 + t183) * qJD(4)) * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t148, 0.2e1 * t147, (-t108 - t109) * t117, -qJ(3) * t109 * t117 + t84 * t171 + t100, 0, 0, 0, 0, 0, 0, t71 * t162 + t43, t125 - t212, -t66 - t209, t41 * t69 - t203 + t73, 0, 0, 0, 0, 0, 0, t129 - t199, -t113 * t65 ^ 2 - t186 + t197, (t21 - t200) * t113 - t213 + t192, t215 * t110 + t214 * t113 - t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, t33, t29, -t195, 0, 0, t41 * qJD(2) - t146 - t194, -t162 * t40 + t119, 0, 0, t7, t1, t10, t8, t9, -t196, t40 * t184 - t15 * t113 - t25 * t65 - t41 * t52 + (-t164 * t65 - t186) * pkin(4) + t150, t211 * pkin(4) + t127 * t41 + t40 * t178 + t26 * t65 + t139, -t25 * t127 + t26 * t52 + (-t188 - t180 + (-t179 + t185) * qJD(5)) * pkin(4) + t123 + t161, pkin(4) * t120 + t130 * t25 - t17 * t26 - t40 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, t127 ^ 2 - t52 ^ 2, t52 * t65 - t21, t198, -t127 * t65 - t22, t43, t127 * t40 + t214, t40 * t52 - t215, 0, 0;];
tauc_reg = t11;
