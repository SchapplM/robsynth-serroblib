% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:46
% EndTime: 2019-12-31 20:55:52
% DurationCPUTime: 1.97s
% Computational Cost: add. (3922->279), mult. (10613->351), div. (0->0), fcn. (7537->6), ass. (0->166)
t140 = qJD(2) + qJD(3);
t146 = cos(qJ(2));
t216 = cos(qJ(3));
t176 = t216 * t146;
t166 = qJD(1) * t176;
t144 = sin(qJ(3));
t145 = sin(qJ(2));
t184 = qJD(1) * t145;
t175 = t144 * t184;
t100 = -t166 + t175;
t112 = t144 * t146 + t216 * t145;
t102 = qJD(1) * t112;
t143 = sin(pkin(8));
t195 = cos(pkin(8));
t62 = t195 * t100 + t143 * t102;
t196 = t62 * t140;
t191 = t144 * t145;
t160 = t140 * t191;
t186 = t140 * t166;
t67 = qJD(1) * t160 - t186;
t77 = t140 * t112;
t68 = t77 * qJD(1);
t38 = -t143 * t68 - t195 * t67;
t21 = t38 + t196;
t155 = -t143 * t100 + t195 * t102;
t221 = t62 * t155;
t209 = t155 ^ 2;
t200 = qJD(2) * pkin(2);
t138 = t145 * t200;
t220 = 0.2e1 * t138;
t182 = qJD(1) * qJD(2);
t219 = -0.2e1 * t182;
t194 = t100 * qJ(4);
t218 = -pkin(7) - pkin(6);
t122 = t218 * t145;
t116 = qJD(1) * t122;
t123 = t218 * t146;
t118 = qJD(1) * t123;
t177 = t216 * t118;
t74 = -t144 * t116 + t177;
t151 = t74 + t194;
t174 = t216 * qJD(3);
t181 = t143 * t144 * pkin(2);
t103 = t144 * t118;
t75 = t216 * t116 + t103;
t97 = t102 * qJ(4);
t55 = -t97 + t75;
t201 = qJD(3) * t181 + t143 * t151 + (-pkin(2) * t174 + t55) * t195;
t107 = t116 + t200;
t69 = t216 * t107 + t103;
t53 = t69 - t97;
t136 = -t146 * pkin(2) - pkin(1);
t121 = qJD(1) * t136;
t78 = t100 * pkin(3) + qJD(4) + t121;
t28 = t62 * pkin(4) - qJ(5) * t155 + t78;
t213 = t28 * t155;
t80 = t144 * t122 - t216 * t123;
t179 = qJD(2) * t218;
t167 = qJD(1) * t179;
t108 = t145 * t167;
t109 = t146 * t167;
t183 = qJD(3) * t144;
t171 = -t107 * t174 - t216 * t108 - t144 * t109 - t118 * t183;
t16 = -t68 * qJ(4) - t100 * qJD(4) - t171;
t170 = -t144 * t108 + t216 * t109;
t70 = t144 * t107 - t177;
t40 = -qJD(3) * t70 + t170;
t17 = t67 * qJ(4) - t102 * qJD(4) + t40;
t2 = t143 * t16 - t195 * t17;
t79 = t216 * t122 + t144 * t123;
t156 = -t112 * qJ(4) + t79;
t111 = -t176 + t191;
t59 = -t111 * qJ(4) + t80;
t35 = t143 * t59 - t195 * t156;
t217 = t2 * t35;
t215 = t102 * pkin(3);
t54 = t70 - t194;
t50 = t195 * t54;
t25 = t143 * t53 + t50;
t214 = t25 * t155;
t212 = t28 * t62;
t210 = t62 ^ 2;
t117 = t145 * t179;
t119 = t146 * t179;
t46 = -t80 * qJD(3) - t144 * t117 + t216 * t119;
t76 = -qJD(2) * t176 - t146 * t174 + t160;
t149 = t76 * qJ(4) - t112 * qJD(4) + t46;
t45 = t216 * t117 + t144 * t119 + t122 * t174 + t123 * t183;
t29 = -t77 * qJ(4) - t111 * qJD(4) + t45;
t7 = t143 * t29 - t195 * t149;
t207 = t7 * t140;
t206 = t78 * t155;
t205 = t78 * t62;
t8 = t143 * t149 + t195 * t29;
t204 = t8 * t140;
t3 = t143 * t17 + t195 * t16;
t172 = t195 * t144;
t203 = t143 * t55 - t195 * t151 - (t216 * t143 + t172) * qJD(3) * pkin(2);
t202 = -qJD(5) + t201;
t47 = t140 * pkin(3) + t53;
t24 = t143 * t47 + t50;
t42 = -t143 * t76 + t195 * t77;
t199 = t140 * t42;
t198 = t140 * t155;
t197 = t143 * t54;
t193 = t102 * t100;
t192 = t121 * t102;
t148 = qJD(1) ^ 2;
t190 = t146 * t148;
t147 = qJD(2) ^ 2;
t189 = t147 * t145;
t188 = t147 * t146;
t26 = t195 * t53 - t197;
t187 = qJD(5) - t26;
t135 = t216 * pkin(2) + pkin(3);
t96 = pkin(2) * t172 + t143 * t135;
t185 = t145 ^ 2 - t146 ^ 2;
t139 = t140 * qJD(5);
t1 = t139 + t3;
t137 = pkin(2) * t184;
t180 = t145 * t190;
t178 = t203 * t155;
t173 = t145 * t182;
t57 = pkin(2) * t173 + t68 * pkin(3);
t66 = t77 * pkin(3) + t138;
t37 = -t143 * t67 + t195 * t68;
t169 = pkin(1) * t219;
t168 = t26 * t140 - t3;
t165 = t146 * t173;
t23 = t195 * t47 - t197;
t18 = -t140 * pkin(4) + qJD(5) - t23;
t22 = t140 * qJ(5) + t24;
t164 = t155 * t22 + t18 * t62;
t163 = t155 * t24 - t23 * t62;
t72 = t195 * t111 + t143 * t112;
t162 = t37 * t72 + t62 * t42;
t161 = -t209 - t210;
t12 = t209 - t210;
t82 = t111 * pkin(3) + t136;
t159 = t25 * t140 - t2;
t158 = t37 + t198;
t20 = -t37 + t198;
t157 = t121 * t100 + t171;
t33 = pkin(4) * t155 + qJ(5) * t62 + t215;
t154 = t203 * t140 - t2;
t95 = t195 * t135 - t181;
t36 = t143 * t156 + t195 * t59;
t73 = -t143 * t111 + t195 * t112;
t153 = t155 * t7 + t2 * t73 + t35 * t38 - t36 * t37 - t8 * t62;
t152 = -t38 + t196;
t43 = -t143 * t77 - t195 * t76;
t150 = t155 * t42 + t73 * t37 + t38 * t72 + t43 * t62;
t5 = t37 * pkin(4) - t38 * qJ(5) - qJD(5) * t155 + t57;
t133 = -t195 * pkin(3) - pkin(4);
t132 = t143 * pkin(3) + qJ(5);
t88 = -pkin(4) - t95;
t87 = qJ(5) + t96;
t81 = t137 + t215;
t56 = -t100 ^ 2 + t102 ^ 2;
t48 = t186 + (t100 - t175) * t140;
t41 = t43 * t140;
t34 = t72 * pkin(4) - t73 * qJ(5) + t82;
t32 = t137 + t33;
t9 = t42 * pkin(4) - t43 * qJ(5) - t73 * qJD(5) + t66;
t6 = t155 * t43 + t38 * t73;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t165, t185 * t219, t188, -0.2e1 * t165, -t189, 0, -pkin(6) * t188 + t145 * t169, pkin(6) * t189 + t146 * t169, 0, 0, -t102 * t76 - t67 * t112, t76 * t100 - t102 * t77 + t67 * t111 - t112 * t68, -t76 * t140, t100 * t77 + t68 * t111, -t77 * t140, 0, t121 * t77 + t136 * t68 + t46 * t140 + (qJD(1) * t111 + t100) * t138, t102 * t220 - t121 * t76 - t136 * t67 - t45 * t140, -t45 * t100 - t46 * t102 + t111 * t171 - t40 * t112 + t79 * t67 - t80 * t68 + t69 * t76 - t70 * t77, t121 * t220 - t171 * t80 + t40 * t79 + t70 * t45 + t69 * t46, t6, -t150, t41, t162, -t199, 0, t82 * t37 + t78 * t42 + t57 * t72 + t66 * t62 - t207, t155 * t66 + t82 * t38 + t78 * t43 + t57 * t73 - t204, -t23 * t43 - t24 * t42 - t3 * t72 + t153, -t23 * t7 + t24 * t8 + t3 * t36 + t57 * t82 + t78 * t66 + t217, t6, t41, t150, 0, t199, t162, t28 * t42 + t34 * t37 + t5 * t72 + t9 * t62 - t207, -t1 * t72 + t18 * t43 - t22 * t42 + t153, -t155 * t9 - t28 * t43 - t34 * t38 - t5 * t73 + t204, t1 * t36 + t18 * t7 + t22 * t8 + t28 * t9 + t5 * t34 + t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180, t185 * t148, 0, t180, 0, 0, t148 * pkin(1) * t145, pkin(1) * t190, 0, 0, t193, t56, t48, -t193, 0, 0, -t100 * t137 - t192 - t74 * t140 + (t177 + (-pkin(2) * t140 - t107) * t144) * qJD(3) + t170, t75 * t140 + (-t102 * t184 - t140 * t174) * pkin(2) + t157, (t70 + t74) * t102 + (-t69 + t75) * t100 + (t216 * t67 - t144 * t68 + (-t216 * t100 + t102 * t144) * qJD(3)) * pkin(2), -t69 * t74 - t70 * t75 + (-t121 * t184 + t216 * t40 - t144 * t171 + (-t144 * t69 + t216 * t70) * qJD(3)) * pkin(2), t221, t12, t21, -t221, t20, 0, -t81 * t62 + t154 - t206, t201 * t140 - t155 * t81 + t205 - t3, t201 * t62 - t96 * t37 - t95 * t38 + t163 - t178, -t2 * t95 - t201 * t24 + t203 * t23 + t3 * t96 - t78 * t81, t221, t21, -t12, 0, -t20, -t221, -t32 * t62 + t154 - t213, t202 * t62 - t87 * t37 + t88 * t38 + t164 - t178, -t202 * t140 + t155 * t32 + t1 - t212, t1 * t87 - t203 * t18 + t2 * t88 - t202 * t22 - t28 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t56, t48, -t193, 0, 0, t70 * t140 - t192 + t40, t69 * t140 + t157, 0, 0, t221, t12, t21, -t221, t20, 0, -t62 * t215 + t159 - t206, -t155 * t215 + t168 + t205, -t214 + t26 * t62 + (-t143 * t37 - t195 * t38) * pkin(3) + t163, t23 * t25 - t24 * t26 + (-t102 * t78 + t143 * t3 - t195 * t2) * pkin(3), t221, t21, -t12, 0, -t20, -t221, -t33 * t62 + t159 - t213, -t132 * t37 + t133 * t38 - t187 * t62 + t164 - t214, t155 * t33 + 0.2e1 * t139 - t168 - t212, t1 * t132 + t2 * t133 - t18 * t25 + t187 * t22 - t28 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, -t152, t161, t155 * t23 + t24 * t62 + t57, 0, 0, 0, 0, 0, 0, t158, t161, t152, -t155 * t18 + t22 * t62 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, t21, -t140 ^ 2 - t209, -t22 * t140 + t2 + t213;];
tauc_reg = t4;
