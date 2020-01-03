% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:26
% EndTime: 2019-12-31 21:20:35
% DurationCPUTime: 2.90s
% Computational Cost: add. (4552->332), mult. (11157->428), div. (0->0), fcn. (7566->6), ass. (0->190)
t134 = sin(qJ(3));
t135 = sin(qJ(2));
t137 = cos(qJ(2));
t227 = cos(qJ(3));
t107 = t134 * t137 + t227 * t135;
t190 = qJD(1) * t107;
t238 = qJD(5) + t190;
t129 = qJD(2) + qJD(3);
t133 = sin(qJ(5));
t136 = cos(qJ(5));
t179 = t227 * t137;
t165 = qJD(1) * t179;
t189 = qJD(1) * t135;
t178 = t134 * t189;
t94 = -t165 + t178;
t76 = t133 * t129 - t136 * t94;
t172 = t238 * t76;
t186 = qJD(5) * t136;
t187 = qJD(5) * t133;
t75 = t129 * t107;
t67 = t75 * qJD(1);
t35 = t129 * t187 - t133 * t67 - t94 * t186;
t242 = t35 - t172;
t78 = t136 * t129 + t133 * t94;
t235 = t238 * t78;
t201 = qJD(5) * t78;
t36 = -t136 * t67 + t201;
t241 = -t36 + t235;
t231 = -pkin(7) - pkin(6);
t112 = t231 * t135;
t109 = qJD(1) * t112;
t220 = qJD(2) * pkin(2);
t102 = t109 + t220;
t113 = t231 * t137;
t110 = qJD(1) * t113;
t200 = t134 * t110;
t70 = -t227 * t102 - t200;
t194 = qJD(4) + t70;
t229 = t190 * pkin(4);
t195 = t229 + t194;
t232 = pkin(3) + pkin(8);
t37 = -t232 * t129 + t195;
t125 = -t137 * pkin(2) - pkin(1);
t111 = qJD(1) * t125;
t144 = -qJ(4) * t190 + t111;
t39 = t232 * t94 + t144;
t14 = t133 * t37 + t136 * t39;
t185 = qJD(1) * qJD(2);
t174 = t135 * t185;
t199 = t134 * t135;
t161 = t129 * t199;
t192 = t129 * t165;
t66 = qJD(1) * t161 - t192;
t148 = pkin(2) * t174 + t66 * qJ(4) - qJD(4) * t190;
t11 = t232 * t67 + t148;
t180 = qJD(2) * t231;
t166 = qJD(1) * t180;
t103 = t135 * t166;
t104 = t137 * t166;
t175 = qJD(3) * t227;
t188 = qJD(3) * t134;
t34 = t102 * t188 + t134 * t103 - t227 * t104 - t110 * t175;
t24 = -t66 * pkin(4) + t34;
t2 = -qJD(5) * t14 - t133 * t11 + t136 * t24;
t240 = t14 * t238 + t2;
t127 = t135 * t220;
t239 = 0.2e1 * t127;
t237 = -0.2e1 * t185;
t171 = -t102 * t175 - t227 * t103 - t134 * t104 - t110 * t188;
t31 = -t129 * qJD(4) + t171;
t18 = -t67 * pkin(4) - t31;
t230 = t94 * pkin(4);
t100 = t227 * t110;
t71 = t134 * t102 - t100;
t62 = -t129 * qJ(4) - t71;
t43 = -t62 - t230;
t236 = t18 * t133 + t43 * t186;
t73 = t227 * t109 + t200;
t202 = pkin(2) * t175 + qJD(4) - t73;
t51 = t66 * t107;
t176 = qJD(2) * t227;
t74 = t161 + (-t175 - t176) * t137;
t234 = -t190 * t74 - t51;
t157 = t133 * t39 - t136 * t37;
t1 = -qJD(5) * t157 + t136 * t11 + t133 * t24;
t233 = t190 ^ 2;
t12 = t157 * t187;
t215 = t157 * t133;
t228 = -t190 * t215 - t12;
t79 = -t227 * t112 - t134 * t113;
t226 = t34 * t79;
t57 = t94 * pkin(3) + t144;
t225 = t57 * t190;
t224 = t78 * t76;
t223 = t238 * t94;
t222 = t94 * t190;
t219 = t111 * t190;
t218 = t129 * t74;
t217 = t129 * t75;
t216 = t129 * t94;
t214 = t133 * t66;
t213 = t133 * t75;
t212 = t133 * t78;
t211 = t136 * t35;
t210 = t136 * t43;
t64 = t136 * t66;
t209 = t136 * t76;
t208 = t136 * t238;
t16 = t18 * t136;
t207 = t36 * t133;
t155 = t231 * t176;
t167 = t134 * t180;
t40 = -t112 * t175 - t113 * t188 - t135 * t155 - t137 * t167;
t206 = t40 * t129;
t80 = t134 * t112 - t227 * t113;
t41 = qJD(3) * t80 + t135 * t167 - t137 * t155;
t205 = t41 * t129;
t203 = t229 + t202;
t140 = qJD(1) ^ 2;
t198 = t137 * t140;
t139 = qJD(2) ^ 2;
t197 = t139 * t135;
t196 = t139 * t137;
t191 = t135 ^ 2 - t137 ^ 2;
t184 = pkin(2) * t188;
t183 = -t94 ^ 2 + t233;
t182 = t135 * t198;
t181 = -t14 * t190 - t2;
t177 = -t14 * t94 + t16;
t173 = t238 * t43;
t170 = -t36 + t201;
t169 = t133 * t238;
t168 = pkin(1) * t237;
t72 = t134 * t109 - t100;
t124 = -t227 * pkin(2) - pkin(3);
t164 = t137 * t174;
t163 = -t72 + t184;
t68 = pkin(3) * t190 + t94 * qJ(4);
t106 = -t179 + t199;
t162 = t67 * t106 + t94 * t75;
t159 = t133 * t14 - t136 * t157;
t158 = t136 * t14 + t215;
t154 = -t107 * qJ(4) + t125;
t52 = t232 * t106 + t154;
t60 = t107 * pkin(4) + t79;
t27 = t133 * t60 + t136 * t52;
t26 = -t133 * t52 + t136 * t60;
t156 = -t157 * t94 + t190 * t210 + t236;
t59 = pkin(2) * t189 + t68;
t153 = t34 + t225;
t152 = -t169 * t238 - t64;
t151 = t111 * t94 + t171;
t150 = t71 * t129 - t34;
t149 = -t187 * t238 - t64;
t147 = t74 * qJ(4) - t107 * qJD(4) + t127;
t146 = -t57 * t94 - t31;
t145 = -t208 * t238 + t214;
t143 = t66 * t106 - t107 * t67 - t190 * t75 + t74 * t94;
t45 = -t129 * t178 + t192 + t216;
t142 = t34 * t107 + t190 * t41 + t40 * t94 - t79 * t66 - t80 * t67;
t141 = qJD(5) * t158 + t1 * t133 + t2 * t136;
t122 = t134 * pkin(2) + qJ(4);
t121 = -pkin(8) + t124;
t90 = t190 * pkin(8);
t69 = t106 * pkin(3) + t154;
t61 = -t106 * pkin(4) + t80;
t58 = -t129 * pkin(3) + t194;
t55 = t72 - t230;
t50 = t71 - t230;
t48 = t68 + t90;
t46 = -t129 * t190 + t67;
t44 = t66 - t216;
t42 = t59 + t90;
t30 = t75 * pkin(3) + t147;
t29 = -t74 * pkin(4) + t41;
t28 = -t75 * pkin(4) - t40;
t25 = t67 * pkin(3) + t148;
t23 = t133 * t50 + t136 * t48;
t22 = -t133 * t48 + t136 * t50;
t21 = t133 * t55 + t136 * t42;
t20 = -t133 * t42 + t136 * t55;
t17 = t232 * t75 + t147;
t10 = -t76 * t94 + t145;
t9 = t78 * t94 + t152;
t8 = t136 * t172 + t207;
t7 = -t169 * t78 - t211;
t5 = (-t36 - t235) * t136 + (t35 + t172) * t133;
t4 = -qJD(5) * t27 - t133 * t17 + t136 * t29;
t3 = qJD(5) * t26 + t133 * t29 + t136 * t17;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t164, t191 * t237, t196, -0.2e1 * t164, -t197, 0, -pkin(6) * t196 + t135 * t168, pkin(6) * t197 + t137 * t168, 0, 0, t234, t143, -t218, t162, -t217, 0, t111 * t75 + t125 * t67 - t205 + (qJD(1) * t106 + t94) * t127, -t111 * t74 - t125 * t66 + t190 * t239 + t206, t106 * t171 - t70 * t74 - t71 * t75 + t142, t111 * t239 - t171 * t80 - t71 * t40 + t70 * t41 + t226, 0, t218, t217, t234, t143, t162, t31 * t106 - t58 * t74 + t62 * t75 + t142, -t25 * t106 - t30 * t94 - t57 * t75 - t69 * t67 + t205, -t25 * t107 - t190 * t30 + t57 * t74 + t69 * t66 - t206, t25 * t69 + t57 * t30 - t31 * t80 + t62 * t40 + t58 * t41 + t226, t75 * t212 + (-t35 * t133 + t186 * t78) * t106, (-t133 * t76 + t136 * t78) * t75 + (-t207 - t211 + (-t209 - t212) * qJD(5)) * t106, t238 * t213 - t35 * t107 - t78 * t74 + (t186 * t238 - t214) * t106, -t75 * t209 + (-t136 * t36 + t187 * t76) * t106, t106 * t149 - t36 * t107 + t208 * t75 + t76 * t74, -t238 * t74 - t51, -t75 * t210 + t2 * t107 + t157 * t74 - t26 * t66 + t28 * t76 + t61 * t36 + t4 * t238 + (t187 * t43 - t16) * t106, -t1 * t107 + t236 * t106 + t14 * t74 + t43 * t213 - t238 * t3 + t27 * t66 + t28 * t78 - t61 * t35, t26 * t35 - t27 * t36 - t3 * t76 - t4 * t78 + t158 * t75 + (-qJD(5) * t159 + t1 * t136 - t133 * t2) * t106, t1 * t27 + t14 * t3 - t157 * t4 + t18 * t61 + t2 * t26 + t43 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, t191 * t140, 0, t182, 0, 0, t140 * pkin(1) * t135, pkin(1) * t198, 0, 0, t222, t183, t45, -t222, 0, 0, -t219 + t72 * t129 + (-t129 * t188 - t94 * t189) * pkin(2) - t34, t73 * t129 + (-t129 * t175 - t189 * t190) * pkin(2) + t151, (t71 - t72) * t190 + (t70 + t73) * t94 + (t227 * t66 - t134 * t67 + (t134 * t190 - t227 * t94) * qJD(3)) * pkin(2), -t70 * t72 - t71 * t73 + (-t111 * t189 - t227 * t34 - t134 * t171 + (t134 * t70 + t227 * t71) * qJD(3)) * pkin(2), 0, t44, t46, t222, t183, -t222, -t122 * t67 - t124 * t66 + (t163 - t62) * t190 + (t58 - t202) * t94, t129 * t163 + t59 * t94 + t153, t129 * t202 + t190 * t59 + t146, -t31 * t122 + t34 * t124 + t163 * t58 - t202 * t62 - t57 * t59, t7, t5, t9, t8, t10, t223, -t121 * t64 + t122 * t36 + t203 * t76 + (-t121 * t187 + t136 * t184 - t20) * t238 + t156, -t122 * t35 + (-t121 * t186 + t21) * t238 + t203 * t78 + (t121 * t66 - t184 * t238 - t173) * t133 + t177, t20 * t78 + t21 * t76 + (t121 * t170 - t184 * t76 - t1) * t133 + (-t78 * t184 + t121 * t35 + (-t121 * t76 - t14) * qJD(5) + t181) * t136 + t228, t121 * t141 + t18 * t122 - t14 * t21 + t157 * t20 + t159 * t184 + t203 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, t183, t45, -t222, 0, 0, t150 - t219, -t70 * t129 + t151, 0, 0, 0, t44, t46, t222, t183, -t222, pkin(3) * t66 - qJ(4) * t67 + (-t62 - t71) * t190 + (t58 - t194) * t94, t68 * t94 - t150 + t225, t129 * t194 + t190 * t68 + t146, -t34 * pkin(3) - t31 * qJ(4) - t194 * t62 - t57 * t68 - t58 * t71, t7, t5, t9, t8, t10, t223, qJ(4) * t36 - t149 * t232 + t195 * t76 - t22 * t238 + t156, -qJ(4) * t35 + (t186 * t232 + t23) * t238 + t195 * t78 + (-t232 * t66 - t173) * t133 + t177, t22 * t78 + t23 * t76 + (-t170 * t232 - t1) * t133 + (-t232 * t35 + (t232 * t76 - t14) * qJD(5) + t181) * t136 + t228, t18 * qJ(4) - t14 * t23 - t141 * t232 + t157 * t22 + t195 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t222, -t129 ^ 2 - t233, t62 * t129 + t153, 0, 0, 0, 0, 0, 0, -t129 * t76 + t152, -t129 * t78 + t145, t241 * t133 + t242 * t136, -t43 * t129 + t12 + (t157 * t190 + t1) * t133 + t240 * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, -t76 ^ 2 + t78 ^ 2, -t242, -t224, t241, -t66, -t43 * t78 + t240, -t157 * t238 + t43 * t76 - t1, 0, 0;];
tauc_reg = t6;
