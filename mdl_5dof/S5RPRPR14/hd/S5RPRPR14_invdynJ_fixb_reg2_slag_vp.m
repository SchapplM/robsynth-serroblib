% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR14_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:21
% EndTime: 2019-12-31 18:35:26
% DurationCPUTime: 2.30s
% Computational Cost: add. (3499->362), mult. (6889->455), div. (0->0), fcn. (4538->10), ass. (0->186)
t120 = sin(pkin(8));
t123 = sin(qJ(3));
t126 = cos(qJ(3));
t210 = cos(pkin(8));
t138 = -t120 * t126 - t210 * t123;
t65 = t138 * qJD(1);
t234 = qJD(5) - t65;
t122 = sin(qJ(5));
t125 = cos(qJ(5));
t128 = -pkin(1) - pkin(6);
t81 = t128 * qJD(1) + qJD(2);
t166 = -qJ(4) * qJD(1) + t81;
t57 = t166 * t123;
t52 = t210 * t57;
t193 = qJD(1) * t126;
t58 = -qJ(4) * t193 + t126 * t81;
t54 = qJD(3) * pkin(3) + t58;
t27 = t120 * t54 + t52;
t23 = qJD(3) * pkin(7) + t27;
t171 = t210 * t126;
t158 = qJD(1) * t171;
t194 = qJD(1) * t123;
t176 = t120 * t194;
t68 = t158 - t176;
t77 = pkin(3) * t194 + qJD(1) * qJ(2) + qJD(4);
t28 = -t65 * pkin(4) - t68 * pkin(7) + t77;
t8 = t122 * t28 + t125 * t23;
t240 = t234 * t8;
t7 = -t122 * t23 + t125 * t28;
t239 = t7 * t234;
t130 = qJD(1) ^ 2;
t124 = sin(qJ(1));
t127 = cos(qJ(1));
t233 = g(1) * t124 - g(2) * t127;
t238 = -t130 * qJ(2) - t233;
t168 = t125 * t234;
t187 = qJD(1) * qJD(3);
t175 = t123 * t187;
t165 = qJDD(1) * t210;
t183 = t126 * qJDD(1);
t182 = qJD(3) * t158 + t120 * t183 + t123 * t165;
t41 = t120 * t175 - t182;
t39 = -qJDD(5) + t41;
t214 = t122 * t39;
t237 = t168 * t234 - t214;
t115 = qJDD(1) * qJ(2);
t118 = t123 ^ 2;
t119 = t126 ^ 2;
t196 = t118 + t119;
t80 = t128 * qJDD(1) + qJDD(2);
t172 = t196 * t80;
t113 = qJ(3) + pkin(8);
t102 = cos(t113);
t156 = g(1) * t127 + g(2) * t124;
t235 = t156 * t102;
t101 = sin(t113);
t133 = g(3) * t101 - t102 * t233;
t186 = qJD(1) * qJD(4);
t192 = qJD(3) * t123;
t75 = t126 * t80;
t25 = -t126 * t186 - t81 * t192 + qJDD(3) * pkin(3) + t75 + (t175 - t183) * qJ(4);
t191 = qJD(3) * t126;
t31 = t166 * t191 + (-qJ(4) * qJDD(1) - t186 + t80) * t123;
t9 = -t120 * t31 + t210 * t25;
t5 = -qJDD(3) * pkin(4) - t9;
t92 = t120 * pkin(3) + pkin(7);
t232 = -qJD(5) * t234 * t92 + t133 - t5;
t116 = qJD(1) * qJD(2);
t180 = 0.2e1 * t116;
t231 = 0.2e1 * t115 + t180 - t156;
t225 = g(3) * t102;
t230 = -t101 * t233 - t225;
t229 = t68 ^ 2;
t67 = -qJD(3) * t171 + t120 * t192;
t227 = t7 * t67;
t224 = g(3) * t123;
t109 = t123 * pkin(3);
t188 = t125 * qJD(3);
t46 = t122 * t68 - t188;
t223 = t46 * t65;
t48 = t122 * qJD(3) + t125 * t68;
t222 = t48 * t46;
t221 = t48 * t68;
t220 = t68 * t46;
t219 = t68 * t65;
t189 = qJD(5) * t125;
t184 = t123 * qJDD(1);
t145 = -t120 * t184 + t126 * t165;
t70 = t138 * qJD(3);
t42 = -qJD(1) * t70 - t145;
t170 = -t125 * qJDD(3) - t122 * t42;
t19 = qJD(5) * t48 + t170;
t218 = -t122 * t19 - t46 * t189;
t10 = t120 * t25 + t210 * t31;
t74 = -t120 * t123 + t171;
t217 = t70 * qJD(3) + t74 * qJDD(3);
t216 = qJD(5) * t8;
t215 = t120 * t57;
t213 = t122 * t48;
t36 = t125 * t39;
t190 = qJD(5) * t122;
t18 = -qJD(5) * t188 - t122 * qJDD(3) + t125 * t42 + t190 * t68;
t212 = t18 * t122;
t211 = t19 * t125;
t209 = pkin(1) * qJDD(1);
t208 = t124 * t122;
t207 = t124 * t125;
t206 = t127 * t122;
t205 = t127 * t125;
t203 = t77 * qJD(1);
t94 = qJ(2) + t109;
t202 = qJ(4) - t128;
t85 = pkin(3) * t191 + qJD(2);
t200 = (t180 + t115) * qJ(2);
t199 = t127 * pkin(1) + t124 * qJ(2);
t197 = t118 - t119;
t129 = qJD(3) ^ 2;
t195 = -t129 - t130;
t185 = qJDD(3) * t123;
t179 = t74 * t190;
t178 = t74 * t189;
t177 = t126 * t130 * t123;
t174 = t126 * t187;
t108 = t127 * qJ(2);
t173 = -t124 * pkin(1) + t108;
t169 = t122 * t234;
t167 = t202 * t126;
t164 = -qJD(5) * t138 + qJD(1);
t163 = t196 * qJDD(1);
t162 = qJDD(2) - t209;
t161 = t123 * t174;
t51 = qJDD(4) + t115 + t116 + (t174 + t184) * pkin(3);
t12 = -t41 * pkin(4) + t42 * pkin(7) + t51;
t6 = qJDD(3) * pkin(7) + t10;
t1 = qJD(5) * t7 + t122 * t12 + t125 * t6;
t160 = -t1 * t138 - t8 * t67;
t26 = t210 * t54 - t215;
t22 = -qJD(3) * pkin(4) - t26;
t159 = t22 * t70 + t5 * t74;
t157 = t101 * pkin(4) - t102 * pkin(7);
t154 = t122 * t8 + t125 * t7;
t153 = t138 * t18 - t48 * t67;
t152 = -t74 * t18 + t70 * t48;
t151 = t138 * t19 + t46 * t67;
t150 = t74 * t19 + t70 * t46;
t149 = -t234 * t70 + t39 * t74;
t148 = t138 * t41 + t65 * t67;
t147 = -t74 * t42 + t68 * t70;
t40 = -pkin(4) * t138 - t74 * pkin(7) + t94;
t78 = t202 * t123;
t44 = -t120 * t167 - t210 * t78;
t16 = -t122 * t44 + t125 * t40;
t17 = t122 * t40 + t125 * t44;
t121 = -qJ(4) - pkin(6);
t144 = t127 * t109 + t124 * t121 + t173;
t143 = t67 * qJD(3) + qJDD(3) * t138;
t142 = -t36 + (t122 * t65 - t190) * t234;
t141 = t124 * t109 - t127 * t121 + t199;
t140 = -qJD(5) * t28 + t225 - t6;
t139 = t22 * t234 + t39 * t92;
t137 = 0.2e1 * qJ(2) * t187 + qJDD(3) * t128;
t134 = -t126 * qJD(4) + t202 * t192;
t132 = -t10 * t138 + t26 * t70 - t27 * t67 + t9 * t74 - t233;
t131 = -t128 * t129 + t231;
t105 = qJDD(3) * t126;
t93 = -t210 * pkin(3) - pkin(4);
t64 = t65 ^ 2;
t63 = t101 * t205 - t208;
t62 = t101 * t206 + t207;
t61 = t101 * t207 + t206;
t60 = -t101 * t208 + t205;
t55 = -qJD(3) * t167 - t123 * qJD(4);
t43 = -t120 * t78 + t210 * t167;
t35 = pkin(3) * t193 + t68 * pkin(4) - t65 * pkin(7);
t34 = -t67 * pkin(4) - t70 * pkin(7) + t85;
t33 = t210 * t58 - t215;
t32 = t120 * t58 + t52;
t30 = t120 * t134 + t210 * t55;
t29 = t120 * t55 - t210 * t134;
t14 = t122 * t35 + t125 * t33;
t13 = -t122 * t33 + t125 * t35;
t11 = t125 * t12;
t4 = -qJD(5) * t17 - t122 * t30 + t125 * t34;
t3 = qJD(5) * t16 + t122 * t34 + t125 * t30;
t2 = -t122 * t6 + t11 - t216;
t15 = [0, 0, 0, 0, 0, qJDD(1), t233, t156, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t233 - 0.2e1 * t209, t231, -t162 * pkin(1) - g(1) * t173 - g(2) * t199 + t200, t119 * qJDD(1) - 0.2e1 * t161, -0.2e1 * t123 * t183 + 0.2e1 * t187 * t197, -t129 * t123 + t105, t118 * qJDD(1) + 0.2e1 * t161, -t129 * t126 - t185, 0, t123 * t131 + t126 * t137, -t123 * t137 + t126 * t131, -t128 * t163 - t172 + t233, -g(1) * (t128 * t124 + t108) - g(2) * (t127 * pkin(6) + t199) + t128 * t172 + t200, t147, -t138 * t42 + t74 * t41 + t70 * t65 + t68 * t67, t217, t148, t143, 0, -t29 * qJD(3) - t43 * qJDD(3) - t101 * t156 - t138 * t51 - t94 * t41 - t85 * t65 - t77 * t67, -t30 * qJD(3) - t44 * qJDD(3) - t94 * t42 + t51 * t74 + t85 * t68 + t77 * t70 - t235, t29 * t68 + t30 * t65 + t44 * t41 - t43 * t42 - t132, -g(1) * t144 - g(2) * t141 + t10 * t44 - t26 * t29 + t27 * t30 - t9 * t43 + t51 * t94 + t77 * t85, t125 * t152 - t179 * t48, -(t125 * t46 + t213) * t70 + (t212 - t211 + (t122 * t46 - t125 * t48) * qJD(5)) * t74, -t125 * t149 - t179 * t234 + t153, t122 * t150 + t178 * t46, t122 * t149 - t178 * t234 + t151, t138 * t39 - t234 * t67, -g(1) * t63 - g(2) * t61 + t122 * t159 - t138 * t2 - t16 * t39 + t178 * t22 + t43 * t19 + t234 * t4 + t29 * t46 - t227, g(1) * t62 - g(2) * t60 + t125 * t159 + t17 * t39 - t179 * t22 - t43 * t18 - t234 * t3 + t29 * t48 - t160, t16 * t18 - t17 * t19 - t3 * t46 - t4 * t48 - t154 * t70 + t235 + (-t1 * t122 - t2 * t125 + (t122 * t7 - t125 * t8) * qJD(5)) * t74, t1 * t17 + t8 * t3 + t2 * t16 + t7 * t4 + t5 * t43 + t22 * t29 - g(1) * (t127 * t157 + t144) - g(2) * (t124 * t157 + t141); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t130, t238 + t162, 0, 0, 0, 0, 0, 0, t123 * t195 + t105, t126 * t195 - t185, -t163, t172 + t238, 0, 0, 0, 0, 0, 0, qJD(1) * t65 + t217, -qJD(1) * t68 + t143, -t147 - t148, t132 - t203, 0, 0, 0, 0, 0, 0, -t138 * t214 + (t122 * t67 - t125 * t164) * t234 - t150, -t138 * t36 + (t122 * t164 + t125 * t67) * t234 - t152, (t164 * t48 + t151) * t125 + (t164 * t46 + t153) * t122, (-t164 * t7 + t160) * t125 + (-qJD(1) * t8 + t227 - (-t2 - t216) * t138) * t122 - t159 - t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, -t197 * t130, t183, -t177, -t184, qJDD(3), t126 * t238 + t224 + t75, g(3) * t126 + (-t238 - t80) * t123, 0, 0, -t219, -t64 + t229, t145, t219, (t68 + t176) * qJD(3) - t182, qJDD(3), t32 * qJD(3) - t77 * t68 + (t210 * qJDD(3) + t65 * t193) * pkin(3) + t133 + t9, t33 * qJD(3) - t77 * t65 + (-qJDD(3) * t120 - t193 * t68) * pkin(3) - t10 - t230, (t27 - t32) * t68 - (-t26 + t33) * t65 + (t120 * t41 + t210 * t42) * pkin(3), t26 * t32 - t27 * t33 + (t210 * t9 + t224 + t10 * t120 + (-t233 - t203) * t126) * pkin(3), t168 * t48 - t212, (-t18 + t223) * t125 - t234 * t213 + t218, -t221 + t237, t169 * t46 - t211, t142 + t220, -t234 * t68, t139 * t122 + t232 * t125 - t13 * t234 + t93 * t19 - t32 * t46 - t7 * t68, -t232 * t122 + t139 * t125 + t14 * t234 - t93 * t18 - t32 * t48 + t8 * t68, t13 * t48 + t14 * t46 + (-t19 * t92 + t65 * t7 + t1 + (t48 * t92 - t7) * qJD(5)) * t125 + (-t18 * t92 + t65 * t8 - t2 + (t46 * t92 - t8) * qJD(5)) * t122 + t230, t5 * t93 - t8 * t14 - t7 * t13 - t22 * t32 - g(3) * (-t109 - t157) + (-qJD(5) * t154 + t1 * t125 - t2 * t122) * t92 - t233 * (pkin(3) * t126 + pkin(4) * t102 + pkin(7) * t101); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t68 - t176) * qJD(3) + t182, 0.2e1 * t65 * qJD(3) + t145, -t64 - t229, t26 * t68 - t27 * t65 - t156 + t51, 0, 0, 0, 0, 0, 0, t142 - t220, -t221 - t237, (t18 + t223) * t125 + t48 * t169 + t218, -t22 * t68 + (t2 + t240) * t125 + (t1 - t239) * t122 - t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, -t46 ^ 2 + t48 ^ 2, t234 * t46 - t18, -t222, -t170 + (-qJD(5) + t234) * t48, -t39, -g(1) * t60 - g(2) * t62 + t122 * t140 - t189 * t23 - t22 * t48 + t11 + t240, g(1) * t61 - g(2) * t63 + t22 * t46 + t239 + (qJD(5) * t23 - t12) * t122 + t140 * t125, 0, 0;];
tau_reg = t15;
