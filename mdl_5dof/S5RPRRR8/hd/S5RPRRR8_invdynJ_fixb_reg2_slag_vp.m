% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:11
% EndTime: 2019-12-31 19:06:15
% DurationCPUTime: 2.10s
% Computational Cost: add. (3707->299), mult. (5474->368), div. (0->0), fcn. (3192->10), ass. (0->176)
t127 = sin(qJ(3));
t133 = qJD(4) ^ 2;
t119 = qJDD(1) - qJDD(3);
t130 = cos(qJ(3));
t206 = t130 * t119;
t193 = qJD(1) - qJD(3);
t248 = t193 ^ 2;
t249 = (t133 + t248) * t127 + t206;
t125 = sin(qJ(5));
t128 = cos(qJ(5));
t129 = cos(qJ(4));
t207 = t129 * t119;
t126 = sin(qJ(4));
t211 = t126 * t119;
t120 = qJD(4) + qJD(5);
t209 = t128 * t126;
t69 = t125 * t129 + t209;
t247 = t120 * t69;
t16 = -t125 * t211 + t128 * t207 - t193 * t247;
t228 = t193 * pkin(3);
t200 = qJ(2) * qJD(1);
t132 = -pkin(1) - pkin(2);
t91 = t132 * qJD(1) + qJD(2);
t60 = -t127 * t200 + t130 * t91;
t43 = -t60 + t228;
t246 = t193 * t43;
t217 = t193 * t60;
t61 = t127 * t91 + t130 * t200;
t216 = t61 * t193;
t245 = qJD(3) * t200 - t132 * qJDD(1) - qJDD(2);
t240 = t130 * t193;
t196 = qJD(5) * t128;
t198 = qJD(4) * t129;
t243 = t128 * t198 + t129 * t196;
t232 = sin(qJ(1));
t233 = cos(qJ(1));
t68 = -t232 * t127 - t233 * t130;
t70 = t233 * t127 - t232 * t130;
t170 = g(1) * t68 + g(2) * t70;
t242 = qJD(4) * t193;
t194 = qJD(1) * qJD(2);
t195 = qJ(2) * qJDD(1);
t241 = qJD(3) * t91 + t194 + t195;
t208 = t128 * t129;
t212 = t125 * t126;
t67 = -t208 + t212;
t63 = t67 * t127;
t78 = -t127 * qJ(2) + t130 * t132;
t79 = t130 * qJ(2) + t127 * t132;
t122 = t126 ^ 2;
t123 = t129 ^ 2;
t201 = t122 + t123;
t124 = qJ(4) + qJ(5);
t107 = sin(t124);
t199 = qJD(4) * t126;
t153 = -t193 * t199 + t207;
t173 = t241 * t127 + t245 * t130;
t230 = t119 * pkin(3);
t23 = t173 + t230;
t14 = t153 * pkin(4) + t23;
t171 = g(1) * t70 - g(2) * t68;
t227 = t129 * pkin(4);
t105 = pkin(3) + t227;
t36 = t105 * t193 - t60;
t161 = t120 * t212;
t37 = t161 - t243;
t239 = -t107 * t171 + t14 * t69 - t36 * t37;
t108 = cos(t124);
t238 = t108 * t171 + t14 * t67 + t247 * t36;
t197 = qJD(5) * t125;
t44 = -pkin(7) * t193 + t61;
t183 = -pkin(8) * t193 + t44;
t34 = t183 * t126;
t33 = qJD(4) * pkin(4) - t34;
t35 = t183 * t129;
t186 = t193 * t198;
t229 = t119 * pkin(7);
t28 = -t245 * t127 + t241 * t130;
t22 = t28 - t229;
t6 = -t44 * t198 + qJDD(4) * pkin(4) - t126 * t22 + (t186 + t211) * pkin(8);
t7 = -t153 * pkin(8) + t129 * t22 - t44 * t199;
t1 = (qJD(5) * t33 + t7) * t128 + t125 * t6 - t35 * t197;
t237 = pkin(8) + pkin(7);
t72 = -pkin(7) + t79;
t234 = pkin(8) - t72;
t231 = g(3) * t129;
t190 = t193 * t212;
t49 = t193 * t208 - t190;
t51 = t69 * t193;
t226 = t51 * t49;
t225 = -t127 * t247 + t67 * t240;
t224 = t120 * t63 + t69 * t240;
t88 = t237 * t126;
t89 = t237 * t129;
t40 = -t125 * t89 - t128 * t88;
t187 = qJD(4) * t237;
t73 = t126 * t187;
t74 = t129 * t187;
t223 = t40 * qJD(5) - t125 * t74 - t128 * t73 + t67 * t60;
t41 = -t125 * t88 + t128 * t89;
t222 = -t41 * qJD(5) + t125 * t73 - t128 * t74 + t69 * t60;
t221 = t125 * t35;
t220 = t128 * t35;
t58 = t130 * qJD(2) + t78 * qJD(3);
t219 = t58 * t193;
t59 = t127 * qJD(2) + t79 * qJD(3);
t218 = t59 * t193;
t215 = pkin(1) * qJDD(1);
t213 = t193 * t126;
t210 = t127 * t119;
t204 = t233 * pkin(1) + t232 * qJ(2);
t203 = g(1) * t232 - g(2) * t233;
t202 = t122 - t123;
t192 = pkin(4) * t199;
t191 = 0.2e1 * t194;
t189 = t126 * t248 * t129;
t188 = t233 * pkin(2) + t204;
t181 = qJD(4) * t234;
t178 = t201 * t22;
t177 = t119 * t209 + t125 * t207 + t243 * t193;
t176 = t127 * t193;
t174 = qJDD(2) - t215;
t172 = -t61 + t192;
t71 = pkin(3) - t78;
t169 = -t232 * pkin(1) + t233 * qJ(2);
t15 = -t161 * t193 + t177;
t165 = t15 * t69 - t51 * t37;
t164 = t16 * t67 + t247 * t49;
t162 = -0.2e1 * t126 * t186;
t160 = -t70 * t105 - t237 * t68;
t159 = t68 * t105 - t237 * t70;
t118 = qJDD(4) + qJDD(5);
t158 = t69 * t118 - t37 * t120;
t157 = t67 * t118 + t120 * t247;
t11 = t125 * t33 + t220;
t47 = t234 * t126;
t48 = t234 * t129;
t26 = t125 * t48 + t128 * t47;
t27 = t125 * t47 - t128 * t48;
t156 = t171 - t23;
t150 = g(1) * t233 + g(2) * t232;
t149 = -t170 - t22 + t246;
t148 = -t232 * pkin(2) + t169;
t147 = t193 * t201;
t146 = t171 - t173;
t145 = -pkin(7) * qJDD(4) + (t43 + t60 + t228) * qJD(4);
t144 = -qJDD(4) * t72 + (-t193 * t71 - t43 - t58) * qJD(4);
t143 = t15 * t67 - t69 * t16 + t247 * t51 + t37 * t49;
t142 = -qJDD(4) * t127 + 0.2e1 * t130 * t242;
t2 = -t11 * qJD(5) - t125 * t7 + t128 * t6;
t141 = pkin(7) * t133 - t156 + t216 + t230;
t140 = -t119 * t71 + t133 * t72 + t156 - t218;
t139 = -t170 - t28;
t10 = t128 * t33 - t221;
t138 = -t1 * t67 + t10 * t37 - t11 * t247 - t2 * t69 + t170;
t136 = -g(3) * t107 - t170 * t108 + t36 * t49 - t1;
t135 = g(3) * t108 - t170 * t107 + t36 * t51 + t2;
t134 = qJD(1) ^ 2;
t83 = qJDD(4) * t129 - t133 * t126;
t82 = qJDD(4) * t126 + t133 * t129;
t64 = t71 + t227;
t62 = t69 * t127;
t53 = t123 * t119 + t162;
t52 = -t122 * t119 + t162;
t42 = t59 - t192;
t39 = -t126 * t207 + t202 * t242;
t31 = -t126 * t58 + t129 * t181;
t30 = t126 * t181 + t129 * t58;
t17 = -t49 ^ 2 + t51 ^ 2;
t13 = -t128 * t34 - t221;
t12 = t125 * t34 - t220;
t9 = -t51 * t120 - t16;
t8 = -t177 + (t190 + t49) * t120;
t4 = -t27 * qJD(5) - t125 * t30 + t128 * t31;
t3 = t26 * qJD(5) + t125 * t31 + t128 * t30;
t5 = [0, 0, 0, 0, 0, qJDD(1), t203, t150, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + t203 + 0.2e1 * t215, 0, -t150 + t191 + 0.2e1 * t195, -t174 * pkin(1) - g(1) * t169 - g(2) * t204 + (t191 + t195) * qJ(2), 0, 0, 0, 0, 0, t119, -t78 * t119 - t146 + t218, t79 * t119 - t139 + t219, 0, -g(1) * t148 - g(2) * t188 - t173 * t78 + t28 * t79 + t61 * t58 - t60 * t59, -t52, -0.2e1 * t39, -t82, t53, -t83, 0, t144 * t126 - t140 * t129, t140 * t126 + t144 * t129, -t170 + t201 * (-t119 * t72 - t219 - t22), t23 * t71 + t43 * t59 - g(1) * (t70 * pkin(3) + t68 * pkin(7) + t148) - g(2) * (-t68 * pkin(3) + t70 * pkin(7) + t188) + t201 * (t22 * t72 + t44 * t58), t165, -t143, -t158, -t164, t157, 0, t26 * t118 + t4 * t120 + t64 * t16 + t42 * t49 - t238, -t27 * t118 - t3 * t120 - t64 * t15 - t42 * t51 - t239, t26 * t15 - t27 * t16 - t3 * t49 + t4 * t51 - t138, t1 * t27 + t11 * t3 + t2 * t26 + t10 * t4 + t14 * t64 + t36 * t42 - g(1) * (t148 - t160) - g(2) * (-t159 + t188); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t134, -t134 * qJ(2) + t174 - t203, 0, 0, 0, 0, 0, 0, -t127 * t248 - t206, -t130 * t248 + t210, 0, (-t173 - t216) * t130 + (t28 + t217) * t127 - t203, 0, 0, 0, 0, 0, 0, t142 * t126 - t249 * t129, t249 * t126 + t142 * t129, t147 * t240 - t201 * t210, (t178 - t246) * t127 + (-t147 * t44 - t23) * t130 - t203, 0, 0, 0, 0, 0, 0, -t62 * t118 + t120 * t224 - t130 * t16 - t176 * t49, t63 * t118 - t120 * t225 + t130 * t15 + t176 * t51, -t62 * t15 + t63 * t16 + t224 * t51 - t225 * t49, -t1 * t63 + t10 * t224 + t11 * t225 - t14 * t130 - t176 * t36 - t2 * t62 - t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, t146 - t216, t139 - t217, 0, 0, t52, 0.2e1 * t39, t82, -t53, t83, 0, t145 * t126 - t141 * t129, t141 * t126 + t129 * t145, t170 + t201 * (t217 + t22 - t229), -t43 * t61 - t201 * t60 * t44 + t156 * pkin(3) + (t178 + t170) * pkin(7), -t165, t143, t158, t164, -t157, 0, -t105 * t16 + t40 * t118 + t120 * t222 + t172 * t49 + t238, t105 * t15 - t41 * t118 - t120 * t223 - t172 * t51 + t239, t40 * t15 - t41 * t16 + t222 * t51 - t223 * t49 + t138, -g(1) * t160 - g(2) * t159 + t1 * t41 + t10 * t222 - t14 * t105 + t11 * t223 + t172 * t36 + t2 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, t202 * t248, -t211, t189, -t207, qJDD(4), t149 * t126 + t231, -g(3) * t126 + t129 * t149, 0, 0, -t226, t17, t8, t226, t9, t118, -t12 * t120 + (t118 * t128 - t120 * t197 + t213 * t49) * pkin(4) + t135, t13 * t120 + (-t118 * t125 - t120 * t196 - t213 * t51) * pkin(4) + t136, -(t11 + t12) * t51 + (-t10 + t13) * t49 + (-t125 * t16 + t128 * t15 + (-t125 * t51 - t128 * t49) * qJD(5)) * pkin(4), -t10 * t12 - t11 * t13 + (t231 + t1 * t125 + t128 * t2 + (-t10 * t125 + t11 * t128) * qJD(5) + (t193 * t36 - t170) * t126) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t226, t17, t8, t226, t9, t118, t11 * t120 + t135, t10 * t120 + t136, 0, 0;];
tau_reg = t5;
