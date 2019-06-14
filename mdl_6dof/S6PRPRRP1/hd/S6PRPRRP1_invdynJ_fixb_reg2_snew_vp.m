% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 23:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:32:32
% EndTime: 2019-05-04 23:32:37
% DurationCPUTime: 2.13s
% Computational Cost: add. (7287->273), mult. (13717->377), div. (0->0), fcn. (9909->12), ass. (0->187)
t180 = cos(qJ(4));
t169 = sin(pkin(10));
t172 = cos(pkin(10));
t146 = g(1) * t169 - g(2) * t172;
t170 = sin(pkin(6));
t165 = -g(3) + qJDD(1);
t173 = cos(pkin(6));
t195 = t165 * t173 + qJDD(3);
t188 = -t146 * t170 + t195;
t114 = t180 * t188;
t177 = sin(qJ(4));
t182 = qJD(2) ^ 2;
t230 = pkin(4) * t180;
t194 = -pkin(9) * t177 - t230;
t168 = sin(pkin(11));
t171 = cos(pkin(11));
t147 = -g(1) * t172 - g(2) * t169;
t178 = sin(qJ(2));
t181 = cos(qJ(2));
t215 = t173 * t146;
t242 = t165 * t170 + t215;
t101 = -t147 * t178 + t242 * t181;
t187 = qJDD(2) * pkin(2) + t101;
t102 = t181 * t147 + t242 * t178;
t95 = -t182 * pkin(2) + t102;
t67 = t168 * t187 + t171 * t95;
t63 = -pkin(3) * t182 + qJDD(2) * pkin(8) + t67;
t197 = t182 * t194 + t63;
t238 = qJD(4) ^ 2;
t40 = -qJDD(4) * pkin(4) - t238 * pkin(9) + t197 * t177 - t114;
t176 = sin(qJ(5));
t179 = cos(qJ(5));
t211 = qJD(2) * t177;
t134 = -t179 * qJD(4) + t176 * t211;
t136 = qJD(4) * t176 + t179 * t211;
t113 = t136 * t134;
t207 = qJD(2) * qJD(4);
t158 = t177 * t207;
t205 = t180 * qJDD(2);
t140 = -t158 + t205;
t133 = -qJDD(5) + t140;
t239 = -t113 - t133;
t245 = pkin(5) * t239;
t199 = t180 * t207;
t206 = t177 * qJDD(2);
t139 = t199 + t206;
t107 = -qJD(5) * t134 + qJDD(4) * t176 + t139 * t179;
t155 = qJD(2) * t180 - qJD(5);
t124 = t134 * t155;
t88 = t107 - t124;
t244 = qJ(6) * t88;
t220 = t239 * t176;
t219 = t239 * t179;
t132 = t136 ^ 2;
t153 = t155 ^ 2;
t110 = -t132 - t153;
t131 = t134 ^ 2;
t196 = -t179 * qJDD(4) + t176 * t139;
t106 = -qJD(5) * t136 - t196;
t119 = -pkin(5) * t155 - qJ(6) * t136;
t185 = t177 * t188;
t41 = -t238 * pkin(4) + qJDD(4) * pkin(9) + t197 * t180 + t185;
t192 = -t140 + t158;
t193 = t139 + t199;
t198 = t168 * t95 - t171 * t187;
t62 = -qJDD(2) * pkin(3) - t182 * pkin(8) + t198;
t48 = t192 * pkin(4) - t193 * pkin(9) + t62;
t23 = t176 * t48 + t179 * t41;
t189 = t106 * qJ(6) - 0.2e1 * qJD(6) * t134 + t155 * t119 + t23;
t243 = -t189 + (t110 + t131) * pkin(5);
t29 = -t106 * pkin(5) - t131 * qJ(6) + t136 * t119 + qJDD(6) + t40;
t240 = t107 + t124;
t84 = (qJD(5) + t155) * t136 + t196;
t108 = -t153 - t131;
t70 = t108 * t176 + t219;
t237 = pkin(4) * t70;
t98 = -t113 + t133;
t223 = t176 * t98;
t74 = t110 * t179 + t223;
t236 = pkin(4) * t74;
t209 = qJD(6) * t136;
t127 = -0.2e1 * t209;
t22 = t176 * t41 - t179 * t48;
t186 = -t22 - t244 + t245;
t17 = t127 + t186;
t235 = pkin(5) * t17;
t234 = pkin(5) * t88;
t59 = -t176 * t84 - t179 * t88;
t233 = pkin(9) * t59;
t232 = pkin(9) * t70;
t231 = pkin(9) * t74;
t60 = t176 * t88 - t179 * t84;
t97 = -t131 - t132;
t229 = -pkin(4) * t97 + pkin(9) * t60;
t71 = t108 * t179 - t220;
t83 = (qJD(5) - t155) * t136 + t196;
t228 = -pkin(4) * t83 + pkin(9) * t71;
t221 = t179 * t98;
t75 = -t110 * t176 + t221;
t227 = -pkin(4) * t240 + pkin(9) * t75;
t226 = t17 * t176;
t225 = t17 * t179;
t224 = t176 * t40;
t222 = t179 * t40;
t218 = t155 * t176;
t217 = t155 * t179;
t154 = t177 * t182 * t180;
t148 = qJDD(4) + t154;
t214 = t177 * t148;
t149 = qJDD(4) - t154;
t213 = t180 * t149;
t37 = t177 * t97 + t180 * t60;
t25 = t168 * t37 - t171 * t59;
t204 = pkin(2) * t25 - pkin(3) * t59 + pkin(8) * t37;
t44 = t177 * t83 + t180 * t71;
t31 = t168 * t44 - t171 * t70;
t203 = pkin(2) * t31 - pkin(3) * t70 + pkin(8) * t44;
t50 = t177 * t240 + t180 * t75;
t33 = t168 * t50 - t171 * t74;
t202 = pkin(2) * t33 - pkin(3) * t74 + pkin(8) * t50;
t201 = t180 * t113;
t13 = t176 * t22 + t179 * t23;
t54 = t177 * t63 - t114;
t55 = t180 * t63 + t185;
t28 = t177 * t54 + t180 * t55;
t12 = t176 * t23 - t179 * t22;
t183 = t186 + t245;
t164 = t180 ^ 2;
t163 = t177 ^ 2;
t162 = t164 * t182;
t160 = t163 * t182;
t152 = -t162 - t238;
t151 = -t160 - t238;
t145 = t160 + t162;
t144 = (t163 + t164) * qJDD(2);
t143 = -qJDD(2) * t168 - t171 * t182;
t142 = qJDD(2) * t171 - t168 * t182;
t141 = -0.2e1 * t158 + t205;
t138 = 0.2e1 * t199 + t206;
t128 = 0.2e1 * t209;
t121 = -t132 + t153;
t120 = t131 - t153;
t118 = -t151 * t177 - t213;
t117 = t152 * t180 - t214;
t116 = -t149 * t177 + t151 * t180;
t115 = t148 * t180 + t152 * t177;
t111 = t132 - t131;
t109 = t144 * t168 + t145 * t171;
t94 = t118 * t168 - t138 * t171;
t93 = t117 * t168 + t141 * t171;
t89 = (t134 * t176 + t136 * t179) * t155;
t80 = t107 * t176 - t136 * t217;
t79 = t106 * t179 - t134 * t218;
t78 = t180 * t133 + t177 * (t134 * t179 - t136 * t176) * t155;
t77 = t120 * t176 - t221;
t76 = t121 * t179 + t220;
t65 = t177 * (t107 * t179 + t136 * t218) - t201;
t64 = t177 * (-t106 * t176 - t134 * t217) + t201;
t61 = -pkin(5) * t240 + qJ(6) * t98;
t58 = -t176 * t83 + t179 * t240;
t53 = t177 * (t120 * t179 + t223) + t180 * t84;
t52 = t177 * (-t121 * t176 + t219) - t180 * t88;
t49 = t177 * t75 - t180 * t240;
t43 = t177 * t71 - t180 * t83;
t39 = t177 * (-t176 * t240 - t179 * t83) - t180 * t111;
t36 = t177 * t60 - t180 * t97;
t34 = t168 * t67 - t171 * t198;
t27 = t177 * t55 - t180 * t54;
t26 = -qJ(6) * t110 + t29;
t20 = -pkin(5) * t83 + qJ(6) * t108 - t29;
t19 = t168 * t28 - t171 * t62;
t18 = -pkin(5) * t131 + t189;
t16 = t128 - t186 + t244;
t15 = -qJ(6) * t84 + (-t131 - t97) * pkin(5) + t189;
t14 = t173 * t49 + (t178 * (t168 * t74 + t171 * t50) + t181 * t33) * t170;
t11 = t173 * t43 + (t178 * (t168 * t70 + t171 * t44) + t181 * t31) * t170;
t10 = -pkin(5) * t29 + qJ(6) * t18;
t9 = t13 * t180 + t177 * t40;
t8 = t13 * t177 - t180 * t40;
t7 = t173 * t36 + (t178 * (t168 * t59 + t171 * t37) + t181 * t25) * t170;
t6 = t179 * t18 - t226;
t5 = t176 * t18 + t225;
t4 = t177 * t29 + t180 * t6;
t3 = t177 * t6 - t180 * t29;
t2 = -t12 * t171 + t168 * t9;
t1 = t168 * t4 - t171 * t5;
t21 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t165, 0, 0, 0, 0, 0, 0, (qJDD(2) * t181 - t178 * t182) * t170, (-qJDD(2) * t178 - t181 * t182) * t170, 0, t173 ^ 2 * t165 + (t101 * t181 + t102 * t178 - t215) * t170, 0, 0, 0, 0, 0, 0, (t142 * t181 + t143 * t178) * t170, (-t142 * t178 + t143 * t181) * t170, 0, t173 * t195 + (t178 * (t168 * t198 + t171 * t67) + t181 * t34 - t215) * t170, 0, 0, 0, 0, 0, 0, t173 * t115 + (t178 * (t117 * t171 - t141 * t168) + t181 * t93) * t170, t173 * t116 + (t178 * (t118 * t171 + t138 * t168) + t181 * t94) * t170, (t178 * (t144 * t171 - t145 * t168) + t181 * t109) * t170, t173 * t27 + (t178 * (t168 * t62 + t171 * t28) + t181 * t19) * t170, 0, 0, 0, 0, 0, 0, t11, t14, t7, t173 * t8 + (t178 * (t12 * t168 + t171 * t9) + t181 * t2) * t170, 0, 0, 0, 0, 0, 0, t11, t14, t7, t173 * t3 + (t178 * (t168 * t5 + t171 * t4) + t181 * t1) * t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t101, -t102, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t142 - t198, pkin(2) * t143 - t67, 0, pkin(2) * t34, t193 * t177, t138 * t180 + t141 * t177, t214 + t180 * (-t160 + t238), -t192 * t180, t177 * (t162 - t238) + t213, 0, pkin(2) * t93 + pkin(3) * t141 + pkin(8) * t117 - t180 * t62, pkin(2) * t94 - pkin(3) * t138 + pkin(8) * t118 + t177 * t62, pkin(2) * t109 + pkin(3) * t145 + pkin(8) * t144 + t28, pkin(2) * t19 - pkin(3) * t62 + pkin(8) * t28, t65, t39, t52, t64, t53, t78, t177 * (t224 - t232) + t180 * (t22 - t237) + t203, t177 * (t222 - t231) + t180 * (t23 - t236) + t202, t177 * (-t12 - t233) - t59 * t230 + t204, pkin(2) * t2 + pkin(8) * t9 + (-pkin(3) + t194) * t12, t65, t39, t52, t64, t53, t78, t177 * (-qJ(6) * t219 - t176 * t20 - t232) + t180 * (t128 - t183 - t237) + t203, t177 * (-t176 * t61 + t179 * t26 - t231) + t180 * (-t236 - t243) + t202, t177 * (-t15 * t176 + t16 * t179 - t233) + t180 * (-pkin(4) * t59 + t234) + t204, t177 * (-pkin(9) * t5 - qJ(6) * t225 - t10 * t176) + t180 * (-pkin(4) * t5 - t235) - pkin(3) * t5 + pkin(8) * t4 + pkin(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, 0, 0, 0, 0, 0, 0, t115, t116, 0, t27, 0, 0, 0, 0, 0, 0, t43, t49, t36, t8, 0, 0, 0, 0, 0, 0, t43, t49, t36, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t160 - t162, t206, t154, t205, qJDD(4), -t54, -t55, 0, 0, t80, t58, t76, t79, t77, t89, -t222 + t228, t224 + t227, t13 + t229, -pkin(4) * t40 + pkin(9) * t13, t80, t58, t76, t79, t77, t89, -qJ(6) * t220 + t179 * t20 + t228, t176 * t26 + t179 * t61 + t227, t15 * t179 + t16 * t176 + t229, -pkin(4) * t29 + pkin(9) * t6 - qJ(6) * t226 + t10 * t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, t111, t88, -t113, -t84, -t133, -t22, -t23, 0, 0, t113, t111, t88, -t113, -t84, -t133, t127 + t183, t243, -t234, t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t240, t97, t29;];
tauJ_reg  = t21;
