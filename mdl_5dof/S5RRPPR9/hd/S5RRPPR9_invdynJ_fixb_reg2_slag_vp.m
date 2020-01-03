% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPR9
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:53
% EndTime: 2019-12-31 19:41:58
% DurationCPUTime: 2.82s
% Computational Cost: add. (2218->401), mult. (4466->469), div. (0->0), fcn. (2380->6), ass. (0->215)
t125 = sin(qJ(5));
t126 = sin(qJ(2));
t128 = cos(qJ(5));
t129 = cos(qJ(2));
t199 = qJD(5) * t129;
t186 = t125 * t199;
t198 = t128 * qJD(2);
t97 = t129 * qJDD(1);
t18 = -t128 * t97 - qJD(5) * t198 - t125 * qJDD(2) + (t126 * t198 + t186) * qJD(1);
t206 = qJD(1) * t129;
t48 = t125 * t206 - t198;
t207 = qJD(1) * t126;
t72 = qJD(5) + t207;
t240 = t48 * t72;
t252 = t18 - t240;
t196 = qJD(1) * qJD(2);
t183 = t129 * t196;
t95 = t126 * qJDD(1);
t254 = t183 + t95;
t184 = t126 * t196;
t205 = qJD(2) * t125;
t49 = t128 * t206 + t205;
t19 = t49 * qJD(5) - t128 * qJDD(2) + (-t184 + t97) * t125;
t239 = t49 * t72;
t146 = t19 - t239;
t247 = pkin(2) + pkin(3);
t253 = t126 * t247;
t187 = t247 * qJD(2);
t185 = t247 * qJDD(2);
t238 = pkin(6) - qJ(4);
t119 = qJD(2) * qJ(3);
t91 = pkin(6) * t206;
t52 = -qJ(4) * t206 + t91;
t38 = -t119 - t52;
t105 = t126 * pkin(4);
t169 = t129 * pkin(7) + t105;
t127 = sin(qJ(1));
t113 = g(2) * t127;
t130 = cos(qJ(1));
t250 = g(1) * t130 + t113;
t204 = qJD(2) * t126;
t149 = pkin(6) * t204 + qJD(4) * t129;
t116 = qJDD(2) * qJ(3);
t117 = qJD(2) * qJD(3);
t88 = pkin(6) * t97;
t190 = t116 + t117 + t88;
t65 = qJ(4) * t184;
t15 = qJ(4) * t97 + t149 * qJD(1) - t190 - t65;
t13 = qJDD(2) * pkin(4) - t15;
t197 = pkin(7) + t247;
t242 = g(3) * t126;
t249 = -qJD(5) * t197 * t72 + t129 * t250 - t13 + t242;
t227 = pkin(6) * qJDD(2);
t39 = -qJD(1) * pkin(1) - pkin(2) * t206 - qJ(3) * t207;
t101 = t126 * qJ(3);
t109 = t129 * pkin(2);
t211 = t109 + t101;
t54 = -pkin(1) - t211;
t248 = (qJD(1) * t54 + t39) * qJD(2) - t227;
t28 = pkin(3) * t206 + qJD(4) - t39;
t22 = t169 * qJD(1) + t28;
t194 = pkin(6) * t207;
t167 = qJD(3) + t194;
t82 = qJ(4) * t207;
t156 = t167 - t82;
t26 = -t197 * qJD(2) + t156;
t7 = -t125 * t26 + t128 * t22;
t246 = t7 * t72;
t8 = t125 * t22 + t128 * t26;
t245 = t72 * t8;
t114 = g(1) * t127;
t244 = g(2) * qJ(4);
t243 = g(2) * t130;
t112 = g(3) * t129;
t108 = t129 * pkin(3);
t241 = t48 * t49;
t203 = qJD(2) * t129;
t99 = t126 * qJD(3);
t237 = qJ(3) * t203 + t99;
t236 = qJD(2) * pkin(2);
t235 = t125 * t18;
t47 = qJDD(5) + t254;
t234 = t125 * t47;
t157 = t125 * t72;
t233 = t128 * t19;
t232 = t128 * t47;
t231 = t128 * t48;
t230 = t128 * t72;
t229 = t129 * t48;
t228 = t129 * t49;
t226 = qJ(3) * t129;
t225 = qJD(2) * t38;
t224 = qJD(2) * t48;
t223 = qJD(2) * t49;
t122 = qJDD(1) * pkin(1);
t222 = qJDD(2) * pkin(2);
t121 = t129 ^ 2;
t133 = qJD(1) ^ 2;
t221 = t121 * t133;
t220 = t126 * t127;
t219 = t126 * t130;
t218 = t126 * t133;
t217 = t127 * t128;
t216 = t127 * t129;
t215 = t128 * t130;
t214 = t129 * t130;
t50 = -t82 + t194;
t213 = qJD(3) + t50;
t212 = -qJD(4) - t28;
t210 = t130 * pkin(1) + t127 * pkin(6);
t120 = t126 ^ 2;
t208 = t120 + t121;
t202 = qJD(4) * t126;
t201 = qJD(5) * t125;
t200 = qJD(5) * t128;
t195 = g(1) * t219 + g(2) * t220 - t112;
t193 = t126 * t157;
t192 = t126 * t230;
t191 = 0.2e1 * t116 + 0.2e1 * t117 + t88;
t71 = pkin(6) * t183;
t87 = pkin(6) * t95;
t189 = qJDD(3) + t71 + t87;
t188 = t108 + t211;
t182 = -pkin(1) - t101;
t181 = t114 - t243;
t180 = g(1) * t197;
t44 = pkin(1) + t188;
t179 = qJD(1) * t44 + t28;
t177 = t212 * t126;
t176 = -t87 + t195;
t175 = pkin(2) * t214 + qJ(3) * t219 + t210;
t174 = t126 * t187;
t173 = t126 * t183;
t172 = pkin(2) * t97 + t254 * qJ(3) + qJD(1) * t99 + t122;
t171 = t208 * qJDD(1) * pkin(6);
t170 = pkin(3) * t214 + t175;
t132 = qJD(2) ^ 2;
t168 = pkin(6) * t132 + t243;
t165 = -t125 * t8 - t128 * t7;
t164 = -t125 * t7 + t128 * t8;
t110 = t130 * pkin(6);
t163 = g(1) * (-qJ(4) * t130 + t110);
t162 = -qJDD(3) + t176;
t160 = t188 + t169;
t25 = pkin(1) + t160;
t60 = t238 * t126;
t16 = -t125 * t60 + t128 * t25;
t17 = t125 * t25 + t128 * t60;
t53 = t167 - t236;
t58 = t91 + t119;
t161 = t126 * t58 - t129 * t53;
t155 = t182 - t109;
t154 = t71 - t162;
t135 = -qJ(4) * t254 - qJD(1) * t202 + t189;
t12 = -t197 * qJDD(2) + t135;
t153 = -qJD(5) * t22 - t112 - t12;
t152 = -t72 * t200 - t234;
t151 = t72 * t201 - t232;
t150 = -0.2e1 * pkin(1) * t196 - t227;
t148 = pkin(4) * t129 - t197 * t126;
t147 = pkin(3) * t97 + qJDD(4) + t172;
t144 = -t168 + 0.2e1 * t122;
t142 = t148 * qJD(2);
t32 = qJD(2) * pkin(4) - t38;
t141 = t197 * t47 - t72 * t32;
t27 = -t174 + t237;
t9 = -qJD(1) * t174 + t147;
t140 = -qJD(1) * t27 - qJDD(1) * t44 + t243 - t9;
t139 = -qJ(4) * t95 + t154;
t20 = pkin(2) * t184 - t172;
t37 = pkin(2) * t204 - t237;
t138 = -qJD(1) * t37 - qJDD(1) * t54 - t168 - t20;
t29 = -pkin(6) * t184 + t190;
t33 = t189 - t222;
t137 = -t161 * qJD(2) + t33 * t126 + t29 * t129;
t6 = qJD(1) * t142 + t169 * qJDD(1) + t147;
t1 = qJD(5) * t7 + t12 * t128 + t125 * t6;
t5 = t128 * t6;
t2 = -qJD(5) * t8 - t12 * t125 + t5;
t136 = t164 * qJD(5) + t1 * t125 + t128 * t2 - t243;
t124 = qJ(3) + pkin(4);
t98 = t120 * t133;
t84 = qJ(3) * t206;
t77 = g(1) * t216;
t76 = g(1) * t220;
t70 = qJ(3) * t214;
t68 = qJ(3) * t216;
t67 = t129 * t218;
t63 = -t98 - t132;
t61 = t238 * t129;
t59 = qJDD(2) + t67;
t57 = -t98 + t221;
t56 = qJDD(2) * t129 - t126 * t132;
t55 = qJDD(2) * t126 + t129 * t132;
t51 = pkin(2) * t207 - t84;
t46 = qJDD(1) * t121 - 0.2e1 * t173;
t45 = qJDD(1) * t120 + 0.2e1 * t173;
t43 = -t125 * t127 + t126 * t215;
t42 = -t125 * t219 - t217;
t41 = -t125 * t130 - t126 * t217;
t40 = t125 * t220 - t215;
t36 = t238 * t203 - t202;
t35 = -qJ(4) * t204 + t149;
t34 = -t247 * t207 + t84;
t31 = -t187 + t156;
t30 = t126 * t97 + (-t120 + t121) * t196;
t24 = 0.2e1 * t30;
t23 = t148 * qJD(1) + t84;
t21 = t142 + t237;
t14 = -t185 + t135;
t11 = t125 * t23 + t128 * t52;
t10 = -t125 * t52 + t128 * t23;
t4 = -t17 * qJD(5) - t125 * t36 + t128 * t21;
t3 = t16 * qJD(5) + t125 * t21 + t128 * t36;
t62 = [0, 0, 0, 0, 0, qJDD(1), t181, t250, 0, 0, t45, t24, t55, t46, t56, 0, t150 * t126 + t144 * t129 + t77, -t144 * t126 + t150 * t129 - t76, 0.2e1 * t171 - t250, -g(1) * (-pkin(1) * t127 + t110) - g(2) * t210 + (t208 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), t45, t55, -0.2e1 * t30, 0, -t56, t46, t248 * t126 + t138 * t129 + t77, t171 + t137 - t250, t138 * t126 - t248 * t129 + t76, pkin(6) * t137 - g(1) * t110 - g(2) * t175 - t155 * t114 + t20 * t54 + t39 * t37, t46, t24, t56, t45, t55, 0, qJDD(2) * t61 + t76 + (t179 * t129 - t35) * qJD(2) - t140 * t126, qJDD(2) * t60 - t77 + (t179 * t126 + t36) * qJD(2) + t140 * t129, (-qJD(2) * t31 - qJDD(1) * t61 + t15 + (-qJD(2) * t60 + t35) * qJD(1)) * t129 + (-t225 - qJDD(1) * t60 - t14 + (qJD(2) * t61 - t36) * qJD(1)) * t126 + t250, t14 * t60 + t31 * t36 - t15 * t61 + t38 * t35 + t9 * t44 + t28 * t27 - t163 - g(2) * t170 + (-g(1) * (t155 - t108) + t244) * t127, -t49 * t186 + (-t129 * t18 - t49 * t204) * t128, (t125 * t49 + t231) * t204 + (t235 - t233 + (t125 * t48 - t128 * t49) * qJD(5)) * t129, (t72 * t198 + t18) * t126 + (t151 - t223) * t129, t199 * t231 + (t129 * t19 - t48 * t204) * t125, (-t72 * t205 + t19) * t126 + (-t152 + t224) * t129, t126 * t47 + t72 * t203, -g(1) * t41 - g(2) * t43 + t16 * t47 - t19 * t61 + t35 * t48 + t4 * t72 + (t32 * t205 + t2) * t126 + (qJD(2) * t7 - t125 * t13 - t32 * t200) * t129, -g(1) * t40 - g(2) * t42 - t17 * t47 + t18 * t61 - t3 * t72 + t35 * t49 + (t198 * t32 - t1) * t126 + (-qJD(2) * t8 - t128 * t13 + t201 * t32) * t129, t129 * t136 - t16 * t18 + t165 * t204 + t17 * t19 + t3 * t48 + t4 * t49 + t77, t1 * t17 + t8 * t3 + t2 * t16 + t7 * t4 + t13 * t61 - t32 * t35 - t163 - g(2) * (pkin(4) * t219 + pkin(7) * t214 + t170) + (-g(1) * (t182 - t105) + t244 + t129 * t180) * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t57, t95, t67, t97, qJDD(2), pkin(1) * t218 + t176, t242 - t88 + (pkin(1) * t133 + t250) * t129, 0, 0, -t67, t95, t57, qJDD(2), -t97, t67, 0.2e1 * t222 + (-t126 * t39 + t129 * t51) * qJD(1) + t162, (-pkin(2) * t126 + t226) * qJDD(1) + ((t58 - t119) * t126 + (qJD(3) - t53 - t236) * t129) * qJD(1), (qJD(1) * t51 - g(3)) * t126 + (qJD(1) * t39 - t250) * t129 + t191, t29 * qJ(3) + t58 * qJD(3) - t33 * pkin(2) - t39 * t51 - g(1) * (-pkin(2) * t219 + t70) - g(2) * (-pkin(2) * t220 + t68) - g(3) * t211 + t161 * qJD(1) * pkin(6), t67, -t57, t97, -t67, t95, qJDD(2), qJD(2) * t50 + t65 + (-g(3) + (-pkin(6) * qJD(2) - t34) * qJD(1)) * t126 + (-qJ(4) * qJDD(1) + t212 * qJD(1) - t250) * t129 + t191, -qJD(2) * t52 - 0.2e1 * t185 + ((-qJ(4) * qJD(2) + t34) * t129 + t177) * qJD(1) + t139, (-t226 + t253) * qJDD(1) + (-t213 + t31 + t187) * t206, -g(1) * t70 - g(2) * t68 - g(3) * t188 - t15 * qJ(3) - t14 * t247 - t213 * t38 + t250 * t253 - t28 * t34 - t31 * t52, t49 * t230 - t235, (-t18 - t240) * t128 + (-t19 - t239) * t125, (-t192 + t228) * qJD(1) + t152, t157 * t48 - t233, (t193 - t229) * qJD(1) + t151, -t72 * t206, -t10 * t72 - t124 * t19 + t141 * t125 - t249 * t128 - t7 * t206 - t213 * t48, t11 * t72 + t124 * t18 + t249 * t125 + t141 * t128 + t8 * t206 - t213 * t49, -t10 * t49 - t11 * t48 + (t7 * t207 - t197 * t19 - t1 + (t197 * t49 + t7) * qJD(5)) * t128 + (t8 * t207 - t197 * t18 + t2 + (t197 * t48 + t8) * qJD(5)) * t125 + t195, t13 * t124 - t8 * t11 - t7 * t10 - g(1) * (pkin(4) * t214 + t70) - g(2) * (pkin(4) * t216 + t68) - g(3) * t160 + t213 * t32 + (t197 * t113 + t130 * t180) * t126 - (qJD(5) * t165 + t1 * t128 - t125 * t2) * t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t95, t63, -qJD(2) * t58 + t39 * t207 + t154 - t222, 0, 0, 0, 0, 0, 0, t63, t59, -t95, t225 - t185 + (-qJ(4) * t203 + t177) * qJD(1) + t139, 0, 0, 0, 0, 0, 0, -t128 * t72 ^ 2 + t224 - t234, t157 * t72 + t223 - t232, t252 * t125 + t146 * t128, -qJD(2) * t32 + (t1 - t246) * t128 + (-t2 - t245) * t125 - t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 + 0.2e1 * t183, -t97 + 0.2e1 * t184, -t98 - t221, (-t129 * t38 + (t31 - t187) * t126) * qJD(1) + t147 + t181, 0, 0, 0, 0, 0, 0, (-t193 - t229) * qJD(1) - t151, (-t192 - t228) * qJD(1) + t152, t146 * t125 - t252 * t128, t114 + (t126 * t164 + t129 * t32) * qJD(1) + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, -t48 ^ 2 + t49 ^ 2, t252, -t241, t146, t47, -g(1) * t42 + g(2) * t40 + t125 * t153 - t200 * t26 + t32 * t49 + t245 + t5, g(1) * t43 - g(2) * t41 - t32 * t48 + t246 + (qJD(5) * t26 - t6) * t125 + t153 * t128, 0, 0;];
tau_reg = t62;
