% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR15
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR15_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:30
% EndTime: 2019-12-31 18:37:36
% DurationCPUTime: 3.11s
% Computational Cost: add. (3442->412), mult. (6883->543), div. (0->0), fcn. (4465->10), ass. (0->209)
t144 = sin(qJ(5));
t147 = cos(qJ(5));
t141 = sin(pkin(8));
t148 = cos(qJ(3));
t218 = qJD(1) * t148;
t193 = t141 * t218;
t142 = cos(pkin(8));
t212 = t142 * qJD(3);
t94 = -t193 + t212;
t192 = t142 * t218;
t213 = t141 * qJD(3);
t95 = t192 + t213;
t39 = t144 * t95 - t147 * t94;
t270 = t39 ^ 2;
t42 = t144 * t94 + t147 * t95;
t269 = t42 ^ 2;
t145 = sin(qJ(3));
t211 = t145 * qJD(1);
t119 = qJD(5) + t211;
t268 = t39 * t119;
t98 = t144 * t141 - t147 * t142;
t74 = t98 * t145;
t149 = cos(qJ(1));
t135 = g(2) * t149;
t146 = sin(qJ(1));
t136 = g(1) * t146;
t266 = t136 - t135;
t267 = t266 * t142;
t150 = -pkin(1) - pkin(6);
t113 = t150 * qJDD(1) + qJDD(2);
t139 = t145 ^ 2;
t140 = t148 ^ 2;
t220 = t139 + t140;
t184 = t220 * t113;
t214 = qJD(5) * t147;
t215 = qJD(5) * t144;
t265 = -t141 * t215 + t142 * t214;
t260 = pkin(4) * t141;
t264 = -t260 + t150;
t221 = t139 - t140;
t180 = qJD(1) * t221;
t259 = g(3) * t145;
t263 = t266 * t148 - t259;
t169 = t145 * pkin(3) - t148 * qJ(4);
t103 = qJ(2) + t169;
t87 = t103 * qJD(1);
t114 = t150 * qJD(1) + qJD(2);
t102 = t145 * t114;
t88 = qJD(3) * qJ(4) + t102;
t34 = -t141 * t88 + t142 * t87;
t18 = pkin(4) * t211 - t95 * pkin(7) + t34;
t35 = t141 * t87 + t142 * t88;
t20 = t94 * pkin(7) + t35;
t164 = t144 * t20 - t147 * t18;
t170 = pkin(3) * t148 + qJ(4) * t145;
t76 = qJD(3) * t170 - t148 * qJD(4) + qJD(2);
t32 = qJD(1) * t76 + qJDD(1) * t103;
t233 = t148 * t114;
t45 = qJDD(3) * qJ(4) + t145 * t113 + (qJD(4) + t233) * qJD(3);
t14 = -t141 * t45 + t142 * t32;
t209 = qJD(1) * qJD(3);
t190 = t148 * t209;
t206 = t145 * qJDD(1);
t156 = t190 + t206;
t191 = t145 * t209;
t205 = t148 * qJDD(1);
t116 = t142 * t205;
t226 = t141 * qJDD(3) + t116;
t59 = t142 * t191 - t226;
t8 = t156 * pkin(4) + t59 * pkin(7) + t14;
t15 = t141 * t32 + t142 * t45;
t115 = t141 * t205;
t177 = -t142 * qJDD(3) + t115;
t58 = t141 * t191 - t177;
t9 = t58 * pkin(7) + t15;
t1 = -qJD(5) * t164 + t144 * t8 + t147 * t9;
t262 = 0.2e1 * qJ(2);
t261 = t58 * pkin(4);
t258 = g(3) * t148;
t257 = t142 * pkin(7);
t256 = t42 * t39;
t255 = pkin(7) + qJ(4);
t204 = t145 * t257;
t101 = t170 * qJD(1);
t48 = t142 * t101 - t141 * t233;
t24 = (pkin(4) * t148 + t204) * qJD(1) + t48;
t200 = t141 * t211;
t49 = t141 * t101 + t142 * t233;
t33 = pkin(7) * t200 + t49;
t109 = t255 * t141;
t110 = t255 * t142;
t51 = -t147 * t109 - t144 * t110;
t254 = -qJD(4) * t98 + qJD(5) * t51 - t144 * t24 - t147 * t33;
t52 = -t144 * t109 + t147 * t110;
t99 = t147 * t141 + t144 * t142;
t253 = -qJD(4) * t99 - qJD(5) * t52 + t144 * t33 - t147 * t24;
t72 = t99 * t145;
t75 = t98 * t148;
t82 = t99 * qJD(1);
t252 = -qJD(3) * t75 - qJD(5) * t72 - t82;
t216 = qJD(3) * t148;
t251 = t98 * qJD(1) + qJD(5) * t74 - t99 * t216;
t84 = t99 * qJD(5);
t250 = t145 * t82 + t84;
t199 = t142 * t211;
t249 = -t144 * t200 + t147 * t199 + t265;
t248 = t148 * t58;
t247 = t148 * t59;
t217 = qJD(3) * t145;
t163 = -qJDD(3) * pkin(3) + t114 * t217 + qJDD(4);
t50 = -t148 * t113 + t163;
t246 = t50 * t141;
t245 = t50 * t142;
t244 = t50 * t148;
t243 = t58 * t142;
t242 = t59 * t141;
t196 = t150 * t216;
t44 = t141 * t76 + t142 * t196;
t237 = t145 * t150;
t57 = t141 * t103 + t142 * t237;
t241 = pkin(1) * qJDD(1);
t152 = qJD(1) ^ 2;
t240 = t141 * t152;
t239 = t142 * t152;
t238 = t145 * t149;
t137 = pkin(8) + qJ(5);
t128 = sin(t137);
t236 = t146 * t128;
t129 = cos(t137);
t235 = t146 * t129;
t234 = t146 * t148;
t232 = t148 * t149;
t231 = t149 * t128;
t230 = t149 * t129;
t151 = qJD(3) ^ 2;
t229 = t150 * t151;
t228 = t152 * qJ(2);
t79 = -qJD(3) * pkin(3) + qJD(4) - t233;
t227 = -qJD(4) + t79;
t225 = g(1) * t232 + g(2) * t234;
t203 = 0.2e1 * qJD(1) * qJD(2);
t224 = (qJDD(1) * qJ(2) + t203) * qJ(2);
t223 = t149 * pkin(1) + t146 * qJ(2);
t219 = -t151 - t152;
t210 = qJ(4) * qJDD(1);
t208 = qJDD(3) * t145;
t207 = t139 * qJDD(1);
t202 = t148 * t152 * t145;
t201 = t149 * pkin(6) + t223;
t198 = t144 * t217;
t197 = t147 * t217;
t189 = -t141 * t150 + pkin(4);
t187 = -g(2) * t238 + t258;
t186 = -g(1) * t234 + t259;
t185 = -t144 * t59 - t147 * t58;
t183 = -t113 + t228;
t182 = -t95 + t213;
t181 = -t94 + t212;
t179 = t220 * qJDD(1);
t178 = qJDD(2) - t241;
t176 = g(2) * t201;
t175 = t145 * t190;
t174 = g(1) * t149 + g(2) * t146;
t171 = -t266 - t228;
t168 = -t14 * t141 + t15 * t142;
t167 = t141 * t35 + t142 * t34;
t166 = -t34 * t141 + t35 * t142;
t165 = t141 * t95 - t142 * t94;
t6 = t144 * t18 + t147 * t20;
t92 = t142 * t103;
t37 = t189 * t145 - t148 * t257 + t92;
t47 = -t141 * t148 * pkin(7) + t57;
t16 = -t144 * t47 + t147 * t37;
t17 = t144 * t37 + t147 * t47;
t125 = t142 * pkin(4) + pkin(3);
t161 = t145 * t125 - t148 * t255;
t160 = qJDD(1) * t262 + t203;
t159 = t266 * t141;
t12 = -t144 * t58 + t147 * t59 - t94 * t214 + t95 * t215;
t158 = -t145 * t136 - t187;
t157 = qJDD(3) * t150 + t209 * t262;
t2 = -qJD(5) * t6 - t144 * t9 + t147 * t8;
t154 = t160 - t174;
t13 = qJD(5) * t42 + t185;
t153 = (-t113 - t135) * t148 + t163 - t186;
t132 = t149 * qJ(2);
t130 = qJDD(3) * t148;
t97 = qJDD(5) + t156;
t93 = t264 * t148;
t89 = 0.2e1 * t175 + t207;
t80 = t264 * t217;
t73 = t99 * t148;
t68 = t145 * t230 - t236;
t67 = t145 * t231 + t235;
t66 = t145 * t235 + t231;
t65 = -t145 * t236 + t230;
t64 = -pkin(4) * t200 + t102;
t61 = t142 * t76;
t56 = -t141 * t237 + t92;
t46 = -t94 * pkin(4) + t79;
t43 = -t141 * t196 + t61;
t31 = t145 * pkin(7) * t213 + t44;
t30 = -t141 * t197 - t142 * t198 + t265 * t148;
t28 = -t141 * t198 + t142 * t197 + t148 * t84;
t21 = t61 + (t189 * t148 + t204) * qJD(3);
t19 = t50 - t261;
t4 = -qJD(5) * t17 - t144 * t31 + t147 * t21;
t3 = qJD(5) * t16 + t144 * t21 + t147 * t31;
t5 = [0, 0, 0, 0, 0, qJDD(1), t266, t174, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t266 - 0.2e1 * t241, t154, -t178 * pkin(1) - g(1) * (-t146 * pkin(1) + t132) - g(2) * t223 + t224, t140 * qJDD(1) - 0.2e1 * t175, 0.2e1 * qJD(3) * t180 - 0.2e1 * t145 * t205, -t151 * t145 + t130, t89, -t151 * t148 - t208, 0, t157 * t148 + (t154 - t229) * t145, -t157 * t145 + (t160 - t229) * t148 - t225, -t150 * t179 - t184 + t266, -g(1) * (t150 * t146 + t132) - t176 + t150 * t184 + t224, (-t95 * t217 - t247) * t142, (t242 + t243) * t148 + t165 * t217, (-t59 + t116) * t145 + (-t142 * t180 + t148 * t95) * qJD(3), (t94 * t217 - t248) * t141, (t58 - t115) * t145 + (t141 * t180 + t148 * t94) * qJD(3), t89, t159 + (t246 + t150 * t58 + (qJD(1) * t56 + t34) * qJD(3)) * t148 + (t43 * qJD(1) + t56 * qJDD(1) + t14 - t174 * t142 + (-t141 * t79 - t150 * t94) * qJD(3)) * t145, t267 + (t245 + t150 * t59 + (-qJD(1) * t57 - t35) * qJD(3)) * t148 + (-t44 * qJD(1) - t57 * qJDD(1) - t15 + t174 * t141 + (-t142 * t79 + t150 * t95) * qJD(3)) * t145, -t43 * t95 + t44 * t94 + t56 * t59 + t57 * t58 + (-t14 * t142 - t141 * t15) * t148 + t167 * t217 + t225, t15 * t57 + t35 * t44 + t14 * t56 + t34 * t43 - g(1) * (pkin(3) * t238 - qJ(4) * t232 + t132) - t176 + (t79 * t217 - t244) * t150 + (-g(1) * t150 - g(2) * t169) * t146, t12 * t75 - t42 * t28, t12 * t73 + t75 * t13 + t28 * t39 - t42 * t30, -t28 * t119 - t12 * t145 + t42 * t216 - t75 * t97, t13 * t73 + t39 * t30, -t30 * t119 - t13 * t145 - t39 * t216 - t73 * t97, t119 * t216 + t97 * t145, -g(1) * t68 - g(2) * t66 + t4 * t119 - t93 * t13 + t2 * t145 + t16 * t97 - t164 * t216 + t19 * t73 + t46 * t30 + t80 * t39, g(1) * t67 - g(2) * t65 - t1 * t145 - t3 * t119 + t93 * t12 - t17 * t97 - t19 * t75 - t216 * t6 - t46 * t28 + t80 * t42, -t1 * t73 + t16 * t12 - t17 * t13 - t164 * t28 + t2 * t75 - t3 * t39 - t6 * t30 - t4 * t42 + t225, t1 * t17 + t6 * t3 + t2 * t16 - t164 * t4 - t19 * t93 + t46 * t80 - g(1) * (t125 * t238 - t232 * t255 + t132) - g(2) * (t149 * t260 + t201) + (-g(1) * t264 - g(2) * t161) * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t152, t171 + t178, 0, 0, 0, 0, 0, 0, t219 * t145 + t130, t219 * t148 - t208, -t179, t184 + t171, 0, 0, 0, 0, 0, 0, -t141 * t207 + t248 + (-t239 + (-t94 - 0.2e1 * t193) * qJD(3)) * t145, -t142 * t207 + t247 + (t240 + (t95 - 0.2e1 * t192) * qJD(3)) * t145, (qJD(1) * t95 + t145 * t58 + t94 * t216) * t142 + (-qJD(1) * t94 - t145 * t59 + t95 * t216) * t141, -t244 + t168 * t145 - t167 * qJD(1) + (t145 * t79 + t148 * t166) * qJD(3) - t266, 0, 0, 0, 0, 0, 0, t119 * t251 - t148 * t13 + t217 * t39 - t72 * t97, -t119 * t252 + t148 * t12 + t217 * t42 + t74 * t97, -t72 * t12 + t74 * t13 - t251 * t42 - t252 * t39, -t1 * t74 - t19 * t148 - t164 * t251 - t2 * t72 + t217 * t46 + t252 * t6 - t266; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, -t221 * t152, t205, -t202, -t206, qJDD(3), (-t183 + t135) * t148 + t186, (t183 + t136) * t145 + t187, 0, 0, t199 * t95 - t242, t141 * t58 - t59 * t142 - t165 * t211, t139 * t239 + t141 * t206 + t182 * t218, -t200 * t94 + t243, -t139 * t240 + t142 * t206 + t181 * t218, -t202, pkin(3) * t58 - t245 + (-t267 + (-qJ(4) * t213 - t34) * qJD(1)) * t148 + (-t141 * t210 + g(3) * t142 + t114 * t94 + (t227 * t141 - t48) * qJD(1)) * t145, pkin(3) * t59 + t246 + (t159 + (-qJ(4) * t212 + t35) * qJD(1)) * t148 + (-t142 * t210 - g(3) * t141 - t114 * t95 + (t227 * t142 + t49) * qJD(1)) * t145, t48 * t95 - t49 * t94 + (qJ(4) * t58 + qJD(4) * t94 - t34 * t211 + t15) * t142 + (-qJ(4) * t59 + qJD(4) * t95 - t35 * t211 - t14) * t141 + t158, -t79 * t102 - t34 * t48 - t35 * t49 + t166 * qJD(4) + (-t50 - t263) * pkin(3) + (-t145 * t266 + t168 - t258) * qJ(4), -t12 * t99 + t249 * t42, t12 * t98 - t99 * t13 - t249 * t39 - t250 * t42, t249 * t119 - t42 * t218 + t99 * t97, t13 * t98 + t250 * t39, -t250 * t119 + t39 * t218 - t98 * t97, -t119 * t218, t119 * t253 - t125 * t13 - t129 * t263 + t164 * t218 + t19 * t98 + t250 * t46 - t64 * t39 + t51 * t97, -t254 * t119 + t125 * t12 + t128 * t263 + t19 * t99 + t6 * t218 + t249 * t46 - t64 * t42 - t52 * t97, -t1 * t98 + t51 * t12 - t52 * t13 + t164 * t249 - t2 * t99 - t250 * t6 - t253 * t42 - t254 * t39 + t158, g(3) * t161 + t1 * t52 - t19 * t125 + t2 * t51 - t253 * t164 + t254 * t6 - t46 * t64 - t266 * (t125 * t148 + t145 * t255); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182 * t211 + t177, -t181 * t211 + t226, -t94 ^ 2 - t95 ^ 2, t34 * t95 - t35 * t94 + t153, 0, 0, 0, 0, 0, 0, t42 * t119 + t13, -t12 - t268, -t269 - t270, -t164 * t42 + t39 * t6 + t153 - t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256, t269 - t270, -t12 + t268, -t256, -t185 + (-qJD(5) + t119) * t42, t97, -g(1) * t65 - g(2) * t67 + t6 * t119 + t128 * t258 - t46 * t42 + t2, g(1) * t66 - g(2) * t68 - t119 * t164 + t129 * t258 + t46 * t39 - t1, 0, 0;];
tau_reg = t5;
