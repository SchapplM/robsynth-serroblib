% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:22:49
% EndTime: 2019-03-09 03:22:54
% DurationCPUTime: 2.31s
% Computational Cost: add. (3830->365), mult. (7594->462), div. (0->0), fcn. (5116->10), ass. (0->188)
t234 = sin(pkin(9));
t190 = qJD(3) * t234;
t235 = cos(pkin(9));
t267 = qJD(1) * t190 - qJDD(1) * t235;
t191 = qJD(3) * t235;
t266 = qJD(1) * t191 + qJDD(1) * t234;
t156 = qJD(1) ^ 2;
t150 = sin(qJ(1));
t153 = cos(qJ(1));
t261 = g(1) * t150 - g(2) * t153;
t164 = -qJ(2) * t156 - t261;
t151 = cos(qJ(5));
t149 = sin(qJ(3));
t152 = cos(qJ(3));
t162 = t235 * t149 + t234 * t152;
t86 = t162 * qJD(1);
t263 = qJD(5) + t86;
t189 = t151 * t263;
t148 = sin(qJ(5));
t59 = t149 * t267 - t266 * t152;
t57 = -qJDD(5) + t59;
t238 = t148 * t57;
t265 = t263 * t189 - t238;
t136 = t149 * pkin(3);
t147 = -qJ(4) - pkin(7);
t262 = t153 * t136 + t150 * t147;
t94 = -t234 * t149 + t235 * t152;
t140 = qJ(3) + pkin(9);
t128 = sin(t140);
t225 = t151 * t153;
t228 = t148 * t150;
t82 = -t128 * t228 + t225;
t226 = t150 * t151;
t227 = t148 * t153;
t84 = t128 * t227 + t226;
t260 = -g(1) * t82 - g(2) * t84;
t119 = t234 * pkin(3) + pkin(8);
t154 = -pkin(1) - pkin(7);
t105 = t154 * qJD(1) + qJD(2);
t204 = t152 * qJDD(1);
t207 = qJD(1) * qJD(4);
t208 = qJD(1) * qJD(3);
t214 = qJD(3) * t149;
t104 = t154 * qJDD(1) + qJDD(2);
t96 = t152 * t104;
t35 = -t152 * t207 - t105 * t214 + qJDD(3) * pkin(3) + t96 + (t149 * t208 - t204) * qJ(4);
t213 = qJD(3) * t152;
t44 = (-qJ(4) * qJD(1) + t105) * t213 + (-qJ(4) * qJDD(1) + t104 - t207) * t149;
t17 = -t234 * t44 + t235 * t35;
t13 = -qJDD(3) * pkin(4) - t17;
t129 = cos(t140);
t160 = -g(3) * t128 + t129 * t261;
t259 = qJD(5) * t119 * t263 + t13 + t160;
t141 = qJDD(1) * qJ(2);
t180 = g(1) * t153 + g(2) * t150;
t142 = qJD(1) * qJD(2);
t202 = 0.2e1 * t142;
t258 = 0.2e1 * t141 + t202 - t180;
t18 = t234 * t35 + t235 * t44;
t14 = qJDD(3) * pkin(8) + t18;
t216 = qJD(1) * t149;
t79 = -qJ(4) * t216 + t105 * t149;
t197 = t235 * t79;
t215 = qJD(1) * t152;
t80 = -qJ(4) * t215 + t152 * t105;
t76 = qJD(3) * pkin(3) + t80;
t38 = t234 * t76 + t197;
t30 = qJD(3) * pkin(8) + t38;
t89 = t94 * qJD(1);
t98 = pkin(3) * t216 + qJD(1) * qJ(2) + qJD(4);
t39 = pkin(4) * t86 - pkin(8) * t89 + t98;
t16 = t148 * t39 + t151 * t30;
t161 = t266 * t149 + t152 * t267;
t196 = t152 * t208;
t205 = t149 * qJDD(1);
t169 = qJDD(4) + t141 + t142 + (t196 + t205) * pkin(3);
t21 = -t59 * pkin(4) + t161 * pkin(8) + t169;
t20 = t151 * t21;
t210 = t151 * qJD(3);
t212 = qJD(5) * t148;
t25 = -qJD(5) * t210 - t148 * qJDD(3) + t151 * t161 + t89 * t212;
t67 = qJD(3) * t148 + t151 * t89;
t1 = -t57 * pkin(5) + t25 * qJ(6) - t16 * qJD(5) - t67 * qJD(6) - t148 * t14 + t20;
t65 = t148 * t89 - t210;
t9 = -qJ(6) * t65 + t16;
t257 = t263 * t9 + t1;
t256 = t67 ^ 2;
t15 = -t148 * t30 + t151 * t39;
t8 = -qJ(6) * t67 + t15;
t6 = pkin(5) * t263 + t8;
t254 = -t8 + t6;
t222 = qJ(6) + t119;
t187 = qJD(5) * t222;
t230 = qJ(6) * t151;
t50 = pkin(3) * t215 + pkin(4) * t89 + pkin(8) * t86;
t46 = t151 * t50;
t73 = t234 * t79;
t48 = t235 * t80 - t73;
t251 = -pkin(5) * t89 - t151 * t187 - t86 * t230 - t46 + (-qJD(6) + t48) * t148;
t250 = pkin(5) * t148;
t248 = g(3) * t129;
t247 = g(3) * t149;
t246 = t151 * pkin(5);
t245 = t65 * t86;
t244 = t65 * t89;
t243 = t67 * t89;
t211 = qJD(5) * t151;
t158 = -t151 * qJDD(3) - t148 * t161;
t26 = t67 * qJD(5) + t158;
t242 = -t148 * t26 - t65 * t211;
t241 = t148 * t50 + t151 * t48;
t224 = qJ(2) + t136;
t58 = pkin(4) * t162 - pkin(8) * t94 + t224;
t223 = qJ(4) - t154;
t100 = t223 * t152;
t99 = t223 * t149;
t62 = -t234 * t100 - t235 * t99;
t60 = t151 * t62;
t240 = t148 * t58 + t60;
t231 = qJ(6) * t148;
t239 = t151 * qJD(6) - t148 * t187 - t86 * t231 - t241;
t237 = t148 * t67;
t52 = t151 * t57;
t236 = t25 * t148;
t233 = pkin(1) * qJDD(1);
t229 = qJD(1) * t98;
t220 = t153 * pkin(1) + t150 * qJ(2);
t145 = t152 ^ 2;
t218 = t149 ^ 2 - t145;
t155 = qJD(3) ^ 2;
t217 = -t155 - t156;
t209 = pkin(3) * t213 + qJD(2);
t206 = qJDD(3) * t149;
t77 = -t152 * qJD(4) + t223 * t214;
t78 = -qJD(3) * t100 - t149 * qJD(4);
t43 = t234 * t77 + t235 * t78;
t87 = t149 * t190 - t152 * t191;
t88 = -t149 * t191 - t152 * t190;
t49 = -pkin(4) * t87 - pkin(8) * t88 + t209;
t203 = t148 * t49 + t151 * t43 + t58 * t211;
t200 = t94 * t212;
t199 = t94 * t211;
t198 = t150 * t136 + t220;
t135 = t153 * qJ(2);
t195 = -t150 * pkin(1) + t135;
t192 = -qJD(5) * t39 - t14;
t188 = t263 * t148;
t184 = -qJD(5) * t162 - qJD(1);
t183 = qJDD(2) - t233;
t120 = -t235 * pkin(3) - pkin(4);
t182 = -t30 * t211 + t20;
t167 = t151 * t14 + t148 * t21 + t39 * t211 - t30 * t212;
t2 = -qJ(6) * t26 - qJD(6) * t65 + t167;
t181 = -t263 * t6 + t2;
t37 = t235 * t76 - t73;
t29 = -qJD(3) * pkin(4) - t37;
t178 = t13 * t94 + t29 * t88;
t177 = -t162 * t25 - t67 * t87;
t176 = -t94 * t25 + t88 * t67;
t175 = -t162 * t26 + t65 * t87;
t174 = -t263 * t88 + t57 * t94;
t42 = t234 * t78 - t235 * t77;
t47 = t234 * t80 + t197;
t61 = t235 * t100 - t234 * t99;
t171 = -qJ(6) * t88 - qJD(6) * t94;
t124 = pkin(4) + t246;
t146 = -qJ(6) - pkin(8);
t170 = t124 * t128 + t129 * t146;
t168 = -t52 + (-t148 * t86 - t212) * t263;
t166 = t119 * t57 + t263 * t29;
t165 = 0.2e1 * qJ(2) * t208 + qJDD(3) * t154;
t5 = t26 * pkin(5) + qJDD(6) + t13;
t159 = -t162 * t18 - t17 * t94 - t37 * t88 + t38 * t87 + t261;
t157 = -t154 * t155 + t258;
t132 = qJDD(3) * t152;
t92 = t222 * t151;
t91 = t222 * t148;
t85 = t128 * t225 - t228;
t83 = t128 * t226 + t227;
t64 = t65 ^ 2;
t54 = t151 * t58;
t41 = t151 * t49;
t23 = t65 * pkin(5) + qJD(6) + t29;
t22 = -t94 * t231 + t240;
t12 = pkin(5) * t162 - t148 * t62 - t94 * t230 + t54;
t4 = -qJ(6) * t199 + (-qJD(5) * t62 + t171) * t148 + t203;
t3 = -t87 * pkin(5) - t148 * t43 + t41 + t171 * t151 + (-t60 + (qJ(6) * t94 - t58) * t148) * qJD(5);
t7 = [qJDD(1), t261, t180, qJDD(2) - t261 - 0.2e1 * t233, t258, -t183 * pkin(1) - g(1) * t195 - g(2) * t220 + (t202 + t141) * qJ(2), qJDD(1) * t145 - 0.2e1 * t149 * t196, -0.2e1 * t149 * t204 + 0.2e1 * t218 * t208, -t149 * t155 + t132, -t152 * t155 - t206, 0, t157 * t149 + t152 * t165, -t149 * t165 + t152 * t157, -t161 * t61 + t42 * t89 - t43 * t86 + t62 * t59 + t159, t18 * t62 + t38 * t43 - t17 * t61 - t37 * t42 + t169 * t224 + t98 * t209 - g(1) * (t195 + t262) - g(2) * (-t153 * t147 + t198) t176 * t151 - t67 * t200 (-t151 * t65 - t237) * t88 + (t236 - t151 * t26 + (t148 * t65 - t151 * t67) * qJD(5)) * t94, -t174 * t151 - t200 * t263 + t177, t174 * t148 - t199 * t263 + t175, -t162 * t57 - t263 * t87 (-t62 * t211 + t41) * t263 - t54 * t57 + t182 * t162 - t15 * t87 + t42 * t65 + t61 * t26 + t29 * t199 - g(1) * t85 - g(2) * t83 + ((-qJD(5) * t58 - t43) * t263 + t62 * t57 + t192 * t162 + t178) * t148 -(-t62 * t212 + t203) * t263 + t240 * t57 - t167 * t162 + t16 * t87 + t42 * t67 - t61 * t25 - t29 * t200 + g(1) * t84 - g(2) * t82 + t178 * t151, t12 * t25 - t22 * t26 - t3 * t67 - t4 * t65 + (-t148 * t9 - t151 * t6) * t88 + t180 * t129 + (-t1 * t151 - t2 * t148 + (t148 * t6 - t151 * t9) * qJD(5)) * t94, t2 * t22 + t9 * t4 + t1 * t12 + t6 * t3 + t5 * (t94 * t250 + t61) + t23 * ((t148 * t88 + t199) * pkin(5) + t42) - g(1) * (t135 + t170 * t153 + (-pkin(1) - t250) * t150 + t262) - g(2) * ((-t147 + t250) * t153 + t170 * t150 + t198); 0, 0, 0, qJDD(1), -t156, t183 + t164, 0, 0, 0, 0, 0, t217 * t149 + t132, t217 * t152 - t206, t161 * t94 + t162 * t59 + t87 * t86 - t88 * t89, -t159 - t229, 0, 0, 0, 0, 0, t162 * t238 - t94 * t26 - t88 * t65 + (t148 * t87 + t151 * t184) * t263, t162 * t52 + (-t148 * t184 + t151 * t87) * t263 - t176 (-t184 * t67 + t175) * t151 + (-t184 * t65 + t177) * t148, -t23 * t88 - t5 * t94 + (t162 * t2 + t184 * t6 - t87 * t9) * t151 + (-t1 * t162 + t184 * t9 + t6 * t87) * t148 - t261; 0, 0, 0, 0, 0, 0, t152 * t156 * t149, -t218 * t156, t204, -t205, qJDD(3), t152 * t164 + t247 + t96, g(3) * t152 + (-t104 - t164) * t149 (t38 - t47) * t89 - (-t48 + t37) * t86 + (t235 * t161 + t234 * t59) * pkin(3), t37 * t47 - t38 * t48 + (t234 * t18 + t235 * t17 + t247 + (-t261 - t229) * t152) * pkin(3), t67 * t189 - t236 (-t25 - t245) * t151 - t263 * t237 + t242, -t243 + t265, t168 + t244, -t263 * t89, t120 * t26 - t15 * t89 - t46 * t263 - t47 * t65 + (t263 * t48 + t166) * t148 - t259 * t151, -t120 * t25 + t259 * t148 + t166 * t151 + t16 * t89 + t241 * t263 - t47 * t67, -t128 * t261 - t257 * t148 + t181 * t151 - t239 * t65 - t91 * t25 - t251 * t67 - t92 * t26 - t248, t2 * t92 - t1 * t91 + t5 * (t120 - t246) - g(3) * (-t136 - t170) + t239 * t9 + t251 * t6 + (pkin(5) * t188 - t47) * t23 - t261 * (pkin(3) * t152 + t124 * t129 - t128 * t146); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86 ^ 2 - t89 ^ 2, t37 * t89 + t38 * t86 + t169 - t180, 0, 0, 0, 0, 0, t168 - t244, -t243 - t265 (t25 - t245) * t151 + t67 * t188 + t242, t181 * t148 + t257 * t151 - t23 * t89 - t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t65, -t64 + t256, t263 * t65 - t25, -t158 + (-qJD(5) + t263) * t67, -t57, t16 * t263 - t29 * t67 + (t192 + t248) * t148 + t182 + t260, g(1) * t83 - g(2) * t85 + t15 * t263 + t151 * t248 + t29 * t65 - t167, pkin(5) * t25 - t254 * t65, t254 * t9 + (t148 * t248 - t23 * t67 + t1 + t260) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64 - t256, t6 * t67 + t9 * t65 + t160 + t5;];
tau_reg  = t7;
