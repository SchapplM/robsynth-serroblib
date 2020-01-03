% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR13_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:33
% EndTime: 2019-12-31 21:46:44
% DurationCPUTime: 3.49s
% Computational Cost: add. (3309->382), mult. (9026->536), div. (0->0), fcn. (6719->8), ass. (0->182)
t145 = sin(qJ(3));
t148 = cos(qJ(3));
t237 = cos(pkin(5));
t194 = t237 * qJD(1);
t169 = t194 + qJD(2);
t146 = sin(qJ(2));
t143 = sin(pkin(5));
t225 = qJD(1) * t143;
t207 = t146 * t225;
t89 = t145 * t169 + t148 * t207;
t82 = qJD(5) + t89;
t87 = t145 * t207 - t148 * t169;
t263 = qJD(3) * t87;
t187 = pkin(1) * t194;
t133 = t146 * t187;
t149 = cos(qJ(2));
t224 = qJD(1) * t149;
t205 = t143 * t224;
t106 = pkin(7) * t205 + t133;
t69 = pkin(8) * t169 + t106;
t102 = (-pkin(2) * t149 - pkin(8) * t146 - pkin(1)) * t143;
t81 = qJD(1) * t102;
t36 = t145 * t69 - t148 * t81;
t228 = -qJD(4) - t36;
t127 = -qJD(3) + t205;
t216 = qJD(1) * qJD(2);
t202 = t143 * t216;
t130 = t146 * t202;
t175 = pkin(3) * t130;
t219 = qJD(3) * t148;
t221 = qJD(3) * t145;
t164 = t143 * (pkin(2) * t146 - pkin(8) * t149);
t105 = qJD(2) * t164;
t97 = qJD(1) * t105;
t200 = t149 * t237;
t232 = t143 * t146;
t261 = pkin(1) * t200 - pkin(7) * t232;
t107 = t261 * qJD(2);
t98 = qJD(1) * t107;
t198 = t145 * t98 - t148 * t97 + t69 * t219 + t81 * t221;
t14 = -t175 + t198;
t37 = t145 * t81 + t148 * t69;
t34 = t127 * qJ(4) - t37;
t262 = -t34 * t127 + t14;
t186 = t149 * t202;
t235 = qJD(3) * t89;
t55 = t145 * t186 + t235;
t189 = t145 * t205;
t260 = -t145 * qJD(4) - t106 + (-t189 + t221) * pkin(3);
t227 = t89 * pkin(4) - t228;
t255 = t87 * pkin(4);
t22 = -t34 - t255;
t257 = pkin(3) + pkin(9);
t54 = -t148 * t186 + t263;
t259 = t257 * t54 + (t22 - t37 + t255) * t82;
t191 = t237 * t146 * pkin(1);
t231 = t143 * t149;
t101 = pkin(7) * t231 + t237 * pkin(8) + t191;
t162 = t101 * t221 - t102 * t219 - t145 * t105 - t148 * t107;
t223 = qJD(2) * t146;
t15 = -t143 * (qJ(4) * t223 - qJD(4) * t149) + t162;
t144 = sin(qJ(5));
t147 = cos(qJ(5));
t174 = t147 * t127 - t144 * t87;
t21 = -qJD(5) * t174 + t130 * t144 - t147 * t55;
t16 = t257 * t127 + t227;
t103 = -pkin(7) * t207 + t149 * t187;
t68 = -pkin(2) * t169 - t103;
t152 = -t89 * qJ(4) + t68;
t19 = t257 * t87 + t152;
t178 = t144 * t19 - t147 * t16;
t206 = t143 * t223;
t173 = t257 * t206;
t7 = -t54 * pkin(4) - qJD(1) * t173 + t198;
t99 = pkin(7) * t186 + qJD(2) * t133;
t157 = t54 * qJ(4) - t89 * qJD(4) + t99;
t8 = t257 * t55 + t157;
t1 = -qJD(5) * t178 + t144 * t7 + t147 * t8;
t258 = t89 ^ 2;
t151 = qJD(1) ^ 2;
t256 = pkin(4) + pkin(8);
t28 = t87 * pkin(3) + t152;
t254 = t28 * t89;
t51 = -t144 * t127 - t147 * t87;
t253 = t51 * t82;
t252 = t174 * t82;
t251 = t89 * t87;
t250 = t148 * t101 + t145 * t102;
t104 = qJD(1) * t164;
t249 = t148 * t103 + t145 * t104;
t248 = t127 * t51;
t247 = t127 * t174;
t246 = t127 * t87;
t245 = t144 * t54;
t244 = t144 * t82;
t48 = t147 * t54;
t217 = qJD(5) * t147;
t218 = qJD(5) * t144;
t20 = t127 * t218 + t147 * t130 + t144 * t55 + t87 * t217;
t243 = t20 * t147;
t241 = t87 * qJ(4);
t240 = t89 * t127;
t236 = qJ(4) * t148;
t239 = -qJ(4) * t219 + t205 * t236 + t260;
t230 = t145 * t149;
t238 = -t256 * t221 - (-pkin(4) * t230 + qJ(4) * t146) * t225 - t249;
t234 = t127 * t148;
t140 = t143 ^ 2;
t233 = t140 * t151;
t229 = t148 * t149;
t204 = qJD(2) * t231;
t108 = pkin(7) * t204 + qJD(2) * t191;
t226 = t146 ^ 2 - t149 ^ 2;
t222 = qJD(2) * t148;
t220 = qJD(3) * t147;
t215 = t145 * pkin(8) * t127;
t214 = pkin(8) * t234;
t213 = pkin(8) * t223;
t212 = pkin(8) * t222;
t129 = t256 * t148;
t211 = t146 * t233;
t210 = t145 * t232;
t209 = t147 * t231;
t203 = t140 * t216;
t201 = -t145 * qJ(4) - pkin(2);
t199 = -t145 * t97 - t148 * t98 - t81 * t219 + t69 * t221;
t197 = -t145 * t101 + t148 * t102;
t196 = -t145 * t103 + t148 * t104;
t192 = 0.2e1 * t203;
t185 = -qJD(4) * t127 - t199;
t128 = t256 * t145;
t184 = -qJD(5) * t128 - t260 + t127 * (pkin(9) * t145 - t236);
t40 = pkin(3) * t231 - t197;
t183 = -0.2e1 * pkin(1) * t203;
t112 = -t257 * t148 + t201;
t181 = qJD(5) * t112 - qJD(3) * t129 + (pkin(4) * t229 - t257 * t146) * t225 - t196;
t4 = t144 * t16 + t147 * t19;
t111 = t237 * t145 + t148 * t232;
t25 = t111 * pkin(4) + pkin(9) * t231 + t40;
t110 = -t237 * t148 + t210;
t100 = -t237 * pkin(2) - t261;
t153 = -t111 * qJ(4) + t100;
t29 = t257 * t110 + t153;
t177 = -t144 * t29 + t147 * t25;
t176 = t144 * t25 + t147 * t29;
t172 = -t101 * t219 - t102 * t221 + t148 * t105 - t145 * t107;
t171 = qJ(4) * t130;
t168 = 0.2e1 * t194 + qJD(2);
t39 = qJ(4) * t231 - t250;
t166 = -t244 * t82 - t48;
t165 = -t37 * t127 - t198;
t62 = t110 * t147 + t144 * t231;
t11 = -t171 - t185;
t6 = -t55 * pkin(4) - t11;
t163 = t6 + (t82 * t257 + t241) * t82;
t61 = -qJD(3) * t210 + (t237 * qJD(3) + t204) * t148;
t156 = -t61 * qJ(4) - t111 * qJD(4) + t108;
t155 = -t82 ^ 2 * t147 + t245;
t2 = -qJD(5) * t4 - t144 * t8 + t147 * t7;
t154 = -t54 - t246;
t122 = -t148 * pkin(3) + t201;
t76 = (t144 * t230 + t146 * t147) * t225;
t75 = t144 * t207 - t147 * t189;
t63 = t110 * t144 - t209;
t60 = qJD(3) * t111 + t145 * t204;
t50 = t54 * t145;
t45 = t89 * pkin(3) + t241;
t44 = -pkin(3) * t207 - t196;
t42 = -qJ(4) * t207 - t249;
t41 = t54 * t111;
t38 = t110 * pkin(3) + t153;
t32 = t127 * pkin(3) - t228;
t30 = -t110 * pkin(4) - t39;
t27 = qJD(5) * t62 + t60 * t144 + t147 * t206;
t26 = -qJD(5) * t209 - t60 * t147 + (qJD(5) * t110 + t206) * t144;
t18 = t60 * pkin(3) + t156;
t17 = -pkin(3) * t206 - t172;
t13 = t55 * pkin(3) + t157;
t12 = t257 * t60 + t156;
t10 = -t60 * pkin(4) - t15;
t9 = t61 * pkin(4) - t172 - t173;
t3 = [0, 0, 0, t146 * t149 * t192, -t226 * t192, t168 * t204, -t168 * t206, 0, -t108 * t169 + t146 * t183 - t99 * t237, -t107 * t169 + t149 * t183 - t98 * t237, t89 * t61 - t41, t54 * t110 - t111 * t55 - t89 * t60 - t61 * t87, -t61 * t127 + (t149 * t54 + (qJD(1) * t111 + t89) * t223) * t143, t60 * t127 + (t149 * t55 + (-qJD(1) * t110 - t87) * t223) * t143, (-t127 * t143 - t140 * t224) * t223, -t172 * t127 + t108 * t87 + t100 * t55 + t99 * t110 + t68 * t60 + (t198 * t149 + (qJD(1) * t197 - t36) * t223) * t143, -t162 * t127 + t108 * t89 - t100 * t54 + t99 * t111 + t68 * t61 + (-t199 * t149 + (-t250 * qJD(1) - t37) * t223) * t143, t11 * t110 + t14 * t111 + t15 * t87 + t17 * t89 + t32 * t61 + t34 * t60 + t39 * t55 - t40 * t54, -t13 * t110 - t17 * t127 - t18 * t87 - t28 * t60 - t38 * t55 + (-t14 * t149 + (qJD(1) * t40 + t32) * t223) * t143, -t13 * t111 + t15 * t127 - t18 * t89 - t28 * t61 + t38 * t54 + (t11 * t149 + (-qJD(1) * t39 - t34) * t223) * t143, t11 * t39 + t13 * t38 + t14 * t40 + t34 * t15 + t32 * t17 + t28 * t18, -t174 * t27 + t20 * t63, t174 * t26 + t20 * t62 - t63 * t21 - t27 * t51, t20 * t111 - t174 * t61 + t27 * t82 - t63 * t54, -t21 * t111 - t26 * t82 - t51 * t61 - t62 * t54, t82 * t61 - t41, (-qJD(5) * t176 - t144 * t12 + t147 * t9) * t82 - t177 * t54 + t2 * t111 - t178 * t61 + t10 * t51 + t30 * t21 - t6 * t62 + t22 * t26, -(qJD(5) * t177 + t147 * t12 + t144 * t9) * t82 + t176 * t54 - t1 * t111 - t4 * t61 - t10 * t174 + t30 * t20 + t6 * t63 + t22 * t27; 0, 0, 0, -t149 * t211, t226 * t233, -t143 * t151 * t200, t169 * t207 - t130, 0, pkin(1) * t211 + t106 * t169 - t99, pkin(7) * t130 + t103 * t169 + (-qJD(2) * t194 + t233) * t149 * pkin(1), -t234 * t89 - t50, (-t54 + t246) * t148 + (-t55 + t240) * t145, -t127 * t219 + (t127 * t229 + (t145 * qJD(2) - t89) * t146) * t225, t127 * t221 + (-t127 * t230 + (t87 + t222) * t146) * t225, t127 * t207, -pkin(2) * t55 - t99 * t148 + t196 * t127 - t106 * t87 + (t68 * t145 + t214) * qJD(3) + (t36 * t146 + (-t149 * t68 - t213) * t145) * t225, pkin(2) * t54 + t99 * t145 - t249 * t127 - t106 * t89 + (t68 * t148 - t215) * qJD(3) + (-t68 * t229 + (t37 - t212) * t146) * t225, -t42 * t87 - t44 * t89 + (-t11 - t127 * t32 + (-t55 + t235) * pkin(8)) * t148 + ((-t54 + t263) * pkin(8) + t262) * t145, -t122 * t55 + t44 * t127 + t13 * t148 - t239 * t87 + (-t145 * t28 - t214) * qJD(3) + (-t146 * t32 + (t149 * t28 + t213) * t145) * t225, t122 * t54 - t42 * t127 - t13 * t145 - t239 * t89 + (-t148 * t28 + t215) * qJD(3) + (t28 * t229 + (t34 + t212) * t146) * t225, t13 * t122 - t32 * t44 - t34 * t42 + t239 * t28 + (-t11 * t148 + t14 * t145 + (t145 * t34 + t148 * t32) * qJD(3)) * pkin(8), -t20 * t144 * t148 - (t144 * t221 - t148 * t217 - t76) * t174, t76 * t51 - t174 * t75 + (-t144 * t51 - t147 * t174) * t221 + (t144 * t21 - t243 + (-t144 * t174 + t147 * t51) * qJD(5)) * t148, -t76 * t82 + (qJD(3) * t244 + t20) * t145 + (-t217 * t82 + t245 + t247) * t148, t75 * t82 + (t220 * t82 - t21) * t145 + (t218 * t82 + t248 + t48) * t148, -t234 * t82 - t50, -(-t144 * t112 + t147 * t128) * t54 + t129 * t21 - t22 * t75 + (t144 * t184 - t147 * t181) * t82 + t238 * t51 + (-t22 * t220 + t2) * t145 + (t127 * t178 + t6 * t147 - t218 * t22) * t148, -t1 * t145 + t129 * t20 - t22 * t76 - t238 * t174 + (t112 * t54 + t184 * t82) * t147 + (t128 * t54 + t181 * t82 + t22 * t221) * t144 + (t127 * t4 - t6 * t144 - t217 * t22) * t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, -t87 ^ 2 + t258, t154, -t240 - t55, t130, -t68 * t89 + t165, t36 * t127 + t68 * t87 + t199, pkin(3) * t54 - qJ(4) * t55 + (-t34 - t37) * t89 + (t32 + t228) * t87, t45 * t87 - t165 - 0.2e1 * t175 + t254, t228 * t127 - t28 * t87 + t45 * t89 + 0.2e1 * t171 + t185, -t14 * pkin(3) - t11 * qJ(4) + t228 * t34 - t28 * t45 - t32 * t37, t174 * t244 + t243, (-t21 + t252) * t147 + (-t20 + t253) * t144, -t174 * t87 + t166, -t51 * t87 + t155, t82 * t87, qJ(4) * t21 + t163 * t144 + t259 * t147 - t178 * t87 + t227 * t51, qJ(4) * t20 - t259 * t144 + t163 * t147 - t174 * t227 - t4 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, t130 - t251, -t127 ^ 2 - t258, t254 + t262, 0, 0, 0, 0, 0, t166 + t248, t155 - t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174 * t51, t174 ^ 2 - t51 ^ 2, t20 + t253, -t21 - t252, -t54, t174 * t22 + t4 * t82 + t2, -t178 * t82 + t22 * t51 - t1;];
tauc_reg = t3;
