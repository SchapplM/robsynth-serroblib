% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRP11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP11_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:38
% EndTime: 2019-12-31 18:54:45
% DurationCPUTime: 3.85s
% Computational Cost: add. (4747->422), mult. (11340->501), div. (0->0), fcn. (8387->10), ass. (0->203)
t143 = sin(pkin(8));
t147 = sin(qJ(3));
t144 = cos(pkin(8));
t259 = cos(qJ(3));
t204 = t259 * t144;
t171 = -t147 * t143 + t204;
t283 = t171 * qJD(1);
t284 = t283 * qJD(3);
t90 = qJD(4) - t283;
t234 = qJDD(1) * pkin(1);
t133 = qJDD(2) - t234;
t148 = sin(qJ(1));
t150 = cos(qJ(1));
t274 = g(1) * t148 - g(2) * t150;
t174 = -t133 + t274;
t106 = t259 * t143 + t147 * t144;
t165 = t106 * qJDD(1);
t154 = t165 + t284;
t282 = -qJD(3) * qJD(4) - t154;
t203 = qJD(3) * t259;
t216 = qJD(3) * t147;
t102 = t143 * t216 - t144 * t203;
t146 = sin(qJ(4));
t149 = cos(qJ(4));
t100 = t106 * qJD(1);
t214 = qJD(4) * t149;
t189 = -t149 * qJDD(3) + t100 * t214;
t217 = qJD(3) * t146;
t43 = t146 * t165 + (qJD(4) + t283) * t217 + t189;
t236 = t43 * t149;
t215 = qJD(4) * t146;
t42 = -t146 * qJDD(3) + t100 * t215 + t282 * t149;
t237 = t42 * t146;
t76 = -t149 * qJD(3) + t100 * t146;
t78 = t100 * t149 + t217;
t281 = ((t146 * t76 - t149 * t78) * qJD(4) - t236 + t237) * t106 + (t146 * t78 + t149 * t76) * t102;
t247 = pkin(6) + qJ(2);
t115 = t247 * t144;
t108 = qJD(1) * t115;
t114 = t247 * t143;
t107 = qJD(1) * t114;
t213 = qJD(1) * qJD(2);
t266 = t247 * qJDD(1) + t213;
t85 = t266 * t143;
t86 = t266 * t144;
t209 = t107 * t203 + t147 * t85 - t259 * t86;
t33 = -t108 * t216 - t209;
t29 = qJDD(3) * pkin(7) + t33;
t210 = t144 * qJDD(1);
t103 = t106 * qJD(3);
t211 = t143 * qJDD(1);
t186 = -qJDD(1) * t204 + t147 * t211;
t68 = qJD(1) * t103 + t186;
t32 = -pkin(2) * t210 + t68 * pkin(3) - pkin(7) * t154 + t133;
t129 = t144 * pkin(2) + pkin(1);
t112 = -t129 * qJD(1) + qJD(2);
t49 = -pkin(3) * t283 - pkin(7) * t100 + t112;
t70 = -t147 * t107 + t259 * t108;
t64 = qJD(3) * pkin(7) + t70;
t200 = t146 * t29 - t149 * t32 + t64 * t214 + t49 * t215;
t62 = qJDD(4) + t68;
t262 = pkin(4) * t62;
t2 = qJDD(5) + t200 - t262;
t21 = t146 * t49 + t149 * t64;
t17 = qJ(5) * t90 + t21;
t253 = t17 * t90;
t280 = -t2 + t253;
t251 = t21 * t90;
t279 = -t200 + t251;
t54 = t146 * t62;
t173 = t90 * t214 + t54;
t55 = t149 * t62;
t278 = -t90 * t215 + t55;
t142 = pkin(8) + qJ(3);
t134 = sin(t142);
t192 = g(1) * t150 + g(2) * t148;
t276 = t192 * t134;
t275 = -t259 * t114 - t147 * t115;
t135 = cos(t142);
t219 = t135 * pkin(3) + t134 * pkin(7);
t272 = t100 * qJD(3);
t271 = qJ(2) * qJDD(1);
t97 = t147 * t108;
t69 = -t259 * t107 - t97;
t63 = -qJD(3) * pkin(3) - t69;
t22 = t76 * pkin(4) - t78 * qJ(5) + t63;
t261 = pkin(7) * t62;
t270 = t90 * t22 - t261;
t233 = t102 * t146;
t269 = t103 * t76 + t106 * t173 - t171 * t43 - t90 * t233;
t255 = g(3) * t134;
t268 = t192 * t135 + t255;
t67 = -pkin(3) * t171 - pkin(7) * t106 - t129;
t75 = -t147 * t114 + t259 * t115;
t245 = t146 * t67 + t149 * t75;
t50 = t171 * qJD(2) + t275 * qJD(3);
t66 = pkin(3) * t103 + pkin(7) * t102;
t11 = -qJD(4) * t245 - t146 * t50 + t149 * t66;
t265 = t78 ^ 2;
t264 = t90 ^ 2;
t263 = t100 ^ 2;
t260 = pkin(7) * t78;
t254 = g(3) * t135;
t20 = -t146 * t64 + t149 * t49;
t252 = t20 * t90;
t250 = t76 * t90;
t249 = t76 * t283;
t248 = t78 * t76;
t246 = -t146 * t43 - t76 * t214;
t65 = pkin(3) * t100 - pkin(7) * t283;
t36 = t146 * t65 + t149 * t69;
t187 = pkin(4) * t146 - qJ(5) * t149;
t244 = -t146 * qJD(5) + t90 * t187 - t70;
t243 = qJ(5) * t62;
t242 = t100 * t76;
t241 = t100 * t78;
t240 = t100 * t283;
t239 = t146 * t283;
t238 = t149 * t90;
t235 = t90 * t100;
t232 = t102 * t149;
t231 = t134 * t148;
t230 = t134 * t150;
t229 = t135 * t150;
t228 = t247 * t150;
t227 = t146 * t148;
t224 = t148 * t149;
t223 = t149 * t150;
t222 = t150 * t146;
t221 = qJD(5) - t20;
t220 = (g(1) * t223 + g(2) * t224) * t134;
t140 = t143 ^ 2;
t141 = t144 ^ 2;
t218 = t140 + t141;
t208 = pkin(7) * qJD(4) * t90;
t207 = t76 ^ 2 - t265;
t72 = t78 * t215;
t116 = t150 * t129;
t205 = pkin(3) * t229 + pkin(7) * t230 + t116;
t34 = t107 * t216 - t108 * t203 - t147 * t86 - t259 * t85;
t199 = t146 * t90;
t198 = t218 * qJD(1) ^ 2;
t197 = 0.2e1 * t218;
t196 = t78 * t239 - t72;
t91 = t135 * t227 + t223;
t93 = t135 * t222 - t224;
t195 = g(1) * t91 - g(2) * t93;
t92 = t135 * t224 - t222;
t94 = t135 * t223 + t227;
t194 = g(1) * t92 - g(2) * t94;
t193 = -g(1) * t231 + g(2) * t230;
t190 = (qJD(4) * t76 - t42) * pkin(7);
t188 = pkin(4) * t149 + qJ(5) * t146;
t16 = -pkin(4) * t90 + t221;
t185 = t146 * t17 - t149 * t16;
t184 = t146 * t21 + t149 * t20;
t35 = -t146 * t69 + t149 * t65;
t40 = -t146 * t75 + t149 * t67;
t178 = -t238 * t283 + t173;
t177 = t90 * t239 + t278;
t176 = t208 + t254;
t175 = pkin(3) + t188;
t30 = -qJDD(3) * pkin(3) - t34;
t172 = t90 * t63 - t261;
t3 = t146 * t32 + t149 * t29 + t49 * t214 - t64 * t215;
t10 = t146 * t66 + t149 * t50 + t67 * t214 - t75 * t215;
t170 = t176 + t30;
t111 = -t129 * qJDD(1) + qJDD(2);
t169 = (t42 + t249) * t149 + t246;
t168 = t76 * t199 - t236;
t167 = -t177 - t242;
t164 = t174 + t234;
t162 = g(1) * t93 + g(2) * t91 + t146 * t255 - t200;
t161 = -t254 + t276;
t160 = -pkin(7) * t236 - t268;
t159 = (-g(1) * (-t129 - t219) - g(2) * t247) * t148;
t158 = -t246 * t106 - t76 * t233;
t157 = t197 * t213 - t192;
t156 = t22 * t78 + qJDD(5) - t162;
t155 = -g(1) * t94 - g(2) * t92 - t149 * t255 + t3;
t51 = qJD(2) * t106 + qJD(3) * t75;
t153 = t282 * t146 + t78 * t90 - t189;
t119 = pkin(7) * t229;
t117 = t148 * t135 * pkin(7);
t96 = t283 ^ 2;
t45 = pkin(4) * t78 + qJ(5) * t76;
t44 = t187 * t106 - t275;
t31 = t103 * t90 - t171 * t62;
t28 = pkin(4) * t171 - t40;
t27 = -qJ(5) * t171 + t245;
t19 = -pkin(4) * t100 - t35;
t18 = qJ(5) * t100 + t36;
t15 = -t42 + t250;
t14 = t178 - t241;
t13 = -t187 * t102 + (t188 * qJD(4) - qJD(5) * t149) * t106 + t51;
t12 = t78 * t238 - t237;
t9 = -t78 * t232 + (-t149 * t42 - t72) * t106;
t8 = -pkin(4) * t103 - t11;
t7 = qJ(5) * t103 - qJD(5) * t171 + t10;
t6 = t103 * t78 + t278 * t106 + t171 * t42 - t90 * t232;
t5 = pkin(4) * t43 + qJ(5) * t42 - qJD(5) * t78 + t30;
t1 = qJD(5) * t90 + t243 + t3;
t4 = [0, 0, 0, 0, 0, qJDD(1), t274, t192, 0, 0, t140 * qJDD(1), 0.2e1 * t143 * t210, 0, t141 * qJDD(1), 0, 0, t164 * t144, -t164 * t143, t197 * t271 + t157, t174 * pkin(1) + (t218 * t271 + t157) * qJ(2), -t100 * t102 + t106 * t154, -t100 * t103 - t102 * t283 - t106 * t68 + t154 * t171, -qJD(3) * t102 + qJDD(3) * t106, -t103 * t283 - t171 * t68, -qJD(3) * t103 + qJDD(3) * t171, 0, -qJD(3) * t51 + qJDD(3) * t275 + t103 * t112 - t111 * t171 - t129 * t68 + t135 * t274, -t50 * qJD(3) - t75 * qJDD(3) - t112 * t102 + t111 * t106 - t129 * t154 + t193, t51 * t100 + t69 * t102 - t70 * t103 - t34 * t106 - t154 * t275 + t33 * t171 + t283 * t50 - t75 * t68 - t192, t33 * t75 + t70 * t50 + t34 * t275 - t69 * t51 - t111 * t129 - g(1) * (-t129 * t148 + t228) - g(2) * (t148 * t247 + t116), t9, t281, t6, t158, -t269, t31, -t63 * t233 + t103 * t20 + t171 * t200 + t11 * t90 + t40 * t62 - t43 * t275 + t51 * t76 + (t146 * t30 + t63 * t214) * t106 + t194, -t63 * t232 - t10 * t90 - t103 * t21 + t171 * t3 - t245 * t62 + t42 * t275 + t51 * t78 + (t149 * t30 - t215 * t63) * t106 - t195, -t10 * t76 - t11 * t78 + t40 * t42 - t245 * t43 + t184 * t102 + (-t146 * t3 + t149 * t200 + (t146 * t20 - t149 * t21) * qJD(4)) * t106 - t193, -g(1) * t228 - g(2) * t205 + t21 * t10 + t20 * t11 - t200 * t40 + t245 * t3 - t275 * t30 + t63 * t51 + t159, t9, t6, -t281, t31, t269, t158, -t22 * t233 - t103 * t16 + t171 * t2 + t13 * t76 - t28 * t62 + t43 * t44 - t8 * t90 + (t146 * t5 + t214 * t22) * t106 + t194, -t27 * t43 - t28 * t42 - t7 * t76 + t78 * t8 + t185 * t102 + (-t1 * t146 + t149 * t2 + (-t146 * t16 - t149 * t17) * qJD(4)) * t106 - t193, t22 * t232 - t1 * t171 + t103 * t17 - t13 * t78 + t27 * t62 + t42 * t44 + t7 * t90 + (-t149 * t5 + t215 * t22) * t106 + t195, t1 * t27 + t17 * t7 + t5 * t44 + t22 * t13 + t2 * t28 + t16 * t8 - g(1) * (-pkin(4) * t92 - qJ(5) * t91 + t228) - g(2) * (pkin(4) * t94 + qJ(5) * t93 + t205) + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210, t211, -t198, -qJ(2) * t198 - t174, 0, 0, 0, 0, 0, 0, t186 + 0.2e1 * t272, t165 + 0.2e1 * t284, -t96 - t263, t100 * t69 - t283 * t70 + t111 - t274, 0, 0, 0, 0, 0, 0, t177 - t242, -t149 * t264 - t241 - t54, t199 * t78 + t169, -t100 * t63 + t279 * t149 + (t3 - t252) * t146 - t274, 0, 0, 0, 0, 0, 0, -t146 * t264 - t242 + t55, t169 - t196, t178 + t241, -t100 * t22 + t280 * t149 + (t90 * t16 + t1) * t146 - t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t240, -t96 + t263, t165, t240, -t186, qJDD(3), qJD(3) * t70 - t100 * t112 + t161 + t34, -t112 * t283 + (t69 + t97) * qJD(3) + t209 + t268, 0, 0, t12, (-t42 + t249) * t149 + t196 + t246, t14, t168, -t167, -t235, -pkin(3) * t43 - t100 * t20 + t146 * t172 - t149 * t170 - t35 * t90 - t70 * t76 + t220, pkin(3) * t42 + t100 * t21 + t36 * t90 - t70 * t78 + t172 * t149 + (t170 - t276) * t146, t35 * t78 + t36 * t76 + (t20 * t283 + t3 + (-t20 + t260) * qJD(4)) * t149 + (t190 - t279) * t146 + t160, -t30 * pkin(3) - t21 * t36 - t20 * t35 - t63 * t70 - g(1) * (-pkin(3) * t230 + t119) - g(2) * (-pkin(3) * t231 + t117) - g(3) * t219 + (-t184 * qJD(4) + t146 * t200 + t3 * t149) * pkin(7), t12, t14, t72 + (-t283 * t78 + t43) * t146 + (t42 + t250) * t149, -t235, t167, t168, t100 * t16 - t175 * t43 + t19 * t90 + t244 * t76 + (-t176 - t5) * t149 + t270 * t146 + t220, t18 * t76 - t19 * t78 + (-t16 * t283 + t1 + (t16 + t260) * qJD(4)) * t149 + (t190 - t280) * t146 + t160, -t100 * t17 - t175 * t42 - t18 * t90 - t244 * t78 - t270 * t149 + (t161 - t5 - t208) * t146, -t17 * t18 - t16 * t19 - g(1) * t119 - g(2) * t117 - g(3) * (t135 * t188 + t219) + t244 * t22 + (-qJD(4) * t185 + t1 * t149 + t2 * t146) * pkin(7) + (-t5 + t276) * t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, -t207, t15, -t248, t153, t62, -t63 * t78 + t162 + t251, t63 * t76 - t155 + t252, 0, 0, t248, t15, t207, t62, -t153, -t248, -t45 * t76 - t156 + t251 + 0.2e1 * t262, pkin(4) * t42 - qJ(5) * t43 + (t17 - t21) * t78 + (t16 - t221) * t76, 0.2e1 * t243 - t22 * t76 + t45 * t78 + (0.2e1 * qJD(5) - t20) * t90 + t155, t1 * qJ(5) - t2 * pkin(4) - t22 * t45 - t16 * t21 - g(1) * (-pkin(4) * t93 + qJ(5) * t94) - g(2) * (-pkin(4) * t91 + qJ(5) * t92) + t221 * t17 + t187 * t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t186 + t248 - t272, t15, -t264 - t265, t156 - t253 - t262;];
tau_reg = t4;
