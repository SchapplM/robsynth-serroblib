% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR11_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:28:00
% EndTime: 2019-12-31 18:28:04
% DurationCPUTime: 2.22s
% Computational Cost: add. (3118->327), mult. (7406->392), div. (0->0), fcn. (5490->10), ass. (0->184)
t153 = sin(pkin(8));
t243 = pkin(6) + qJ(2);
t117 = t243 * t153;
t110 = qJD(1) * t117;
t249 = cos(qJ(3));
t154 = cos(pkin(8));
t118 = t243 * t154;
t111 = qJD(1) * t118;
t157 = sin(qJ(3));
t96 = t157 * t111;
t65 = -t249 * t110 - t96;
t220 = qJD(4) - t65;
t156 = sin(qJ(5));
t159 = cos(qJ(5));
t109 = t249 * t153 + t157 * t154;
t263 = t109 * qJD(1);
t206 = t249 * t154;
t192 = qJD(1) * t206;
t223 = t157 * t153;
t205 = qJD(1) * t223;
t98 = -t192 + t205;
t48 = -t156 * t98 - t159 * t263;
t200 = qJDD(1) * t249;
t210 = t154 * qJDD(1);
t208 = qJD(3) * t192 + t153 * t200 + t157 * t210;
t60 = qJD(3) * t205 - t208;
t103 = t109 * qJD(3);
t211 = t153 * qJDD(1);
t188 = -t154 * t200 + t157 * t211;
t61 = qJD(1) * t103 + t188;
t167 = t48 * qJD(5) + t156 * t60 + t159 * t61;
t149 = qJD(3) - qJD(5);
t276 = t48 * t149;
t278 = t167 + t276;
t245 = t48 ^ 2;
t185 = t156 * t263 - t159 * t98;
t246 = t185 ^ 2;
t277 = t245 - t246;
t244 = t185 * t48;
t138 = t154 * pkin(2) + pkin(1);
t112 = -t138 * qJDD(1) + qJDD(2);
t232 = t60 * qJ(4);
t275 = -t112 - t232;
t113 = -t138 * qJD(1) + qJD(2);
t274 = -t263 * qJ(4) + t113;
t230 = qJDD(1) * pkin(1);
t158 = sin(qJ(1));
t160 = cos(qJ(1));
t265 = g(1) * t158 - g(2) * t160;
t179 = -qJDD(2) + t230 + t265;
t146 = t153 ^ 2;
t147 = t154 ^ 2;
t217 = t146 + t147;
t194 = 0.2e1 * t217;
t212 = qJD(1) * qJD(2);
t264 = g(1) * t160 + g(2) * t158;
t273 = t194 * t212 - t264;
t214 = qJD(5) * t159;
t215 = qJD(5) * t156;
t177 = t156 * t61 - t159 * t60 + t98 * t214 - t215 * t263;
t266 = t185 * t149;
t272 = t177 - t266;
t271 = -pkin(7) * t263 + t220;
t148 = pkin(8) + qJ(3);
t140 = sin(t148);
t141 = cos(t148);
t184 = t140 * t156 + t141 * t159;
t254 = pkin(3) + pkin(4);
t28 = -t254 * t98 - t274;
t91 = t140 * t159 - t141 * t156;
t73 = t91 * t158;
t224 = t141 * t160;
t226 = t140 * t160;
t75 = t156 * t224 - t159 * t226;
t270 = g(1) * t75 - g(2) * t73 + g(3) * t184 + t28 * t48;
t34 = -t254 * qJD(3) + t271;
t152 = qJD(3) * qJ(4);
t66 = -t157 * t110 + t249 * t111;
t42 = t98 * pkin(7) + t66;
t35 = t152 + t42;
t12 = t156 * t34 + t159 * t35;
t150 = qJDD(3) * qJ(4);
t151 = qJD(3) * qJD(4);
t203 = qJD(3) * t249;
t257 = t243 * qJDD(1) + t212;
t78 = t257 * t153;
t79 = t257 * t154;
t209 = -t110 * t203 - t157 * t78 + t249 * t79;
t216 = qJD(3) * t157;
t22 = -t111 * t216 + t209;
t20 = t150 + t151 + t22;
t10 = t61 * pkin(7) + t20;
t199 = -t110 * t216 + t111 * t203 + t157 * t79 + t249 * t78;
t189 = qJDD(4) + t199;
t9 = t60 * pkin(7) - t254 * qJDD(3) + t189;
t2 = -t12 * qJD(5) - t156 * t10 + t159 * t9;
t269 = -t12 * t149 + t2;
t255 = t263 ^ 2;
t92 = t98 ^ 2;
t268 = -t92 - t255;
t267 = -t92 + t255;
t68 = -t157 * t117 + t249 * t118;
t219 = t141 * pkin(3) + t140 * qJ(4);
t262 = qJ(2) * qJDD(1);
t1 = t159 * t10 + t156 * t9 + t34 * t214 - t35 * t215;
t74 = t184 * t158;
t76 = t184 * t160;
t260 = g(1) * t76 + g(2) * t74 + g(3) * t91 + t185 * t28 - t1;
t213 = t263 * qJD(4);
t259 = -t254 * t61 + t213;
t44 = -t117 * t203 + qJD(2) * t206 + (-qJD(2) * t153 - qJD(3) * t118) * t157;
t258 = -t44 * qJD(3) - t68 * qJDD(3) - t140 * t265;
t256 = g(3) * t140 + (t65 + t96) * qJD(3) + t264 * t141 - t209;
t250 = t61 * pkin(3);
t247 = t141 * pkin(4);
t242 = pkin(7) - t243;
t114 = -t156 * qJ(4) - t159 * t254;
t241 = t114 * qJD(5) - t156 * t42 + t271 * t159;
t115 = t159 * qJ(4) - t156 * t254;
t240 = -t115 * qJD(5) - t271 * t156 - t159 * t42;
t239 = t263 * t98;
t11 = -t156 * t35 + t159 * t34;
t238 = t11 * t149;
t231 = t98 * qJ(4);
t229 = qJDD(3) * pkin(3);
t227 = t140 * t158;
t225 = t141 * t158;
t222 = t160 * t243;
t221 = t66 * qJD(3);
t207 = -g(1) * t226 - g(2) * t227 + g(3) * t141;
t197 = t217 * qJD(1) ^ 2;
t67 = t249 * t117 + t157 * t118;
t196 = t149 ^ 2;
t122 = t160 * t138;
t193 = g(2) * (pkin(3) * t224 + qJ(4) * t226 + t122);
t108 = -t206 + t223;
t187 = t98 * t103 + t61 * t108;
t46 = -t109 * pkin(7) + t67;
t47 = t108 * pkin(7) + t68;
t18 = -t156 * t47 + t159 * t46;
t19 = t156 * t46 + t159 * t47;
t64 = t156 * t108 + t159 * t109;
t102 = t153 * t216 - t154 * t203;
t183 = -t102 * qJ(4) + t109 * qJD(4);
t181 = qJD(3) * t103 + qJDD(3) * t108;
t180 = t109 * qJ(4) + t138;
t178 = -t138 - t219;
t175 = -t199 - t207;
t173 = t179 + t230;
t45 = t109 * qJD(2) + t68 * qJD(3);
t172 = t263 * t45 - t44 * t98 - t67 * t60 - t68 * t61 - t264;
t171 = g(1) * t225 - g(2) * t224 - t45 * qJD(3) - t67 * qJDD(3);
t170 = t102 * t98 - t103 * t263 + t60 * t108 - t109 * t61;
t169 = t112 - t265;
t43 = t98 * pkin(3) + t274;
t168 = t263 * t43 + qJDD(4) - t175;
t165 = t169 + t232;
t164 = 0.2e1 * t263 * qJD(3) + t188;
t145 = qJDD(3) - qJDD(5);
t121 = qJ(4) * t224;
t119 = qJ(4) * t225;
t63 = -t159 * t108 + t156 * t109;
t62 = -t102 * qJD(3) + t109 * qJDD(3);
t59 = t108 * pkin(3) - t180;
t58 = pkin(3) * t263 + t231;
t57 = t152 + t66;
t56 = -qJD(3) * pkin(3) + t220;
t40 = -t254 * t108 + t180;
t39 = t103 * pkin(3) - t183;
t38 = (t98 - t205) * qJD(3) + t208;
t37 = (t98 + t205) * qJD(3) - t208;
t36 = -t254 * t263 - t231;
t30 = t102 * pkin(7) + t45;
t29 = t103 * pkin(7) + t44;
t27 = -t254 * t103 + t183;
t26 = t64 * qJD(5) - t156 * t102 - t159 * t103;
t25 = t159 * t102 - t156 * t103 - t108 * t214 + t109 * t215;
t24 = -t102 * t263 - t60 * t109;
t21 = t189 - t229;
t17 = -t213 + t250 - t275;
t5 = t259 + t275;
t4 = -t19 * qJD(5) - t156 * t29 + t159 * t30;
t3 = t18 * qJD(5) + t156 * t30 + t159 * t29;
t6 = [0, 0, 0, 0, 0, qJDD(1), t265, t264, 0, 0, t146 * qJDD(1), 0.2e1 * t153 * t210, 0, t147 * qJDD(1), 0, 0, t173 * t154, -t173 * t153, t194 * t262 + t273, t179 * pkin(1) + (t217 * t262 + t273) * qJ(2), t24, t170, t62, t187, -t181, 0, t113 * t103 + t112 * t108 - t138 * t61 + t171, -t113 * t102 + t112 * t109 + t138 * t60 + t258, t65 * t102 - t66 * t103 - t22 * t108 + t109 * t199 + t172, t22 * t68 + t66 * t44 + t199 * t67 - t65 * t45 - t112 * t138 - g(1) * (-t158 * t138 + t222) - g(2) * (t158 * t243 + t122), t24, t62, -t170, 0, t181, t187, t43 * t103 + t17 * t108 + t39 * t98 + t59 * t61 + t171, -t56 * t102 - t57 * t103 - t20 * t108 + t21 * t109 + t172, t43 * t102 - t17 * t109 - t263 * t39 + t59 * t60 - t258, t20 * t68 + t57 * t44 + t17 * t59 + t43 * t39 + t21 * t67 + t56 * t45 - g(1) * t222 - t193 + (-g(1) * t178 - g(2) * t243) * t158, t177 * t64 + t25 * t48, t167 * t64 - t177 * t63 + t185 * t25 + t26 * t48, -t64 * t145 + t25 * t149, -t167 * t63 + t185 * t26, t63 * t145 + t26 * t149, 0, g(1) * t74 - g(2) * t76 - t18 * t145 - t4 * t149 - t167 * t40 + t185 * t27 + t28 * t26 + t5 * t63, g(1) * t73 + g(2) * t75 + t19 * t145 + t3 * t149 + t177 * t40 - t28 * t25 - t27 * t48 + t5 * t64, -t1 * t63 + t11 * t25 - t12 * t26 + t167 * t19 - t177 * t18 - t185 * t3 - t2 * t64 + t4 * t48 + t264, t1 * t19 + t12 * t3 + t2 * t18 + t11 * t4 + t5 * t40 + t28 * t27 - t193 + (g(1) * t242 - g(2) * t247) * t160 + (-g(1) * (t178 - t247) + g(2) * t242) * t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210, t211, -t197, -qJ(2) * t197 - t179, 0, 0, 0, 0, 0, 0, t164, -t37, t268, t263 * t65 + t66 * t98 + t169, 0, 0, 0, 0, 0, 0, t164, t268, t37, t250 + t57 * t98 + (-qJD(4) - t56) * t263 + t165, 0, 0, 0, 0, 0, 0, t167 - t276, -t177 - t266, t245 + t246, t11 * t48 - t12 * t185 + t165 - t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t239, t267, t38, -t239, -t188, qJDD(3), -t113 * t263 + t175 + t221, t113 * t98 + t256, 0, 0, t239, t38, -t267, qJDD(3), t188, -t239, -t58 * t98 - t168 + t221 + 0.2e1 * t229, pkin(3) * t60 - t61 * qJ(4) + (t57 - t66) * t263 + (t56 - t220) * t98, t263 * t58 - t43 * t98 + 0.2e1 * t150 + 0.2e1 * t151 - t256, t20 * qJ(4) - t21 * pkin(3) - t43 * t58 - t56 * t66 - g(1) * (-pkin(3) * t226 + t121) - g(2) * (-pkin(3) * t227 + t119) - g(3) * t219 + t220 * t57, t244, -t277, -t272, -t244, -t278, t145, -t114 * t145 - t149 * t240 - t185 * t36 - t2 - t270, t115 * t145 + t149 * t241 + t36 * t48 - t260, -t114 * t177 + t115 * t167 + (t11 - t241) * t185 + (t12 + t240) * t48, t1 * t115 + t2 * t114 - t28 * t36 - g(1) * t121 - g(2) * t119 - g(3) * (t219 + t247) + t241 * t12 + t240 * t11 + t264 * t140 * t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) + t239, t38, -qJD(3) ^ 2 - t255, -t57 * qJD(3) + t168 - t229, 0, 0, 0, 0, 0, 0, -t159 * t145 - t156 * t196 - t185 * t263, t156 * t145 - t159 * t196 + t263 * t48, t278 * t156 - t272 * t159, -t28 * t263 + t269 * t159 + (t1 + t238) * t156 + t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244, t277, t272, t244, t278, -t145, t269 + t270, -t238 + t260, 0, 0;];
tau_reg = t6;
