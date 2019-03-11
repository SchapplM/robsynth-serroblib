% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:56:17
% EndTime: 2019-03-09 01:56:23
% DurationCPUTime: 3.40s
% Computational Cost: add. (5159->423), mult. (10126->487), div. (0->0), fcn. (7056->10), ass. (0->220)
t138 = sin(pkin(9));
t269 = sin(qJ(4));
t210 = t269 * t138;
t111 = qJD(1) * t210;
t139 = cos(pkin(9));
t145 = cos(qJ(4));
t94 = t138 * t145 + t139 * t269;
t154 = -qJD(4) * t111 + qJDD(1) * t94;
t232 = t145 * t139;
t209 = qJD(1) * t232;
t86 = -t111 + t209;
t297 = t154 + (-t86 + t209) * qJD(4);
t296 = t154 + (t86 + t209) * qJD(4);
t88 = t94 * qJD(4);
t95 = -t210 + t232;
t251 = -qJD(4) * t88 + qJDD(4) * t95;
t84 = t94 * qJD(1);
t295 = qJD(1) * t84 - t251;
t144 = cos(qJ(6));
t75 = qJD(6) + t86;
t198 = t144 * t75;
t142 = sin(qJ(6));
t60 = qJD(4) * t144 + t142 * t84;
t294 = t60 * t198;
t134 = pkin(9) + qJ(4);
t123 = sin(t134);
t124 = cos(t134);
t146 = cos(qJ(1));
t130 = g(2) * t146;
t143 = sin(qJ(1));
t284 = g(1) * t143 - t130;
t158 = -g(3) * t124 - t123 * t284;
t208 = qJD(4) * t269;
t222 = qJD(4) * t145;
t141 = -pkin(1) - qJ(3);
t278 = -qJD(1) * qJD(3) + qJDD(1) * t141;
t96 = qJDD(2) + t278;
t201 = -pkin(7) * qJDD(1) + t96;
t65 = t201 * t138;
t66 = t201 * t139;
t108 = qJD(1) * t141 + qJD(2);
t202 = -pkin(7) * qJD(1) + t108;
t70 = t202 * t138;
t71 = t202 * t139;
t204 = -t145 * t65 + t208 * t70 - t222 * t71 - t269 * t66;
t293 = t158 - t204;
t189 = g(1) * t146 + g(2) * t143;
t166 = t189 * t123;
t252 = -pkin(7) + t141;
t97 = t252 * t138;
t98 = t252 * t139;
t53 = t145 * t97 + t269 * t98;
t37 = qJD(3) * t95 + qJD(4) * t53;
t52 = -t145 * t98 + t269 * t97;
t292 = -t37 * qJD(4) - t52 * qJDD(4) - t166;
t42 = -t145 * t71 + t269 * t70;
t228 = -qJD(5) - t42;
t227 = pkin(5) * t86 - t228;
t273 = pkin(4) + pkin(8);
t24 = -qJD(4) * t273 + t227;
t122 = qJD(1) * qJ(2) + qJD(3);
t126 = t138 * pkin(3);
t99 = qJD(1) * t126 + t122;
t170 = -t86 * qJ(5) + t99;
t28 = t273 * t84 + t170;
t10 = t142 * t24 + t144 * t28;
t203 = -t145 * t66 + t208 * t71 + t222 * t70 + t269 * t65;
t178 = qJDD(5) + t203;
t215 = t139 * qJDD(1);
t167 = -qJDD(1) * t210 + t145 * t215;
t50 = qJD(1) * t88 - t167;
t6 = -t50 * pkin(5) - qJDD(4) * t273 + t178;
t245 = t50 * qJ(5);
t135 = qJDD(1) * qJ(2);
t136 = qJD(1) * qJD(2);
t283 = t135 + t136;
t101 = qJDD(3) + t283;
t216 = t138 * qJDD(1);
t119 = pkin(3) * t216;
t92 = t119 + t101;
t155 = -t86 * qJD(5) + t245 + t92;
t51 = qJD(4) * t209 + t154;
t8 = t273 * t51 + t155;
t2 = -qJD(6) * t10 - t142 * t8 + t144 * t6;
t291 = t10 * t75 + t2;
t220 = t142 * qJD(4);
t58 = -t144 * t84 + t220;
t206 = t58 * t75;
t221 = qJD(6) * t144;
t22 = qJD(6) * t220 - qJDD(4) * t144 - t142 * t51 - t221 * t84;
t290 = t22 - t206;
t23 = qJD(6) * t60 + t142 * qJDD(4) - t144 * t51;
t289 = t60 * t75 - t23;
t275 = t84 ^ 2;
t82 = t86 ^ 2;
t288 = -t275 - t82;
t287 = -t275 + t82;
t132 = t138 ^ 2;
t133 = t139 ^ 2;
t223 = t132 + t133;
t286 = t108 * t223;
t140 = -pkin(7) - qJ(3);
t285 = t126 * t146 + t140 * t143;
t160 = -t189 + t101;
t271 = t84 * pkin(5);
t43 = t145 * t70 + t269 * t71;
t38 = -qJD(4) * qJ(5) - t43;
t27 = -t38 - t271;
t49 = -qJDD(6) + t50;
t280 = t27 * t75 + t273 * t49;
t36 = qJD(3) * t94 + t208 * t97 - t222 * t98;
t279 = t36 * qJD(4) - t53 * qJDD(4) - t124 * t189;
t217 = qJDD(4) * qJ(5);
t13 = -qJD(4) * qJD(5) + t204 - t217;
t238 = qJDD(4) * pkin(4);
t15 = t178 - t238;
t35 = -qJD(4) * pkin(4) - t228;
t89 = -t138 * t208 + t139 * t222;
t277 = t13 * t94 + t15 * t95 - t35 * t88 + t38 * t89;
t276 = t203 * t95 + t204 * t94 - t42 * t88 - t43 * t89;
t274 = 0.2e1 * t136;
t272 = t51 * pkin(4);
t9 = -t142 * t28 + t144 * t24;
t270 = t9 * t75;
t117 = g(3) * t123;
t266 = t123 * pkin(4);
t257 = t60 * t58;
t255 = t60 * t84;
t254 = t84 * t58;
t253 = t84 * t86;
t249 = t142 * t49;
t248 = t144 * t22;
t44 = t144 * t49;
t246 = t23 * t142;
t244 = t84 * qJ(5);
t243 = pkin(1) * qJDD(1);
t242 = qJ(5) * t123;
t240 = qJD(4) * t84;
t237 = t123 * t143;
t115 = t124 * qJ(5);
t236 = t124 * t143;
t235 = t124 * t146;
t234 = t143 * t142;
t233 = t143 * t144;
t231 = t146 * t142;
t230 = t146 * t144;
t229 = t43 * qJD(4);
t116 = qJ(2) + t126;
t226 = pkin(4) * t236 + qJ(5) * t237;
t225 = pkin(1) * t146 + qJ(2) * t143;
t214 = qJD(6) * t142 * t94;
t213 = t94 * t221;
t212 = g(1) * t236 - g(2) * t235 - t117;
t211 = t126 * t143 + t225;
t128 = t146 * qJ(2);
t207 = -t143 * pkin(1) + t128;
t205 = t223 * t96;
t200 = t142 * t75;
t197 = qJD(6) * t95 + qJD(1);
t195 = qJDD(2) - t243;
t194 = -t2 * t95 + t9 * t88;
t1 = qJD(6) * t9 + t142 * t6 + t144 * t8;
t193 = -t1 * t95 + t10 * t88;
t7 = -pkin(5) * t51 - t13;
t192 = t27 * t89 + t7 * t94;
t191 = -t95 * qJ(5) + t116;
t190 = t1 - t270;
t185 = -t22 * t95 - t60 * t88;
t184 = -t94 * t22 + t89 * t60;
t183 = t23 * t95 - t58 * t88;
t182 = -t94 * t23 - t89 * t58;
t181 = -t49 * t94 + t75 * t89;
t180 = -t50 * t95 - t86 * t88;
t179 = t51 * t94 + t84 * t89;
t177 = t10 * t144 - t142 * t9;
t176 = pkin(8) * t123 - t115;
t34 = t273 * t94 + t191;
t40 = pkin(5) * t95 + t52;
t17 = t142 * t40 + t144 * t34;
t16 = -t142 * t34 + t144 * t40;
t173 = qJD(4) * t89 + qJDD(4) * t94;
t172 = t207 + t285;
t169 = -t146 * t140 + t211;
t168 = -t200 * t75 - t44;
t165 = t88 * qJ(5) - t95 * qJD(5) + qJD(2);
t164 = -t203 - t212;
t162 = -t198 * t75 + t249;
t161 = qJD(1) * t86 + t173;
t159 = t119 + t160;
t157 = -t179 - t180;
t156 = t50 * t94 - t51 * t95 + t84 * t88 - t86 * t89;
t39 = pkin(4) * t84 + t170;
t152 = t39 * t86 + qJDD(5) - t164;
t151 = t36 * t84 + t37 * t86 - t50 * t52 - t51 * t53 + t284;
t150 = t160 + t283;
t149 = qJD(6) * t273 * t75 + t158 + t7;
t148 = qJD(1) ^ 2;
t105 = t146 * t266;
t103 = pkin(4) * t237;
t79 = -t124 * t234 + t230;
t78 = -t124 * t233 - t231;
t77 = -t124 * t231 - t233;
t76 = -t124 * t230 + t234;
t48 = pkin(4) * t86 + t244;
t47 = pkin(4) * t94 + t191;
t41 = -pkin(5) * t94 + t53;
t33 = -t240 + t50;
t32 = pkin(4) * t89 + t165;
t31 = t273 * t86 + t244;
t30 = t43 - t271;
t26 = -pkin(5) * t88 + t37;
t25 = -t89 * pkin(5) - t36;
t21 = t273 * t89 + t165;
t20 = t144 * t23;
t14 = t155 + t272;
t12 = t142 * t30 + t144 * t31;
t11 = -t142 * t31 + t144 * t30;
t4 = -qJD(6) * t17 - t142 * t21 + t144 * t26;
t3 = qJD(6) * t16 + t142 * t26 + t144 * t21;
t5 = [0, 0, 0, 0, 0, qJDD(1), t284, t189, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t284 - 0.2e1 * t243, 0.2e1 * t135 + t274 - t189, -t195 * pkin(1) - g(1) * t207 - g(2) * t225 + (t135 + t274) * qJ(2), t133 * qJDD(1), -0.2e1 * t138 * t215, 0, t132 * qJDD(1), 0, 0, t150 * t138, t150 * t139, t284 + t223 * (-t278 - t96) t101 * qJ(2) + t122 * qJD(2) - g(1) * (t141 * t143 + t128) - g(2) * (qJ(3) * t146 + t225) + t141 * t205 - qJD(3) * t286, t180, t156, t251, t179, -t173, 0, qJD(2) * t84 + t116 * t51 + t99 * t89 + t92 * t94 + t292, qJD(2) * t86 - t116 * t50 - t99 * t88 + t92 * t95 + t279, t151 + t276, -g(1) * t172 - g(2) * t169 + qJD(2) * t99 + t116 * t92 + t203 * t52 - t204 * t53 - t36 * t43 + t37 * t42, 0, -t251, t173, t180, t156, t179, t151 + t277, -t14 * t94 - t32 * t84 - t39 * t89 - t47 * t51 - t292, -t14 * t95 - t32 * t86 + t39 * t88 + t47 * t50 - t279, t14 * t47 + t39 * t32 - t13 * t53 + t38 * t36 + t15 * t52 + t35 * t37 - g(1) * (-qJ(5) * t235 + t105 + t172) - g(2) * (-qJ(5) * t236 + t103 + t169) t142 * t184 + t213 * t60 (-t142 * t58 + t144 * t60) * t89 + (-t246 - t248 + (-t142 * t60 - t144 * t58) * qJD(6)) * t94, t142 * t181 + t213 * t75 + t185, t144 * t182 + t214 * t58, t144 * t181 - t214 * t75 - t183, -t49 * t95 - t75 * t88, -g(1) * t77 - g(2) * t79 - t144 * t192 - t16 * t49 + t214 * t27 + t41 * t23 + t25 * t58 + t4 * t75 - t194, -g(1) * t76 - g(2) * t78 + t142 * t192 + t17 * t49 + t213 * t27 - t41 * t22 + t25 * t60 - t3 * t75 + t193, t16 * t22 - t17 * t23 - t3 * t58 - t4 * t60 + t177 * t89 - t166 + (t1 * t144 - t2 * t142 + (-t10 * t142 - t144 * t9) * qJD(6)) * t94, t1 * t17 + t10 * t3 + t2 * t16 + t9 * t4 + t7 * t41 + t27 * t25 - g(1) * (t105 + t128 + t285) - g(2) * (t103 + t211) + (-g(1) * t176 - g(2) * (pkin(5) - t140)) * t146 + (-g(1) * (-pkin(1) - pkin(5)) - g(2) * t176) * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t148, -qJ(2) * t148 + t195 - t284, 0, 0, 0, 0, 0, 0, -t148 * t138, -t148 * t139, -t223 * qJDD(1), -t122 * qJD(1) + t205 - t284, 0, 0, 0, 0, 0, 0, -t295, -t161, t157, -qJD(1) * t99 - t276 - t284, 0, 0, 0, 0, 0, 0, t157, t295, t161, -qJD(1) * t39 - t277 - t284, 0, 0, 0, 0, 0, 0, t95 * t44 + (t142 * t197 + t144 * t88) * t75 - t182, -t95 * t249 + (-t142 * t88 + t144 * t197) * t75 + t184 (t197 * t58 + t185) * t144 + (-t197 * t60 + t183) * t142 (-t10 * t197 + t194) * t144 + (t197 * t9 + t193) * t142 + t192 - t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216, t215, -t223 * t148, qJD(1) * t286 + t160, 0, 0, 0, 0, 0, 0, t296, t167 - 0.2e1 * t240, t288, -t42 * t86 + t43 * t84 + t159, 0, 0, 0, 0, 0, 0, t288, -t296, t240 + t50, t272 + t245 - t38 * t84 + (-qJD(5) - t35) * t86 + t159, 0, 0, 0, 0, 0, 0, t162 + t254, t142 * t75 ^ 2 + t255 + t44, -t142 * t290 - t20 + t294, -t142 * t291 + t144 * t190 + t27 * t84 - t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, t287, t167, -t253, -t297, qJDD(4), -t86 * t99 + t164 + t229, -t42 * qJD(4) + t99 * t84 - t293, 0, 0, qJDD(4), t33, t297, t253, t287, -t253, pkin(4) * t50 - qJ(5) * t51 + (-t38 - t43) * t86 + (t35 + t228) * t84, t48 * t84 + t152 - t229 - 0.2e1 * t238, 0.2e1 * t217 - t39 * t84 + t48 * t86 + (0.2e1 * qJD(5) + t42) * qJD(4) + t293, -t13 * qJ(5) - t15 * pkin(4) - t39 * t48 - t35 * t43 - g(1) * t226 - g(3) * (t115 - t266) + t228 * t38 - (-pkin(4) * t124 - t242) * t130, -t200 * t60 - t248, -t20 - t294 + (t22 + t206) * t142, t168 + t255, t198 * t58 + t246, t162 - t254, t75 * t84, qJ(5) * t23 - t11 * t75 + t149 * t142 + t144 * t280 + t227 * t58 + t9 * t84, -qJ(5) * t22 - t10 * t84 + t12 * t75 - t142 * t280 + t149 * t144 + t227 * t60, t11 * t60 + t12 * t58 + (-t10 * t86 - t273 * t22 - t2 + (t273 * t58 - t10) * qJD(6)) * t144 + (t273 * t23 + t86 * t9 - t1 + (-t273 * t60 + t9) * qJD(6)) * t142 - t212, t7 * qJ(5) - t10 * t12 - t9 * t11 - g(1) * (pkin(8) * t236 + t226) - g(3) * (-t123 * t273 + t115) + t227 * t27 - (-t124 * t273 - t242) * t130 - (qJD(6) * t177 + t1 * t142 + t2 * t144) * t273; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, qJDD(4) - t253, -qJD(4) ^ 2 - t82, qJD(4) * t38 + t152 - t238, 0, 0, 0, 0, 0, 0, -qJD(4) * t58 + t168, -qJD(4) * t60 + t162, t142 * t289 + t144 * t290, -t27 * qJD(4) + t142 * t190 + t144 * t291 + t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t257, -t58 ^ 2 + t60 ^ 2, -t290, -t257, t289, -t49, -g(1) * t78 + g(2) * t76 - t117 * t144 - t27 * t60 + t291, g(1) * t79 - g(2) * t77 + t27 * t58 + t270 + (-qJD(6) * t24 - t8) * t144 + (qJD(6) * t28 + t117 - t6) * t142, 0, 0;];
tau_reg  = t5;
