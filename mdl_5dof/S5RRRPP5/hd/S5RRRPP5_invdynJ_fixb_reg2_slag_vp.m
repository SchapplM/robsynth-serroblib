% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:34
% EndTime: 2019-12-31 20:58:38
% DurationCPUTime: 2.24s
% Computational Cost: add. (2963->372), mult. (6771->399), div. (0->0), fcn. (4532->8), ass. (0->194)
t161 = cos(qJ(2));
t258 = t161 * pkin(2);
t139 = pkin(1) + t258;
t103 = t139 * qJD(1);
t159 = sin(qJ(2));
t158 = sin(qJ(3));
t236 = t158 * t161;
t262 = cos(qJ(3));
t94 = t262 * t159 + t236;
t82 = t94 * qJD(1);
t278 = -t82 * qJ(4) - t103;
t152 = qJDD(2) + qJDD(3);
t153 = qJD(2) + qJD(3);
t273 = t152 * qJ(4) + t153 * qJD(4);
t223 = t161 * qJDD(1);
t225 = qJD(1) * qJD(2);
t277 = t159 * t225 - t223;
t264 = pkin(7) + pkin(6);
t210 = qJDD(1) * t262;
t224 = t159 * qJDD(1);
t195 = t158 * t224 - t161 * t210;
t59 = t153 * t94;
t38 = t59 * qJD(1) + t195;
t217 = t262 * t161;
t205 = qJD(1) * t217;
t227 = qJD(1) * t159;
t215 = t158 * t227;
t80 = -t205 + t215;
t276 = t38 * qJ(5) + t80 * qJD(5);
t246 = t82 * t153;
t105 = t264 * t161;
t100 = qJD(1) * t105;
t238 = t158 * t100;
t254 = qJD(2) * pkin(2);
t104 = t264 * t159;
t98 = qJD(1) * t104;
t91 = -t98 + t254;
t52 = t262 * t91 - t238;
t275 = qJD(4) - t52;
t274 = 0.2e1 * t273;
t156 = qJ(2) + qJ(3);
t146 = sin(t156);
t147 = cos(t156);
t231 = t147 * pkin(3) + t146 * qJ(4);
t160 = sin(qJ(1));
t162 = cos(qJ(1));
t211 = g(1) * t160 - g(2) * t162;
t272 = g(1) * t162 + g(2) * t160;
t266 = t82 ^ 2;
t271 = -t153 ^ 2 - t266;
t267 = t80 ^ 2;
t43 = t266 - t267;
t144 = t152 * pkin(3);
t270 = qJDD(4) - t144;
t237 = t158 * t159;
t196 = t153 * t237;
t208 = -t153 * t205 - t158 * t223 - t159 * t210;
t37 = qJD(1) * t196 + t208;
t269 = -t152 * pkin(4) + t37 * qJ(5) - t82 * qJD(5);
t214 = t262 * qJD(3);
t120 = pkin(2) * t214 + qJD(4);
t132 = t158 * pkin(2) + qJ(4);
t268 = t120 * t153 + t132 * t152 + t273;
t265 = pkin(3) + pkin(4);
t263 = t82 * pkin(4);
t261 = pkin(2) * t159;
t259 = g(3) * t161;
t133 = t147 * pkin(4);
t44 = t80 * pkin(3) + t278;
t257 = t44 * t80;
t256 = t82 * t80;
t50 = t82 * pkin(3) + t80 * qJ(4);
t57 = -t262 * t98 - t238;
t65 = -t158 * t104 + t262 * t105;
t255 = qJ(4) * t38;
t253 = t103 * t80;
t252 = t132 * t38;
t251 = t153 * t80;
t89 = t262 * t100;
t53 = t158 * t91 + t89;
t250 = t53 * t153;
t56 = -t158 * t98 + t89;
t249 = t56 * t153;
t248 = t80 * qJ(5);
t73 = t82 * qJ(5);
t41 = t73 + t57;
t245 = t120 - t41;
t244 = t120 - t57;
t243 = pkin(6) * qJDD(1);
t157 = qJDD(1) * pkin(1);
t242 = t146 * t160;
t241 = t146 * t162;
t240 = t147 * t160;
t239 = t147 * t162;
t235 = t162 * t264;
t234 = qJ(5) - t264;
t35 = t73 + t52;
t233 = qJD(4) - t35;
t154 = t159 ^ 2;
t155 = t161 ^ 2;
t229 = t154 - t155;
t228 = t154 + t155;
t226 = qJD(3) * t158;
t222 = pkin(2) * t226;
t221 = t159 * t254;
t166 = qJD(1) ^ 2;
t219 = t159 * t166 * t161;
t75 = t277 * pkin(2) - t157;
t218 = qJD(2) * t264;
t216 = t153 * t226;
t212 = t161 * t225;
t62 = qJDD(2) * pkin(2) + t264 * (-t212 - t224);
t63 = t264 * t277;
t13 = -t100 * t226 + t158 * t62 + t91 * t214 - t262 * t63;
t209 = t100 * t214 - t158 * t63 + t91 * t226 - t262 * t62;
t64 = t262 * t104 + t158 * t105;
t138 = -t262 * pkin(2) - pkin(3);
t110 = t162 * t139;
t207 = g(2) * (pkin(3) * t239 + qJ(4) * t241 + t110);
t206 = qJD(2) * t217;
t204 = t159 * t212;
t40 = t56 + t248;
t203 = -t40 + t222;
t202 = -t56 + t222;
t201 = g(1) * t242 - g(2) * t241;
t200 = g(1) * t240 - g(2) * t239;
t199 = t231 + t258;
t198 = -pkin(3) * t146 - t261;
t93 = -t217 + t237;
t8 = t38 * t93 + t80 * t59;
t58 = -t161 * t214 + t196 - t206;
t33 = t94 * t152 - t58 * t153;
t194 = t152 * t93 + t153 * t59;
t193 = -t152 + t256;
t46 = pkin(2) * t227 + t50;
t9 = t13 + t273;
t192 = t94 * qJ(4) + t139;
t5 = t38 * pkin(3) + t37 * qJ(4) - t82 * qJD(4) + t75;
t10 = t209 + t270;
t191 = -t139 - t231;
t190 = -0.2e1 * pkin(1) * t225 - pkin(6) * qJDD(2);
t99 = t159 * t218;
t25 = -t104 * t214 - t105 * t226 - t218 * t236 - t262 * t99;
t189 = g(1) * t241 + g(2) * t242 - g(3) * t147 - t209;
t188 = -g(1) * t239 - g(2) * t240 - g(3) * t146 + t13;
t26 = t65 * qJD(3) - t158 * t99 + t264 * t206;
t185 = -t25 * t80 + t26 * t82 - t64 * t37 - t65 * t38 - t272;
t184 = -t58 * qJ(4) + t94 * qJD(4) - t221;
t36 = t53 + t248;
t183 = t189 - t270;
t182 = t103 * t82 + t189;
t181 = t52 * t153 - t188;
t180 = t57 * t153 - t188;
t2 = -t38 * pkin(4) + qJDD(5) - t5;
t179 = t65 * t152 + t25 * t153 + t201;
t178 = -t64 * t152 - t26 * t153 + t200;
t177 = t37 * t93 - t94 * t38 + t58 * t80 - t82 * t59;
t165 = qJD(2) ^ 2;
t176 = -pkin(6) * t165 + 0.2e1 * t157 + t211;
t175 = pkin(1) * t166 - t243 + t272;
t174 = t272 * t265 * t146;
t24 = -t265 * t80 + qJD(5) - t278;
t173 = t24 * t80 + t188 + t276;
t172 = -t153 * t215 - t208;
t171 = -t44 * t82 + t183;
t18 = t37 - t251;
t170 = -t24 * t82 - t183 + t269;
t169 = -t195 - t246;
t168 = t169 + t246;
t145 = t153 * qJ(4);
t128 = -pkin(4) + t138;
t109 = qJ(4) * t239;
t107 = qJ(4) * t240;
t106 = pkin(2) * t216;
t51 = t93 * pkin(3) - t192;
t49 = t145 + t53;
t48 = t93 * qJ(5) + t65;
t47 = -t94 * qJ(5) + t64;
t45 = -t153 * pkin(3) + t275;
t39 = -t265 * t93 + t192;
t34 = -t50 - t263;
t28 = t145 + t36;
t27 = -t46 - t263;
t23 = -t265 * t153 + t275 - t73;
t20 = t38 - t246;
t19 = t172 + t251;
t15 = t59 * pkin(3) - t184;
t12 = t58 * qJ(5) - t94 * qJD(5) + t26;
t11 = t59 * qJ(5) + t93 * qJD(5) + t25;
t7 = -t37 * t94 - t82 * t58;
t6 = -t265 * t59 + t184;
t4 = t9 + t276;
t3 = t10 + t269;
t1 = [0, 0, 0, 0, 0, qJDD(1), t211, t272, 0, 0, t154 * qJDD(1) + 0.2e1 * t204, 0.2e1 * t159 * t223 - 0.2e1 * t229 * t225, qJDD(2) * t159 + t165 * t161, t155 * qJDD(1) - 0.2e1 * t204, qJDD(2) * t161 - t165 * t159, 0, t159 * t190 + t161 * t176, -t159 * t176 + t161 * t190, 0.2e1 * t228 * t243 - t272, -g(1) * (-t160 * pkin(1) + t162 * pkin(6)) - g(2) * (t162 * pkin(1) + t160 * pkin(6)) + (pkin(6) ^ 2 * t228 + pkin(1) ^ 2) * qJDD(1), t7, t177, t33, t8, -t194, 0, -t103 * t59 - t139 * t38 + t221 * t80 + t75 * t93 + t178, t103 * t58 + t139 * t37 + t221 * t82 + t75 * t94 - t179, -t13 * t93 + t209 * t94 + t52 * t58 - t53 * t59 + t185, t13 * t65 + t53 * t25 + t209 * t64 - t52 * t26 - t75 * t139 - t103 * t221 - g(1) * (-t160 * t139 + t235) - g(2) * (t160 * t264 + t110), t7, t33, -t177, 0, t194, t8, t15 * t80 + t51 * t38 + t44 * t59 + t5 * t93 + t178, t10 * t94 - t45 * t58 - t49 * t59 - t9 * t93 + t185, -t15 * t82 + t51 * t37 + t44 * t58 - t5 * t94 + t179, t9 * t65 + t49 * t25 + t5 * t51 + t44 * t15 + t10 * t64 + t45 * t26 - g(1) * t235 - t207 + (-g(1) * t191 - g(2) * t264) * t160, t7, -t177, -t33, t8, -t194, 0, -t12 * t153 - t47 * t152 - t2 * t93 - t24 * t59 - t39 * t38 - t6 * t80 + t200, t11 * t153 + t48 * t152 + t2 * t94 - t24 * t58 - t39 * t37 + t6 * t82 + t201, t11 * t80 - t12 * t82 + t23 * t58 + t28 * t59 - t3 * t94 + t47 * t37 + t48 * t38 + t4 * t93 + t272, t4 * t48 + t28 * t11 + t3 * t47 + t23 * t12 + t2 * t39 + t24 * t6 - t207 + (g(1) * t234 - g(2) * t133) * t162 + (-g(1) * (t191 - t133) + g(2) * t234) * t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t219, t229 * t166, t224, t219, t223, qJDD(2), t159 * t175 - t259, g(3) * t159 + t161 * t175, 0, 0, t256, t43, t19, -t256, t168, t152, t249 + (t262 * t152 - t80 * t227 - t216) * pkin(2) + t182, -t253 + (-t152 * t158 - t153 * t214 - t227 * t82) * pkin(2) + t180, (t53 - t56) * t82 + (-t52 + t57) * t80 + (t262 * t37 - t158 * t38 + (t158 * t82 - t262 * t80) * qJD(3)) * pkin(2), t52 * t56 - t53 * t57 + (-t262 * t209 - t259 + t13 * t158 + (-t158 * t52 + t262 * t53) * qJD(3) + (qJD(1) * t103 + t272) * t159) * pkin(2), t256, t19, -t43, t152, t20, -t256, -t138 * t152 - t46 * t80 - t106 + t171 + t249, -t252 - t138 * t37 + (t202 + t49) * t82 + (t45 - t244) * t80, t46 * t82 - t180 - t257 + t268, t9 * t132 + t10 * t138 - t44 * t46 - g(1) * (t162 * t198 + t109) - g(2) * (t160 * t198 + t107) - g(3) * t199 + t244 * t49 + t202 * t45, t256, -t43, t18, -t256, t168, t152, -t128 * t152 + t40 * t153 + t27 * t80 - t106 - t170, -t41 * t153 - t27 * t82 + t173 + t268, t128 * t37 + t252 + (-t203 - t28) * t82 + (-t23 + t245) * t80, t4 * t132 + t3 * t128 - t24 * t27 - g(1) * (-t162 * t261 + t109) - g(2) * (-t160 * t261 + t107) - g(3) * (t133 + t199) + t245 * t28 + t203 * t23 + t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256, t43, t19, -t256, t168, t152, t182 + t250, t181 - t253, 0, 0, t256, t19, -t43, t152, t20, -t256, -t50 * t80 + t144 + t171 + t250, pkin(3) * t37 - t255 + (t49 - t53) * t82 + (t45 - t275) * t80, t50 * t82 - t181 - t257 + t274, t9 * qJ(4) - t10 * pkin(3) - t44 * t50 - t45 * t53 - g(1) * (-pkin(3) * t241 + t109) - g(2) * (-pkin(3) * t242 + t107) - g(3) * t231 + t275 * t49, t256, -t43, t18, -t256, t168, t152, t152 * t265 + t36 * t153 + t34 * t80 - t170, -t35 * t153 - t34 * t82 + t173 + t274, t255 - t265 * t37 + (-t28 + t36) * t82 + (-t23 + t233) * t80, t4 * qJ(4) - t3 * t265 - t23 * t36 - t24 * t34 - g(1) * t109 - g(2) * t107 - g(3) * (t133 + t231) + t233 * t28 + t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t19, t271, -t49 * t153 - t171, 0, 0, 0, 0, 0, 0, t193, t271, t18, -t28 * t153 + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169 - t246, t172 - t251, -t266 - t267, t23 * t82 - t28 * t80 + t2 + t211;];
tau_reg = t1;
