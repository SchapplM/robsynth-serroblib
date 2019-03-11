% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:56:07
% EndTime: 2019-03-09 03:56:16
% DurationCPUTime: 3.53s
% Computational Cost: add. (5069->357), mult. (10527->456), div. (0->0), fcn. (7979->12), ass. (0->206)
t189 = cos(qJ(6));
t247 = qJD(6) * t189;
t182 = sin(pkin(10));
t183 = cos(pkin(10));
t187 = sin(qJ(3));
t191 = cos(qJ(3));
t213 = t182 * t191 + t183 * t187;
t123 = t213 * qJD(1);
t190 = cos(qJ(5));
t107 = t190 * t123;
t252 = qJD(1) * t191;
t253 = qJD(1) * t187;
t126 = -t182 * t253 + t183 * t252;
t186 = sin(qJ(5));
t75 = t126 * t186 + t107;
t311 = t189 * t75;
t321 = t247 + t311;
t176 = qJD(3) + qJD(5);
t279 = t176 * t75;
t249 = qJD(5) * t186;
t212 = t182 * t187 - t183 * t191;
t245 = qJD(1) * qJD(3);
t82 = -qJDD(1) * t213 + t212 * t245;
t241 = t191 * qJDD(1);
t242 = t187 * qJDD(1);
t83 = t182 * t242 - t183 * t241 + t213 * t245;
t31 = -qJD(5) * t107 - t126 * t249 + t186 * t82 - t190 * t83;
t320 = t31 + t279;
t195 = qJD(1) ^ 2;
t192 = cos(qJ(1));
t173 = g(2) * t192;
t188 = sin(qJ(1));
t174 = g(1) * t188;
t257 = t174 - t173;
t203 = -qJ(2) * t195 - t257;
t175 = qJDD(3) + qJDD(5);
t185 = sin(qJ(6));
t215 = -t123 * t186 + t190 * t126;
t248 = qJD(6) * t185;
t13 = t185 * t175 + t176 * t247 + t189 * t31 - t215 * t248;
t250 = qJD(3) * t191;
t251 = qJD(3) * t187;
t124 = t182 * t251 - t183 * t250;
t125 = t213 * qJD(3);
t214 = -t186 * t212 + t190 * t213;
t45 = qJD(5) * t214 - t186 * t124 + t190 * t125;
t60 = t176 * t185 + t189 * t215;
t86 = -t186 * t213 - t190 * t212;
t319 = t13 * t86 - t45 * t60;
t318 = t175 * t86 - t176 * t45;
t14 = qJD(6) * t60 - t189 * t175 + t185 * t31;
t58 = -t189 * t176 + t185 * t215;
t317 = t13 * t189 - t185 * t14 - t321 * t58;
t11 = t13 * t185;
t316 = t321 * t60 + t11;
t32 = qJD(5) * t215 - t186 * t83 - t190 * t82;
t30 = qJDD(6) + t32;
t27 = t185 * t30;
t282 = t60 * t215;
t309 = qJD(6) + t75;
t65 = t309 * t247;
t315 = t309 * t311 + t27 - t282 + t65;
t287 = pkin(8) * t126;
t193 = -pkin(1) - pkin(7);
t146 = t193 * qJD(1) + qJD(2);
t120 = -qJ(4) * t252 + t191 * t146;
t102 = qJD(3) * pkin(3) + t120;
t119 = -qJ(4) * t253 + t146 * t187;
t98 = t182 * t119;
t63 = t183 * t102 - t98;
t47 = qJD(3) * pkin(4) - t287 + t63;
t288 = pkin(8) * t123;
t268 = t183 * t119;
t64 = t182 * t102 + t268;
t48 = t64 - t288;
t22 = -t186 * t48 + t190 * t47;
t18 = -pkin(5) * t176 - t22;
t313 = t18 * t75;
t168 = qJ(3) + pkin(10) + qJ(5);
t160 = cos(t168);
t238 = t160 * t174;
t144 = t193 * qJDD(1) + qJDD(2);
t135 = t191 * t144;
t244 = qJD(1) * qJD(4);
t62 = -t191 * t244 - t146 * t251 + qJDD(3) * pkin(3) + t135 + (t187 * t245 - t241) * qJ(4);
t68 = (-qJ(4) * qJD(1) + t146) * t250 + (-qJ(4) * qJDD(1) + t144 - t244) * t187;
t35 = -t182 * t68 + t183 * t62;
t21 = qJDD(3) * pkin(4) + pkin(8) * t83 + t35;
t23 = t186 * t47 + t190 * t48;
t36 = t182 * t62 + t183 * t68;
t24 = pkin(8) * t82 + t36;
t294 = qJD(5) * t23 + t186 * t24 - t190 * t21;
t3 = -t175 * pkin(5) + t294;
t219 = -t3 - t238;
t310 = t215 * t75;
t159 = sin(t168);
t308 = g(3) * t159 + t160 * t173;
t275 = t215 * t176;
t307 = -t32 + t275;
t305 = t215 ^ 2 - t75 ^ 2;
t151 = g(3) * t160;
t293 = (qJD(5) * t47 + t24) * t190 + t186 * t21 - t48 * t249;
t137 = pkin(3) * t253 + qJD(1) * qJ(2) + qJD(4);
t90 = pkin(4) * t123 + t137;
t304 = t257 * t159 + t75 * t90 + t151 - t293;
t302 = pkin(5) * t215;
t161 = pkin(3) * t183 + pkin(4);
t289 = pkin(3) * t182;
t259 = t186 * t161 + t190 * t289;
t118 = pkin(9) + t259;
t94 = pkin(3) * t252 + pkin(4) * t126;
t301 = (pkin(9) * t75 + qJD(6) * t118 + t302 + t94) * t309;
t300 = (t309 * pkin(9) + t302) * t309;
t283 = t58 * t215;
t299 = t309 * t215;
t177 = qJDD(1) * qJ(2);
t218 = g(1) * t192 + g(2) * t188;
t178 = qJD(1) * qJD(2);
t237 = 0.2e1 * t178;
t298 = 0.2e1 * t177 + t237 - t218;
t19 = pkin(9) * t176 + t23;
t33 = pkin(5) * t75 - pkin(9) * t215 + t90;
t6 = -t185 * t19 + t189 * t33;
t297 = t18 * t248 + t308 * t189 - t6 * t215;
t7 = t185 * t33 + t189 * t19;
t296 = t18 * t247 - t219 * t185 + t7 * t215;
t295 = -t215 * t90 - t238 - t294 + t308;
t199 = t86 * qJD(5) - t190 * t124 - t186 * t125;
t292 = -t175 * t214 - t176 * t199;
t231 = t175 * pkin(9) + qJD(6) * t33 + t293;
t262 = qJ(4) - t193;
t138 = t262 * t187;
t139 = t262 * t191;
t88 = t138 * t182 - t183 * t139;
t69 = pkin(8) * t212 + t88;
t89 = -t183 * t138 - t182 * t139;
t70 = -pkin(8) * t213 + t89;
t38 = t186 * t69 + t190 * t70;
t171 = t187 * pkin(3);
t263 = qJ(2) + t171;
t103 = pkin(4) * t213 + t263;
t39 = pkin(5) * t214 - pkin(9) * t86 + t103;
t37 = t186 * t70 - t190 * t69;
t115 = -t191 * qJD(4) + t262 * t251;
t116 = -qJD(3) * t139 - t187 * qJD(4);
t66 = t183 * t115 - t116 * t182;
t49 = pkin(8) * t125 + t66;
t67 = t182 * t115 + t183 * t116;
t50 = pkin(8) * t124 + t67;
t8 = -qJD(5) * t37 + t186 * t49 + t190 * t50;
t291 = -t18 * t45 - t214 * t231 - (qJD(6) * t39 + t8) * t309 + t3 * t86 - t38 * t30;
t286 = g(3) * t187;
t285 = t18 * t86;
t284 = t39 * t30;
t278 = t185 * t60;
t28 = t189 * t30;
t72 = t183 * t120 - t98;
t207 = t161 * t190 - t186 * t289;
t71 = -t120 * t182 - t268;
t51 = t71 + t288;
t52 = t72 - t287;
t274 = -t207 * qJD(5) + t186 * t51 + t190 * t52;
t273 = t259 * qJD(5) - t186 * t52 + t190 * t51;
t272 = pkin(1) * qJDD(1);
t267 = t185 * t188;
t266 = t185 * t192;
t265 = t188 * t189;
t264 = t189 * t192;
t258 = t192 * pkin(1) + t188 * qJ(2);
t181 = t191 ^ 2;
t256 = t187 ^ 2 - t181;
t194 = qJD(3) ^ 2;
t255 = -t194 - t195;
t254 = qJD(1) * t137;
t246 = pkin(3) * t250 + qJD(2);
t243 = qJDD(3) * t187;
t235 = t86 * t248;
t234 = t191 * t245;
t211 = qJDD(4) + t177 + t178 + (t234 + t242) * pkin(3);
t53 = -pkin(4) * t82 + t211;
t5 = pkin(5) * t32 - pkin(9) * t31 + t53;
t232 = qJD(6) * t19 - t5;
t224 = t185 * t309;
t221 = qJD(6) * t214 + qJD(1);
t220 = qJDD(2) - t272;
t91 = -pkin(4) * t124 + t246;
t216 = t30 * t86 - t309 * t45;
t210 = t28 - (t185 * t75 + t248) * t309;
t208 = -t231 + t151;
t206 = 0.2e1 * qJ(2) * t245 + qJDD(3) * t193;
t205 = -pkin(9) * t30 + t22 * t309 + t313;
t202 = -t118 * t30 + t274 * t309 + t313;
t200 = -t124 * t64 - t125 * t63 - t212 * t35 + t213 * t36 - t257;
t198 = -t193 * t194 + t298;
t184 = -qJ(4) - pkin(7);
t170 = t192 * qJ(2);
t167 = qJDD(3) * t191;
t117 = -pkin(5) - t207;
t114 = t159 * t264 - t267;
t113 = t159 * t266 + t265;
t112 = t159 * t265 + t266;
t111 = -t159 * t267 + t264;
t15 = pkin(5) * t199 + pkin(9) * t45 + t91;
t9 = qJD(5) * t38 + t186 * t50 - t190 * t49;
t4 = t189 * t5;
t1 = [qJDD(1), t257, t218, qJDD(2) - t257 - 0.2e1 * t272, t298, -t220 * pkin(1) - g(1) * (-pkin(1) * t188 + t170) - g(2) * t258 + (t237 + t177) * qJ(2), qJDD(1) * t181 - 0.2e1 * t187 * t234, -0.2e1 * t187 * t241 + 0.2e1 * t256 * t245, -t187 * t194 + t167, -t191 * t194 - t243, 0, t187 * t198 + t191 * t206, -t187 * t206 + t191 * t198, -t123 * t67 - t126 * t66 + t82 * t89 + t83 * t88 - t200, t36 * t89 + t64 * t67 + t35 * t88 + t63 * t66 + t211 * t263 + t137 * t246 - g(1) * (t192 * t171 + t170 + (-pkin(1) + t184) * t188) - g(2) * (t188 * t171 - t184 * t192 + t258) -t215 * t45 + t31 * t86, -t199 * t215 - t214 * t31 - t32 * t86 + t45 * t75, t318, t292, 0, t103 * t32 - t159 * t218 - t37 * t175 - t9 * t176 + t199 * t90 + t214 * t53 + t91 * t75, t103 * t31 - t160 * t218 - t38 * t175 - t8 * t176 + t215 * t91 - t45 * t90 + t53 * t86, t319 * t189 - t60 * t235 -(-t189 * t58 - t278) * t45 + (-t11 - t14 * t189 + (t185 * t58 - t189 * t60) * qJD(6)) * t86, t13 * t214 + t189 * t216 + t199 * t60 - t235 * t309, -t14 * t214 - t185 * t216 - t199 * t58 - t65 * t86, t199 * t309 + t214 * t30, -g(1) * t114 - g(2) * t112 + t37 * t14 + t4 * t214 + t6 * t199 + t9 * t58 + (t15 * t309 + t284 + (-t19 * t214 - t309 * t38 + t285) * qJD(6)) * t189 + t291 * t185, g(1) * t113 - g(2) * t111 + t37 * t13 - t7 * t199 + t9 * t60 + (-(-qJD(6) * t38 + t15) * t309 - t284 + t232 * t214 - qJD(6) * t285) * t185 + t291 * t189; 0, 0, 0, qJDD(1), -t195, t220 + t203, 0, 0, 0, 0, 0, t255 * t187 + t167, t255 * t191 - t243, t123 * t124 + t125 * t126 - t212 * t83 + t213 * t82, t200 - t254, 0, 0, 0, 0, 0, -qJD(1) * t75 + t318, -qJD(1) * t215 + t292, 0, 0, 0, 0, 0, -t214 * t27 - t86 * t14 + t45 * t58 + (-t185 * t199 - t189 * t221) * t309, -t214 * t28 + (t185 * t221 - t189 * t199) * t309 - t319; 0, 0, 0, 0, 0, 0, t191 * t195 * t187, -t256 * t195, t241, -t242, qJDD(3), t191 * t203 + t135 + t286, g(3) * t191 + (-t144 - t203) * t187 (t64 + t71) * t126 - (t63 - t72) * t123 + (t182 * t82 + t183 * t83) * pkin(3), -t63 * t71 - t64 * t72 + (t286 + t182 * t36 + t183 * t35 + (-t257 - t254) * t191) * pkin(3), t310, t305, t320, t307, t175, t207 * t175 - t273 * t176 - t94 * t75 + t295, -t259 * t175 + t176 * t274 - t215 * t94 + t304, t316, -t278 * t309 + t317, t315, t210 + t283, -t299, t117 * t14 + t273 * t58 + (t219 - t301) * t189 + t202 * t185 + t297, t117 * t13 + t273 * t60 + t202 * t189 + (-t308 + t301) * t185 + t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123 ^ 2 - t126 ^ 2, t123 * t64 + t126 * t63 + t211 - t218, 0, 0, 0, 0, 0, t32 + t275, t31 - t279, 0, 0, 0, 0, 0, t210 - t283, -t189 * t309 ^ 2 - t27 - t282; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t310, t305, t320, t307, t175, t176 * t23 + t295, t176 * t22 + t304, t316, -t224 * t60 + t317, t315, -t224 * t309 + t28 + t283, -t299, -pkin(5) * t14 - t23 * t58 + t205 * t185 + (t219 - t300) * t189 + t297, -pkin(5) * t13 - t23 * t60 + t205 * t189 + (-t308 + t300) * t185 + t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t58, -t58 ^ 2 + t60 ^ 2, t309 * t58 + t13, t309 * t60 - t14, t30, -g(1) * t111 - g(2) * t113 - t18 * t60 + t185 * t208 - t19 * t247 + t309 * t7 + t4, g(1) * t112 - g(2) * t114 + t18 * t58 + t185 * t232 + t189 * t208 + t309 * t6;];
tau_reg  = t1;
