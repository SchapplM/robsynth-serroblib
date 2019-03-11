% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:40:43
% EndTime: 2019-03-09 03:41:02
% DurationCPUTime: 10.10s
% Computational Cost: add. (9389->586), mult. (22606->834), div. (0->0), fcn. (15859->10), ass. (0->251)
t231 = sin(pkin(11));
t233 = cos(pkin(11));
t236 = sin(qJ(5));
t239 = cos(qJ(5));
t202 = t231 * t239 + t233 * t236;
t240 = cos(qJ(3));
t247 = t202 * t240;
t165 = qJD(1) * t247;
t188 = t202 * qJD(5);
t281 = -t165 + t188;
t284 = t233 * t239;
t250 = t231 * t236 - t284;
t246 = t250 * t240;
t166 = qJD(1) * t246;
t187 = t250 * qJD(5);
t280 = -t166 + t187;
t226 = sin(pkin(10)) * pkin(1) + pkin(7);
t213 = t226 * qJD(1);
t237 = sin(qJ(3));
t203 = t237 * t213;
t180 = qJD(2) * t240 - t203;
t254 = pkin(3) * t237 - qJ(4) * t240;
t206 = t254 * qJD(1);
t128 = t233 * t180 + t231 * t206;
t278 = qJD(1) * t240;
t265 = t231 * t278;
t109 = -pkin(8) * t265 + t128;
t310 = pkin(8) + qJ(4);
t208 = t310 * t231;
t209 = t310 * t233;
t152 = -t236 * t208 + t239 * t209;
t127 = -t231 * t180 + t233 * t206;
t283 = t233 * t240;
t249 = pkin(4) * t237 - pkin(8) * t283;
t98 = qJD(1) * t249 + t127;
t358 = -t202 * qJD(4) - qJD(5) * t152 + t109 * t236 - t239 * t98;
t274 = qJD(5) * t239;
t357 = qJD(4) * t284 - t239 * t109 - t208 * t274 + (-qJD(4) * t231 - qJD(5) * t209 - t98) * t236;
t279 = qJD(1) * t237;
t364 = -pkin(5) * t279 + t280 * pkin(9) + t358;
t363 = t281 * pkin(9) - t357;
t235 = sin(qJ(6));
t238 = cos(qJ(6));
t196 = t233 * qJD(3) - t231 * t279;
t273 = t231 * qJD(3);
t197 = t233 * t279 + t273;
t259 = t239 * t196 - t197 * t236;
t181 = t237 * qJD(2) + t240 * t213;
t162 = qJD(3) * qJ(4) + t181;
t266 = -cos(pkin(10)) * pkin(1) - pkin(2);
t193 = -pkin(3) * t240 - t237 * qJ(4) + t266;
t167 = t193 * qJD(1);
t101 = -t231 * t162 + t233 * t167;
t79 = -pkin(4) * t278 - t197 * pkin(8) + t101;
t102 = t233 * t162 + t231 * t167;
t84 = pkin(8) * t196 + t102;
t42 = t236 * t79 + t239 * t84;
t36 = pkin(9) * t259 + t42;
t295 = t238 * t36;
t225 = qJD(5) - t278;
t141 = t196 * t236 + t197 * t239;
t41 = -t236 * t84 + t239 * t79;
t35 = -pkin(9) * t141 + t41;
t34 = pkin(5) * t225 + t35;
t10 = t235 * t34 + t295;
t272 = qJD(1) * qJD(3);
t262 = t237 * t272;
t355 = -t141 * t235 + t238 * t259;
t243 = qJD(3) * t246;
t92 = -qJD(1) * t243 + qJD(5) * t259;
t244 = qJD(3) * t247;
t93 = -qJD(1) * t244 - qJD(5) * t141;
t30 = qJD(6) * t355 + t235 * t93 + t238 * t92;
t71 = t141 * t238 + t235 * t259;
t31 = -qJD(6) * t71 - t235 * t92 + t238 * t93;
t271 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t262;
t317 = Ifges(7,4) * t71;
t220 = qJD(6) + t225;
t324 = -t220 / 0.2e1;
t337 = -t71 / 0.2e1;
t339 = -t355 / 0.2e1;
t245 = t249 * qJD(3);
t276 = qJD(3) * t240;
t228 = qJD(2) * t276;
t159 = t228 + (qJD(4) - t203) * qJD(3);
t185 = qJD(3) * t254 - t237 * qJD(4);
t168 = t185 * qJD(1);
t96 = -t231 * t159 + t233 * t168;
t80 = qJD(1) * t245 + t96;
t261 = t240 * t272;
t258 = t231 * t261;
t97 = t233 * t159 + t231 * t168;
t85 = -pkin(8) * t258 + t97;
t16 = -qJD(5) * t42 - t236 * t85 + t239 * t80;
t13 = pkin(5) * t262 - pkin(9) * t92 + t16;
t275 = qJD(5) * t236;
t15 = t236 * t80 + t239 * t85 + t79 * t274 - t275 * t84;
t14 = pkin(9) * t93 + t15;
t296 = t235 * t36;
t9 = t238 * t34 - t296;
t2 = qJD(6) * t9 + t13 * t235 + t14 * t238;
t3 = -qJD(6) * t10 + t13 * t238 - t14 * t235;
t349 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t64 = Ifges(7,4) * t355;
t39 = Ifges(7,1) * t71 + Ifges(7,5) * t220 + t64;
t160 = -qJD(3) * pkin(3) + qJD(4) - t180;
t129 = -t196 * pkin(4) + t160;
t77 = -pkin(5) * t259 + t129;
t362 = t271 + t349 + (Ifges(7,5) * t355 - Ifges(7,6) * t71) * t324 + (t10 * t71 + t355 * t9) * mrSges(7,3) + (-Ifges(7,2) * t71 + t39 + t64) * t339 - t77 * (mrSges(7,1) * t71 + mrSges(7,2) * t355) + (Ifges(7,1) * t355 - t317) * t337;
t361 = Ifges(4,5) / 0.2e1;
t151 = -t239 * t208 - t209 * t236;
t121 = -pkin(9) * t202 + t151;
t122 = -pkin(9) * t250 + t152;
t56 = t121 * t235 + t122 * t238;
t360 = -qJD(6) * t56 + t363 * t235 + t364 * t238;
t55 = t121 * t238 - t122 * t235;
t359 = qJD(6) * t55 + t364 * t235 - t363 * t238;
t153 = pkin(4) * t265 + t181;
t356 = t281 * pkin(5) - t153;
t38 = Ifges(7,2) * t355 + Ifges(7,6) * t220 + t317;
t353 = t38 / 0.2e1;
t352 = t279 / 0.2e1;
t263 = -Ifges(4,6) * qJD(3) / 0.2e1;
t351 = qJD(3) * t361;
t179 = t233 * t193;
t120 = -t233 * t237 * pkin(8) + t179 + (-t226 * t231 - pkin(4)) * t240;
t145 = t231 * t193 + t226 * t283;
t286 = t231 * t237;
t126 = -pkin(8) * t286 + t145;
t58 = t236 * t120 + t239 * t126;
t344 = -t16 * mrSges(6,1) + t15 * mrSges(6,2);
t277 = qJD(3) * t237;
t309 = mrSges(5,2) * t233;
t172 = mrSges(5,1) * t258 + t261 * t309;
t49 = -t93 * mrSges(6,1) + t92 * mrSges(6,2);
t8 = -t31 * mrSges(7,1) + t30 * mrSges(7,2);
t343 = -t172 - t49 - t8;
t215 = t266 * qJD(1);
t229 = Ifges(4,4) * t278;
t303 = Ifges(5,2) * t231;
t307 = Ifges(5,4) * t233;
t255 = -t303 + t307;
t308 = Ifges(5,4) * t231;
t256 = Ifges(5,1) * t233 - t308;
t257 = mrSges(5,1) * t231 + t309;
t319 = t233 / 0.2e1;
t320 = -t231 / 0.2e1;
t342 = -(t101 * t233 + t102 * t231) * mrSges(5,3) + t160 * t257 + t215 * mrSges(4,2) + Ifges(4,1) * t352 + t229 / 0.2e1 + t351 - t180 * mrSges(4,3) + t196 * t255 / 0.2e1 + t197 * t256 / 0.2e1 + (Ifges(5,4) * t197 + Ifges(5,2) * t196 - Ifges(5,6) * t278) * t320 + (Ifges(5,1) * t197 + Ifges(5,4) * t196 - Ifges(5,5) * t278) * t319;
t341 = t30 / 0.2e1;
t340 = t31 / 0.2e1;
t338 = t355 / 0.2e1;
t336 = t71 / 0.2e1;
t335 = t92 / 0.2e1;
t334 = t93 / 0.2e1;
t176 = t202 * t237;
t177 = t250 * t237;
t117 = -t176 * t238 + t177 * t235;
t332 = t117 / 0.2e1;
t118 = -t176 * t235 - t177 * t238;
t331 = t118 / 0.2e1;
t330 = -t259 / 0.2e1;
t329 = t259 / 0.2e1;
t328 = -t141 / 0.2e1;
t327 = t141 / 0.2e1;
t326 = -t176 / 0.2e1;
t325 = -t177 / 0.2e1;
t323 = t220 / 0.2e1;
t322 = -t225 / 0.2e1;
t321 = t225 / 0.2e1;
t306 = Ifges(6,4) * t141;
t304 = Ifges(5,5) * t233;
t301 = Ifges(5,6) * t231;
t104 = -t165 * t238 + t166 * t235;
t143 = t202 * t238 - t235 * t250;
t73 = -qJD(6) * t143 + t187 * t235 - t188 * t238;
t294 = t104 - t73;
t105 = -t165 * t235 - t166 * t238;
t142 = -t202 * t235 - t238 * t250;
t72 = qJD(6) * t142 - t187 * t238 - t188 * t235;
t293 = t105 - t72;
t169 = -t213 * t277 + t228;
t289 = t169 * t240;
t170 = t181 * qJD(3);
t288 = t170 * t240;
t285 = t231 * t240;
t217 = t237 * t226;
t269 = mrSges(4,3) * t279;
t282 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t196 + mrSges(5,2) * t197 + t269;
t264 = t226 * t277;
t136 = t233 * t185 + t231 * t264;
t175 = t240 * pkin(4) * t273 + t226 * t276;
t184 = pkin(4) * t286 + t217;
t270 = Ifges(6,5) * t92 + Ifges(6,6) * t93 + Ifges(6,3) * t262;
t268 = mrSges(4,3) * t278;
t267 = t170 * t217;
t227 = -pkin(4) * t233 - pkin(3);
t57 = t239 * t120 - t236 * t126;
t253 = -t231 * t96 + t233 * t97;
t50 = -pkin(5) * t240 + t177 * pkin(9) + t57;
t51 = -pkin(9) * t176 + t58;
t21 = -t235 * t51 + t238 * t50;
t22 = t235 * t50 + t238 * t51;
t251 = -t101 * t231 + t102 * t233;
t107 = t245 + t136;
t173 = t231 * t185;
t119 = t173 + (-pkin(8) * t285 - t217 * t233) * qJD(3);
t32 = t236 * t107 + t239 * t119 + t120 * t274 - t126 * t275;
t146 = pkin(4) * t258 + t170;
t33 = -qJD(5) * t58 + t239 * t107 - t119 * t236;
t241 = t101 * mrSges(5,1) + t215 * mrSges(4,1) + t41 * mrSges(6,1) + t9 * mrSges(7,1) - Ifges(5,3) * t278 / 0.2e1 + Ifges(5,6) * t196 + Ifges(5,5) * t197 + t263 - (Ifges(4,4) * t237 + Ifges(4,2) * t240) * qJD(1) / 0.2e1 + t220 * Ifges(7,3) + t71 * Ifges(7,5) + t355 * Ifges(7,6) + t225 * Ifges(6,3) + t141 * Ifges(6,5) + t259 * Ifges(6,6) - t10 * mrSges(7,2) - t102 * mrSges(5,2) - t42 * mrSges(6,2);
t216 = -qJD(3) * mrSges(4,2) + t268;
t183 = (mrSges(5,1) * t237 - mrSges(5,3) * t283) * t272;
t182 = (-mrSges(5,2) * t237 - mrSges(5,3) * t285) * t272;
t171 = pkin(5) * t250 + t227;
t164 = -mrSges(5,1) * t278 - t197 * mrSges(5,3);
t163 = mrSges(5,2) * t278 + t196 * mrSges(5,3);
t150 = (Ifges(5,5) * t237 + t240 * t256) * t272;
t149 = (Ifges(5,6) * t237 + t240 * t255) * t272;
t144 = -t226 * t285 + t179;
t137 = -t233 * t264 + t173;
t135 = Ifges(6,4) * t259;
t131 = pkin(5) * t176 + t184;
t125 = t187 * t237 - t244;
t124 = -t188 * t237 - t243;
t116 = mrSges(6,1) * t225 - mrSges(6,3) * t141;
t115 = -mrSges(6,2) * t225 + mrSges(6,3) * t259;
t86 = -pkin(5) * t125 + t175;
t82 = -mrSges(6,2) * t262 + mrSges(6,3) * t93;
t81 = mrSges(6,1) * t262 - mrSges(6,3) * t92;
t78 = -mrSges(6,1) * t259 + mrSges(6,2) * t141;
t63 = Ifges(6,1) * t141 + Ifges(6,5) * t225 + t135;
t62 = Ifges(6,2) * t259 + Ifges(6,6) * t225 + t306;
t60 = mrSges(7,1) * t220 - mrSges(7,3) * t71;
t59 = -mrSges(7,2) * t220 + mrSges(7,3) * t355;
t54 = -t93 * pkin(5) + t146;
t47 = t92 * Ifges(6,1) + t93 * Ifges(6,4) + Ifges(6,5) * t262;
t46 = t92 * Ifges(6,4) + t93 * Ifges(6,2) + Ifges(6,6) * t262;
t44 = -qJD(6) * t118 - t124 * t235 + t125 * t238;
t43 = qJD(6) * t117 + t124 * t238 + t125 * t235;
t40 = -mrSges(7,1) * t355 + mrSges(7,2) * t71;
t26 = -mrSges(7,2) * t262 + mrSges(7,3) * t31;
t25 = mrSges(7,1) * t262 - mrSges(7,3) * t30;
t20 = pkin(9) * t125 + t32;
t19 = pkin(5) * t277 - pkin(9) * t124 + t33;
t12 = t238 * t35 - t296;
t11 = -t235 * t35 - t295;
t7 = t30 * Ifges(7,1) + t31 * Ifges(7,4) + Ifges(7,5) * t262;
t6 = t30 * Ifges(7,4) + t31 * Ifges(7,2) + Ifges(7,6) * t262;
t5 = -qJD(6) * t22 + t19 * t238 - t20 * t235;
t4 = qJD(6) * t21 + t19 * t235 + t20 * t238;
t1 = [(t226 * t172 + t149 * t320 + t150 * t319 + (mrSges(4,3) + t257) * t170 + (-t231 * t97 - t233 * t96) * mrSges(5,3) + ((-m(4) * t226 - mrSges(4,3)) * t181 + t241 + t263 - t226 * t216) * qJD(3) + (Ifges(6,5) * t325 + Ifges(6,6) * t326 + Ifges(7,5) * t331 + Ifges(7,6) * t332 + t266 * mrSges(4,1) + (t304 / 0.2e1 - t301 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,4)) * t237) * t272) * t237 + m(4) * (t226 * t289 + t267) + (Ifges(7,4) * t118 + Ifges(7,2) * t117) * t340 + m(5) * (t101 * t136 + t102 * t137 + t96 * t144 + t97 * t145 + t267) + t44 * t353 + ((t351 + (-m(4) * t180 + m(5) * t160 + t282) * t226 + t342) * qJD(3) + t97 * mrSges(5,2) - t96 * mrSges(5,1) - Ifges(6,6) * t334 - Ifges(6,5) * t335 - Ifges(7,6) * t340 - Ifges(7,5) * t341 + (t266 * mrSges(4,2) + (0.3e1 / 0.2e1 * Ifges(4,4) + 0.3e1 / 0.2e1 * t301 - 0.3e1 / 0.2e1 * t304) * t240 + (-0.3e1 / 0.2e1 * Ifges(5,3) - Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1) + Ifges(5,1) * t233 ^ 2 / 0.2e1 + (-t307 + t303 / 0.2e1) * t231) * t237) * t272 + t344 - t349) * t240 + m(7) * (t10 * t4 + t131 * t54 + t2 * t22 + t21 * t3 + t5 * t9 + t77 * t86) + m(6) * (t129 * t175 + t146 * t184 + t15 * t58 + t16 * t57 + t32 * t42 + t33 * t41) - (t271 + t270) * t240 / 0.2e1 + (-Ifges(6,4) * t177 - Ifges(6,2) * t176) * t334 + (-Ifges(6,1) * t177 - Ifges(6,4) * t176) * t335 + (-t124 * t41 + t125 * t42 - t15 * t176 + t16 * t177) * mrSges(6,3) + t146 * (mrSges(6,1) * t176 - mrSges(6,2) * t177) + mrSges(4,3) * t289 + (Ifges(6,5) * t124 + Ifges(6,6) * t125) * t321 + (Ifges(7,5) * t43 + Ifges(7,6) * t44) * t323 + t21 * t25 + t22 * t26 + t43 * t39 / 0.2e1 + t4 * t59 + t5 * t60 + t77 * (-mrSges(7,1) * t44 + mrSges(7,2) * t43) + t57 * t81 + (Ifges(7,1) * t118 + Ifges(7,4) * t117) * t341 + t58 * t82 + t86 * t40 + t32 * t115 + t33 * t116 + t54 * (-mrSges(7,1) * t117 + mrSges(7,2) * t118) + t124 * t63 / 0.2e1 + t125 * t62 / 0.2e1 + t129 * (-mrSges(6,1) * t125 + mrSges(6,2) * t124) + t131 * t8 + t137 * t163 + t136 * t164 + t175 * t78 + t145 * t182 + t144 * t183 + t184 * t49 + t47 * t325 + t46 * t326 + (Ifges(6,1) * t124 + Ifges(6,4) * t125) * t327 + (Ifges(6,4) * t124 + Ifges(6,2) * t125) * t329 + t7 * t331 + t6 * t332 + (Ifges(7,1) * t43 + Ifges(7,4) * t44) * t336 + (Ifges(7,4) * t43 + Ifges(7,2) * t44) * t338 + (t10 * t44 + t117 * t2 - t118 * t3 - t43 * t9) * mrSges(7,3); t124 * t115 + t125 * t116 + t117 * t25 + t118 * t26 - t176 * t81 - t177 * t82 + t43 * t59 + t44 * t60 + (t182 * t233 - t183 * t231) * t237 + t343 * t240 + ((t163 * t233 - t164 * t231 + t216 - t268) * t240 + (t40 + t78 - t269 + t282) * t237) * qJD(3) + m(7) * (t10 * t43 + t3 * t117 + t2 * t118 - t240 * t54 + t277 * t77 + t9 * t44) + m(6) * (t42 * t124 + t41 * t125 + t129 * t277 - t146 * t240 - t15 * t177 - t16 * t176) + m(4) * (t169 * t237 - t288 + (-t180 * t237 + t181 * t240) * qJD(3)) + m(5) * (-t288 + t253 * t237 + (t160 * t237 + t240 * t251) * qJD(3)); (mrSges(7,1) * t294 - mrSges(7,2) * t293) * t77 + (-t10 * t294 + t142 * t2 - t143 * t3 + t293 * t9) * mrSges(7,3) - t250 * t46 / 0.2e1 + (t165 / 0.2e1 - t188 / 0.2e1) * t62 + (-Ifges(6,5) * t166 - Ifges(6,6) * t165) * t322 + (-Ifges(6,1) * t166 - Ifges(6,4) * t165) * t328 + (-Ifges(6,4) * t166 - Ifges(6,2) * t165) * t330 + (t166 / 0.2e1 - t187 / 0.2e1) * t63 + (t97 * mrSges(5,3) + qJD(4) * t163 + qJ(4) * t182 + t149 / 0.2e1 - t170 * mrSges(5,1)) * t233 + (mrSges(6,1) * t281 - mrSges(6,2) * t280) * t129 - t282 * t181 + (-t96 * mrSges(5,3) - qJD(4) * t164 - qJ(4) * t183 + t150 / 0.2e1 + t170 * mrSges(5,2)) * t231 + (-t15 * t250 - t16 * t202 + t280 * t41 - t281 * t42) * mrSges(6,3) + t146 * (mrSges(6,1) * t250 + mrSges(6,2) * t202) + (Ifges(6,4) * t202 - Ifges(6,2) * t250) * t334 + (Ifges(6,1) * t202 - Ifges(6,4) * t250) * t335 + t359 * t59 + (t10 * t359 + t171 * t54 + t2 * t56 + t3 * t55 + t356 * t77 + t360 * t9) * m(7) + t360 * t60 + ((-t229 / 0.2e1 + (-t301 + t304) * t278 / 0.2e1 + (t361 + (Ifges(5,1) * t231 + t307) * t319 + (Ifges(5,2) * t233 + t308) * t320) * qJD(3) - t342) * t240 + ((-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t278 - t241 + t263 + t181 * mrSges(4,3) + Ifges(4,4) * t352) * t237 + (Ifges(5,5) * t231 + Ifges(6,5) * t202 + Ifges(7,5) * t143 + Ifges(5,6) * t233 - Ifges(6,6) * t250 + Ifges(7,6) * t142) * t277 / 0.2e1) * qJD(1) + t357 * t115 + (-t129 * t153 + t146 * t227 + t15 * t152 + t151 * t16 + t357 * t42 + t358 * t41) * m(6) + t358 * t116 + (-pkin(3) * t170 + qJ(4) * t253 + qJD(4) * t251 - t101 * t127 - t102 * t128 - t160 * t181) * m(5) + t356 * t40 + (t73 / 0.2e1 - t104 / 0.2e1) * t38 + (t72 / 0.2e1 - t105 / 0.2e1) * t39 + (-Ifges(6,5) * t187 - Ifges(6,6) * t188) * t321 + (-Ifges(6,1) * t187 - Ifges(6,4) * t188) * t327 + (-Ifges(6,4) * t187 - Ifges(6,2) * t188) * t329 + (Ifges(7,5) * t72 + Ifges(7,6) * t73) * t323 + (Ifges(7,5) * t105 + Ifges(7,6) * t104) * t324 + t55 * t25 + t56 * t26 + t142 * t6 / 0.2e1 + t143 * t7 / 0.2e1 + t54 * (-mrSges(7,1) * t142 + mrSges(7,2) * t143) + t151 * t81 + t152 * t82 - t153 * t78 - t128 * t163 - t127 * t164 - t169 * mrSges(4,2) - t170 * mrSges(4,1) + t171 * t8 - pkin(3) * t172 + t202 * t47 / 0.2e1 - t180 * t216 + (Ifges(7,1) * t72 + Ifges(7,4) * t73) * t336 + (Ifges(7,1) * t105 + Ifges(7,4) * t104) * t337 + (Ifges(7,4) * t72 + Ifges(7,2) * t73) * t338 + (Ifges(7,4) * t105 + Ifges(7,2) * t104) * t339 + (Ifges(7,4) * t143 + Ifges(7,2) * t142) * t340 + (Ifges(7,1) * t143 + Ifges(7,4) * t142) * t341 + t227 * t49; -t259 * t115 + t141 * t116 - t196 * t163 + t197 * t164 - t355 * t59 + t71 * t60 + (-t10 * t355 + t71 * t9 + t54) * m(7) + (t141 * t41 - t259 * t42 + t146) * m(6) + (t101 * t197 - t102 * t196 + t170) * m(5) - t343; -m(7) * (t10 * t12 + t11 * t9) + t270 + t71 * t353 - t344 + (-Ifges(6,2) * t141 + t135 + t63) * t330 + (t141 * t42 + t259 * t41) * mrSges(6,3) - t129 * (mrSges(6,1) * t141 + mrSges(6,2) * t259) + (Ifges(6,1) * t259 - t306) * t328 + (Ifges(6,5) * t259 - Ifges(6,6) * t141) * t322 + (-t141 * t40 + t235 * t26 + t238 * t25 + (-t235 * t60 + t238 * t59) * qJD(6) + (-t141 * t77 + t2 * t235 + t238 * t3 + (t10 * t238 - t235 * t9) * qJD(6)) * m(7)) * pkin(5) - t12 * t59 - t11 * t60 - t41 * t115 + t42 * t116 + t62 * t327 + t362; t10 * t60 + t38 * t336 - t9 * t59 + t362;];
tauc  = t1(:);
