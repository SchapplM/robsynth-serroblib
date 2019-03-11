% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:07:39
% EndTime: 2019-03-09 04:07:59
% DurationCPUTime: 10.18s
% Computational Cost: add. (9252->590), mult. (20945->846), div. (0->0), fcn. (14494->8), ass. (0->257)
t230 = sin(pkin(10));
t231 = cos(pkin(10));
t233 = sin(qJ(5));
t236 = cos(qJ(5));
t208 = t230 * t236 + t231 * t233;
t188 = t208 * qJD(1);
t234 = sin(qJ(3));
t167 = t234 * t188;
t190 = t208 * qJD(5);
t279 = t167 + t190;
t282 = t231 * t236;
t250 = t230 * t233 - t282;
t350 = t250 * t234;
t168 = qJD(1) * t350;
t345 = qJD(5) * t250;
t278 = t168 + t345;
t237 = cos(qJ(3));
t254 = pkin(3) * t237 + qJ(4) * t234;
t210 = t254 * qJD(1);
t238 = -pkin(1) - pkin(7);
t221 = qJD(1) * t238 + qJD(2);
t284 = t230 * t237;
t151 = t231 * t210 - t221 * t284;
t283 = t231 * t234;
t270 = pkin(8) * t283;
t244 = (pkin(4) * t237 + t270) * qJD(1);
t109 = t151 + t244;
t281 = t231 * t237;
t152 = t230 * t210 + t221 * t281;
t229 = t234 * qJD(1);
t267 = t230 * t229;
t125 = pkin(8) * t267 + t152;
t310 = pkin(8) + qJ(4);
t214 = t310 * t230;
t215 = t310 * t231;
t158 = -t233 * t214 + t236 * t215;
t368 = -t208 * qJD(4) - qJD(5) * t158 - t236 * t109 + t125 * t233;
t272 = qJD(5) * t236;
t367 = qJD(4) * t282 - t236 * t125 - t214 * t272 + (-qJD(4) * t230 - qJD(5) * t215 - t109) * t233;
t276 = qJD(1) * t237;
t374 = -pkin(5) * t276 + pkin(9) * t278 + t368;
t373 = pkin(9) * t279 - t367;
t179 = t250 * t237;
t349 = -qJD(3) * t179 - t190 * t234 - t188;
t177 = t208 * t237;
t348 = t250 * qJD(1) - qJD(3) * t177 + t234 * t345;
t232 = sin(qJ(6));
t235 = cos(qJ(6));
t202 = qJD(3) * t231 - t230 * t276;
t203 = qJD(3) * t230 + t231 * t276;
t260 = t236 * t202 - t203 * t233;
t213 = pkin(3) * t234 - qJ(4) * t237 + qJ(2);
t193 = t213 * qJD(1);
t212 = t234 * t221;
t196 = qJD(3) * qJ(4) + t212;
t128 = t231 * t193 - t196 * t230;
t93 = pkin(4) * t229 - pkin(8) * t203 + t128;
t129 = t230 * t193 + t231 * t196;
t95 = pkin(8) * t202 + t129;
t49 = t233 * t93 + t236 * t95;
t39 = pkin(9) * t260 + t49;
t291 = t235 * t39;
t227 = t229 + qJD(5);
t141 = t202 * t233 + t203 * t236;
t48 = -t233 * t95 + t236 * t93;
t38 = -pkin(9) * t141 + t48;
t37 = pkin(5) * t227 + t38;
t10 = t232 * t37 + t291;
t271 = qJD(1) * qJD(3);
t263 = t237 * t271;
t365 = -t141 * t232 + t235 * t260;
t242 = qJD(3) * t350;
t90 = qJD(1) * t242 + qJD(5) * t260;
t275 = qJD(3) * t234;
t243 = t208 * t275;
t91 = qJD(1) * t243 - qJD(5) * t141;
t28 = qJD(6) * t365 + t232 * t91 + t235 * t90;
t70 = t141 * t235 + t232 * t260;
t29 = -qJD(6) * t70 - t232 * t90 + t235 * t91;
t269 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t263;
t317 = Ifges(7,4) * t70;
t220 = qJD(6) + t227;
t324 = -t220 / 0.2e1;
t339 = -t70 / 0.2e1;
t341 = -t365 / 0.2e1;
t180 = qJD(3) * t254 - qJD(4) * t237 + qJD(2);
t159 = t180 * qJD(1);
t285 = t221 * t237;
t183 = (qJD(4) + t285) * qJD(3);
t100 = t231 * t159 - t183 * t230;
t81 = qJD(3) * t244 + t100;
t101 = t230 * t159 + t231 * t183;
t266 = t230 * t275;
t259 = qJD(1) * t266;
t94 = pkin(8) * t259 + t101;
t16 = -qJD(5) * t49 - t233 * t94 + t236 * t81;
t11 = pkin(5) * t263 - pkin(9) * t90 + t16;
t273 = qJD(5) * t233;
t15 = t233 * t81 + t236 * t94 + t93 * t272 - t273 * t95;
t14 = pkin(9) * t91 + t15;
t292 = t232 * t39;
t9 = t235 * t37 - t292;
t2 = qJD(6) * t9 + t11 * t232 + t14 * t235;
t3 = -qJD(6) * t10 + t11 * t235 - t14 * t232;
t358 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t64 = Ifges(7,4) * t365;
t36 = Ifges(7,1) * t70 + Ifges(7,5) * t220 + t64;
t185 = -qJD(3) * pkin(3) + qJD(4) - t285;
t148 = -pkin(4) * t202 + t185;
t80 = -pkin(5) * t260 + t148;
t372 = t269 + t358 + (Ifges(7,5) * t365 - Ifges(7,6) * t70) * t324 + (t10 * t70 + t365 * t9) * mrSges(7,3) + (-Ifges(7,2) * t70 + t36 + t64) * t341 - t80 * (mrSges(7,1) * t70 + mrSges(7,2) * t365) + (Ifges(7,1) * t365 - t317) * t339;
t371 = -qJD(3) / 0.2e1;
t157 = -t236 * t214 - t215 * t233;
t117 = -pkin(9) * t208 + t157;
t118 = -pkin(9) * t250 + t158;
t55 = t117 * t232 + t118 * t235;
t370 = -qJD(6) * t55 + t373 * t232 + t235 * t374;
t54 = t117 * t235 - t118 * t232;
t369 = qJD(6) * t54 + t232 * t374 - t373 * t235;
t171 = -pkin(4) * t267 + t212;
t366 = pkin(5) * t279 - t171;
t35 = Ifges(7,2) * t365 + Ifges(7,6) * t220 + t317;
t363 = t35 / 0.2e1;
t264 = Ifges(4,6) * t371;
t362 = Ifges(4,5) * t371;
t47 = -t91 * mrSges(6,1) + t90 * mrSges(6,2);
t8 = -t29 * mrSges(7,1) + t28 * mrSges(7,2);
t361 = -t47 - t8;
t359 = -qJD(1) / 0.2e1;
t176 = t208 * t234;
t111 = -t176 * t235 + t232 * t350;
t354 = qJD(6) * t111 + t348 * t232 + t235 * t349;
t113 = -t176 * t232 - t235 * t350;
t353 = -qJD(6) * t113 - t232 * t349 + t348 * t235;
t352 = qJ(2) * (m(3) + m(4));
t200 = t231 * t213;
t262 = -t230 * t238 + pkin(4);
t137 = -pkin(8) * t281 + t234 * t262 + t200;
t280 = t234 * t238;
t161 = t230 * t213 + t231 * t280;
t150 = -pkin(8) * t284 + t161;
t74 = t233 * t137 + t236 * t150;
t347 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t346 = qJD(3) * t237;
t252 = t128 * t231 + t129 * t230;
t302 = Ifges(5,2) * t230;
t305 = Ifges(5,4) * t231;
t255 = t302 - t305;
t306 = Ifges(5,4) * t230;
t256 = -Ifges(5,1) * t231 + t306;
t257 = -mrSges(5,1) * t230 - mrSges(5,2) * t231;
t319 = -t231 / 0.2e1;
t320 = t230 / 0.2e1;
t344 = -t252 * mrSges(5,3) - t185 * t257 - t362 - (Ifges(4,1) * t237 - Ifges(4,4) * t234) * t359 - t202 * t255 / 0.2e1 - t203 * t256 / 0.2e1 - (t203 * Ifges(5,4) + t202 * Ifges(5,2) + Ifges(5,6) * t229) * t320 - (t203 * Ifges(5,1) + t202 * Ifges(5,4) + Ifges(5,5) * t229) * t319;
t343 = t28 / 0.2e1;
t342 = t29 / 0.2e1;
t340 = t365 / 0.2e1;
t338 = t70 / 0.2e1;
t337 = t90 / 0.2e1;
t336 = t91 / 0.2e1;
t112 = -t177 * t235 + t179 * t232;
t333 = t112 / 0.2e1;
t114 = -t177 * t232 - t179 * t235;
t332 = t114 / 0.2e1;
t331 = -t260 / 0.2e1;
t330 = t260 / 0.2e1;
t329 = -t141 / 0.2e1;
t328 = t141 / 0.2e1;
t327 = (Ifges(5,5) * t237 + t234 * t256) * t271 / 0.2e1;
t326 = -t177 / 0.2e1;
t325 = -t179 / 0.2e1;
t323 = t220 / 0.2e1;
t322 = -t227 / 0.2e1;
t321 = t227 / 0.2e1;
t143 = -t208 * t232 - t235 * t250;
t71 = qJD(6) * t143 - t190 * t232 - t235 * t345;
t99 = t167 * t232 + t168 * t235;
t309 = t71 - t99;
t144 = t208 * t235 - t232 * t250;
t72 = -qJD(6) * t144 - t190 * t235 + t232 * t345;
t98 = t167 * t235 - t168 * t232;
t308 = t72 - t98;
t307 = Ifges(4,4) * t237;
t304 = Ifges(6,4) * t141;
t303 = Ifges(5,5) * t231;
t301 = Ifges(5,6) * t230;
t300 = qJ(2) * mrSges(4,1);
t299 = qJ(2) * mrSges(4,2);
t288 = qJD(3) * mrSges(4,2);
t274 = qJD(3) * t238;
t265 = t237 * t274;
t147 = t230 * t180 + t231 * t265;
t277 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t202 - mrSges(5,2) * t203 - mrSges(4,3) * t276;
t268 = Ifges(6,5) * t90 + Ifges(6,6) * t91 + Ifges(6,3) * t263;
t228 = -pkin(4) * t231 - pkin(3);
t209 = t221 * t275;
t73 = t236 * t137 - t150 * t233;
t201 = pkin(4) * t284 - t237 * t238;
t258 = mrSges(4,1) * t234 + mrSges(4,2) * t237;
t52 = pkin(5) * t234 + pkin(9) * t179 + t73;
t53 = -pkin(9) * t177 + t74;
t30 = -t232 * t53 + t235 * t52;
t31 = t232 * t52 + t235 * t53;
t253 = -t100 * t230 + t101 * t231;
t251 = -t128 * t230 + t129 * t231;
t186 = -pkin(4) * t266 + t234 * t274;
t248 = t303 / 0.2e1 - t301 / 0.2e1;
t247 = t257 * qJD(1);
t162 = -pkin(4) * t259 + t209;
t164 = t231 * t180;
t104 = t164 + (t237 * t262 + t270) * qJD(3);
t123 = pkin(8) * t266 + t147;
t32 = t233 * t104 + t236 * t123 + t137 * t272 - t150 * t273;
t33 = -qJD(5) * t74 + t236 * t104 - t123 * t233;
t240 = t128 * mrSges(5,1) + t48 * mrSges(6,1) + t9 * mrSges(7,1) + Ifges(5,3) * t229 / 0.2e1 + t203 * Ifges(5,5) + t202 * Ifges(5,6) + t264 + (-t234 * Ifges(4,2) + t307) * t359 + t220 * Ifges(7,3) + t70 * Ifges(7,5) + t365 * Ifges(7,6) + t227 * Ifges(6,3) + t141 * Ifges(6,5) + t260 * Ifges(6,6) - t10 * mrSges(7,2) - t129 * mrSges(5,2) - t49 * mrSges(6,2);
t217 = -mrSges(4,3) * t229 - t288;
t211 = t258 * qJD(1);
t182 = (mrSges(5,1) * t237 + mrSges(5,3) * t283) * t271;
t181 = (mrSges(5,3) * t230 * t234 - mrSges(5,2) * t237) * t271;
t170 = t247 * t275;
t169 = pkin(5) * t250 + t228;
t166 = mrSges(5,1) * t229 - mrSges(5,3) * t203;
t165 = -mrSges(5,2) * t229 + mrSges(5,3) * t202;
t160 = -t230 * t280 + t200;
t155 = (Ifges(5,6) * t237 + t234 * t255) * t271;
t146 = -t230 * t265 + t164;
t145 = pkin(5) * t177 + t201;
t136 = Ifges(6,4) * t260;
t122 = t237 * t345 + t243;
t120 = -t190 * t237 + t242;
t108 = mrSges(6,1) * t227 - mrSges(6,3) * t141;
t107 = -mrSges(6,2) * t227 + mrSges(6,3) * t260;
t92 = -pkin(5) * t122 + t186;
t77 = -mrSges(6,2) * t263 + mrSges(6,3) * t91;
t76 = mrSges(6,1) * t263 - mrSges(6,3) * t90;
t75 = -mrSges(6,1) * t260 + mrSges(6,2) * t141;
t63 = Ifges(6,1) * t141 + Ifges(6,5) * t227 + t136;
t62 = Ifges(6,2) * t260 + Ifges(6,6) * t227 + t304;
t60 = -pkin(5) * t91 + t162;
t59 = mrSges(7,1) * t220 - mrSges(7,3) * t70;
t58 = -mrSges(7,2) * t220 + mrSges(7,3) * t365;
t46 = t90 * Ifges(6,1) + t91 * Ifges(6,4) + Ifges(6,5) * t263;
t45 = t90 * Ifges(6,4) + t91 * Ifges(6,2) + Ifges(6,6) * t263;
t44 = -qJD(6) * t114 - t120 * t232 + t122 * t235;
t42 = qJD(6) * t112 + t120 * t235 + t122 * t232;
t40 = -mrSges(7,1) * t365 + mrSges(7,2) * t70;
t24 = pkin(9) * t122 + t32;
t23 = -mrSges(7,2) * t263 + mrSges(7,3) * t29;
t22 = mrSges(7,1) * t263 - mrSges(7,3) * t28;
t21 = pkin(5) * t346 - pkin(9) * t120 + t33;
t13 = t235 * t38 - t292;
t12 = -t232 * t38 - t291;
t7 = t28 * Ifges(7,1) + t29 * Ifges(7,4) + Ifges(7,5) * t263;
t6 = t28 * Ifges(7,4) + t29 * Ifges(7,2) + Ifges(7,6) * t263;
t5 = -qJD(6) * t31 + t21 * t235 - t232 * t24;
t4 = qJD(6) * t30 + t21 * t232 + t235 * t24;
t1 = [(-t120 * t48 + t122 * t49 - t15 * t177 + t16 * t179) * mrSges(6,3) + (-Ifges(6,4) * t179 - Ifges(6,2) * t177) * t336 + (-Ifges(6,1) * t179 - Ifges(6,4) * t177) * t337 + t162 * (mrSges(6,1) * t177 - mrSges(6,2) * t179) + (t100 * mrSges(5,1) - t101 * mrSges(5,2) + Ifges(6,5) * t337 + Ifges(7,5) * t343 + Ifges(6,6) * t336 + Ifges(7,6) * t342 + (t362 + (m(5) * t185 - t277) * t238 + (-0.2e1 * t299 + (0.3e1 / 0.2e1 * Ifges(4,4) - 0.3e1 / 0.2e1 * t303 + 0.3e1 / 0.2e1 * t301) * t234 + (0.3e1 / 0.2e1 * Ifges(5,3) + Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2) - Ifges(5,1) * t231 ^ 2 / 0.2e1 + (t305 - t302 / 0.2e1) * t230) * t237) * qJD(1) - t344) * qJD(3) + t347 + t358) * t234 + (t231 * t327 - t230 * t155 / 0.2e1 - t238 * t170 + (-t100 * t231 - t101 * t230) * mrSges(5,3)) * t237 + m(5) * (t100 * t160 + t101 * t161 + t128 * t146 + t129 * t147) + (t10 * t44 + t112 * t2 - t114 * t3 - t42 * t9) * mrSges(7,3) + (t269 + t268) * t234 / 0.2e1 + (((2 * mrSges(3,3)) + t258 + 0.2e1 * t352) * qJD(2) + (Ifges(6,5) * t325 + Ifges(6,6) * t326 + Ifges(7,5) * t332 + Ifges(7,6) * t333 + 0.2e1 * t300 + (-0.3e1 / 0.2e1 * Ifges(4,4) + t248) * t237) * t346) * qJD(1) + (Ifges(7,4) * t114 + Ifges(7,2) * t112) * t342 + (Ifges(7,1) * t114 + Ifges(7,4) * t112) * t343 + t30 * t22 + t31 * t23 + t42 * t36 / 0.2e1 + t4 * t58 + t5 * t59 + (Ifges(6,5) * t120 + Ifges(6,6) * t122) * t321 + (Ifges(7,5) * t42 + Ifges(7,6) * t44) * t323 + t46 * t325 + t45 * t326 + (Ifges(6,1) * t120 + Ifges(6,4) * t122) * t328 + (Ifges(6,4) * t120 + Ifges(6,2) * t122) * t330 + t7 * t332 + t6 * t333 + (Ifges(7,1) * t42 + Ifges(7,4) * t44) * t338 + (Ifges(7,4) * t42 + Ifges(7,2) * t44) * t340 + t73 * t76 + t74 * t77 + t80 * (-mrSges(7,1) * t44 + mrSges(7,2) * t42) + t92 * t40 + t32 * t107 + t33 * t108 + t60 * (-mrSges(7,1) * t112 + mrSges(7,2) * t114) + t120 * t63 / 0.2e1 + t122 * t62 / 0.2e1 + t145 * t8 + t148 * (-mrSges(6,1) * t122 + mrSges(6,2) * t120) + t147 * t165 + t146 * t166 + t161 * t181 + t160 * t182 + t186 * t75 + t201 * t47 + m(6) * (t148 * t186 + t15 * t74 + t16 * t73 + t162 * t201 + t32 * t49 + t33 * t48) + m(7) * (t10 * t4 + t145 * t60 + t2 * t31 + t3 * t30 + t5 * t9 + t80 * t92) + qJD(2) * t211 + ((-m(5) * t238 - t257) * t212 + t238 * t217 + t240 + t264) * t346 + t44 * t363; t111 * t22 + t113 * t23 - t176 * t76 - t350 * t77 + t353 * t59 + t354 * t58 + (t181 * t231 - t182 * t230) * t234 + t348 * t108 + t349 * t107 + (-t170 + t361) * t237 + ((t231 * t165 - t230 * t166 + t217) * t237 + (t40 + t75 - t277) * t234) * qJD(3) + m(5) * (t253 * t234 + (t185 * t234 + (t251 - t212) * t237) * qJD(3)) + (t10 * t354 + t111 * t3 + t113 * t2 - t237 * t60 + t80 * t275 + t353 * t9) * m(7) + (t148 * t275 - t15 * t350 - t16 * t176 - t162 * t237 + t348 * t48 + t349 * t49) * m(6) + (-m(5) * t252 - t165 * t230 - t166 * t231 - t211 + (-mrSges(3,3) - t352) * qJD(1)) * qJD(1); (t71 / 0.2e1 - t99 / 0.2e1) * t36 + (mrSges(6,1) * t279 - mrSges(6,2) * t278) * t148 + (-t100 * mrSges(5,3) - qJ(4) * t182 - qJD(4) * t166 + t327) * t230 + (-mrSges(7,1) * t308 + mrSges(7,2) * t309) * t80 + (t10 * t308 + t143 * t2 - t144 * t3 - t309 * t9) * mrSges(7,3) + t366 * t40 + ((-t217 - t288) * t237 + ((-mrSges(5,1) * t231 + mrSges(5,2) * t230 - mrSges(4,1)) * qJD(3) + t277) * t234) * t221 + (-t167 / 0.2e1 - t190 / 0.2e1) * t62 + (t72 / 0.2e1 - t98 / 0.2e1) * t35 + (-t168 / 0.2e1 - t345 / 0.2e1) * t63 + (-Ifges(6,5) * t345 - Ifges(6,6) * t190) * t321 + (-Ifges(6,1) * t345 - Ifges(6,4) * t190) * t328 + (-Ifges(6,4) * t345 - Ifges(6,2) * t190) * t330 + (-t15 * t250 - t16 * t208 + t278 * t48 - t279 * t49) * mrSges(6,3) + (Ifges(6,4) * t208 - Ifges(6,2) * t250) * t336 + (Ifges(6,1) * t208 - Ifges(6,4) * t250) * t337 + t162 * (mrSges(6,1) * t250 + mrSges(6,2) * t208) + (((t299 + (-Ifges(4,4) / 0.2e1 + t248) * t234) * qJD(1) + ((Ifges(5,2) * t231 + t306) * t320 + (Ifges(5,1) * t230 + t305) * t319 - Ifges(4,5) / 0.2e1) * qJD(3) + t344) * t234 + (-t240 + t264 + (-t300 + t307 / 0.2e1 + (-Ifges(5,3) / 0.2e1 + Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t234) * qJD(1)) * t237 + (Ifges(5,5) * t230 + Ifges(6,5) * t208 + Ifges(7,5) * t144 + Ifges(5,6) * t231 - Ifges(6,6) * t250 + Ifges(7,6) * t143) * t346 / 0.2e1) * qJD(1) - t250 * t45 / 0.2e1 + t367 * t107 + (-t148 * t171 + t15 * t158 + t157 * t16 + t162 * t228 + t367 * t49 + t368 * t48) * m(6) + t368 * t108 + (qJ(4) * t181 + qJD(4) * t165 + t101 * mrSges(5,3) + t155 / 0.2e1) * t231 + (-pkin(3) * t209 + qJ(4) * t253 + qJD(4) * t251 - t128 * t151 - t129 * t152 - t185 * t212) * m(5) + t369 * t58 + t370 * t59 + (t10 * t369 + t169 * t60 + t2 * t55 + t3 * t54 + t366 * t80 + t370 * t9) * m(7) + t54 * t22 + t55 * t23 + (Ifges(6,5) * t168 + Ifges(6,6) * t167) * t322 + (Ifges(7,5) * t71 + Ifges(7,6) * t72) * t323 + (Ifges(7,5) * t99 + Ifges(7,6) * t98) * t324 + (Ifges(6,1) * t168 + Ifges(6,4) * t167) * t329 + (Ifges(6,4) * t168 + Ifges(6,2) * t167) * t331 + (Ifges(7,1) * t71 + Ifges(7,4) * t72) * t338 + (Ifges(7,1) * t99 + Ifges(7,4) * t98) * t339 + (Ifges(7,4) * t71 + Ifges(7,2) * t72) * t340 + (Ifges(7,4) * t99 + Ifges(7,2) * t98) * t341 + (Ifges(7,4) * t144 + Ifges(7,2) * t143) * t342 + (Ifges(7,1) * t144 + Ifges(7,4) * t143) * t343 + t143 * t6 / 0.2e1 + t144 * t7 / 0.2e1 + t60 * (-mrSges(7,1) * t143 + mrSges(7,2) * t144) + t157 * t76 + t158 * t77 - t152 * t165 - t151 * t166 + t169 * t8 - pkin(3) * t170 - t171 * t75 + t208 * t46 / 0.2e1 + t228 * t47; -t260 * t107 + t141 * t108 - t202 * t165 + t203 * t166 - t365 * t58 + t70 * t59 + (m(5) * t221 + t247) * t275 - m(5) * (-t128 * t203 + t129 * t202) + (-t10 * t365 + t70 * t9 + t60) * m(7) + (t141 * t48 - t260 * t49 + t162) * m(6) - t361; -m(7) * (t10 * t13 + t12 * t9) + (-Ifges(6,2) * t141 + t136 + t63) * t331 + t347 + (t141 * t49 + t260 * t48) * mrSges(6,3) + (Ifges(6,1) * t260 - t304) * t329 - t148 * (mrSges(6,1) * t141 + mrSges(6,2) * t260) + (Ifges(6,5) * t260 - Ifges(6,6) * t141) * t322 + t70 * t363 + t268 + (-t141 * t40 + t235 * t22 + t232 * t23 + (-t232 * t59 + t235 * t58) * qJD(6) + (-t141 * t80 + t2 * t232 + t235 * t3 + (t10 * t235 - t232 * t9) * qJD(6)) * m(7)) * pkin(5) - t13 * t58 - t12 * t59 + t62 * t328 - t48 * t107 + t49 * t108 + t372; t10 * t59 + t35 * t338 - t9 * t58 + t372;];
tauc  = t1(:);
