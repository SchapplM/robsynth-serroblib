% Calculate vector of centrifugal and coriolis load on the joints for
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:09:00
% EndTime: 2018-11-23 16:09:09
% DurationCPUTime: 9.37s
% Computational Cost: add. (9252->590), mult. (20945->846), div. (0->0), fcn. (14494->8), ass. (0->258)
t230 = sin(pkin(10));
t231 = cos(pkin(10));
t233 = sin(qJ(5));
t236 = cos(qJ(5));
t208 = t230 * t236 + t231 * t233;
t188 = t208 * qJD(1);
t234 = sin(qJ(3));
t167 = t234 * t188;
t190 = t208 * qJD(5);
t280 = t167 + t190;
t283 = t231 * t236;
t250 = t230 * t233 - t283;
t351 = t250 * t234;
t168 = qJD(1) * t351;
t346 = qJD(5) * t250;
t279 = t168 + t346;
t237 = cos(qJ(3));
t254 = pkin(3) * t237 + qJ(4) * t234;
t210 = t254 * qJD(1);
t238 = -pkin(1) - pkin(7);
t221 = qJD(1) * t238 + qJD(2);
t285 = t230 * t237;
t151 = t231 * t210 - t221 * t285;
t284 = t231 * t234;
t270 = pkin(8) * t284;
t244 = (pkin(4) * t237 + t270) * qJD(1);
t109 = t151 + t244;
t282 = t231 * t237;
t152 = t230 * t210 + t221 * t282;
t229 = t234 * qJD(1);
t267 = t230 * t229;
t125 = pkin(8) * t267 + t152;
t311 = pkin(8) + qJ(4);
t214 = t311 * t230;
t215 = t311 * t231;
t158 = -t233 * t214 + t236 * t215;
t369 = -t208 * qJD(4) - qJD(5) * t158 - t236 * t109 + t125 * t233;
t273 = qJD(5) * t236;
t368 = qJD(4) * t283 - t236 * t125 - t214 * t273 + (-qJD(4) * t230 - qJD(5) * t215 - t109) * t233;
t277 = qJD(1) * t237;
t375 = -pkin(5) * t277 + pkin(9) * t279 + t369;
t374 = pkin(9) * t280 - t368;
t179 = t250 * t237;
t350 = -qJD(3) * t179 - t190 * t234 - t188;
t177 = t208 * t237;
t349 = t250 * qJD(1) - qJD(3) * t177 + t234 * t346;
t232 = sin(qJ(6));
t235 = cos(qJ(6));
t202 = t231 * qJD(3) - t230 * t277;
t272 = t230 * qJD(3);
t203 = t231 * t277 + t272;
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
t292 = t235 * t39;
t227 = t229 + qJD(5);
t141 = t202 * t233 + t203 * t236;
t48 = -t233 * t95 + t236 * t93;
t38 = -pkin(9) * t141 + t48;
t37 = pkin(5) * t227 + t38;
t10 = t232 * t37 + t292;
t271 = qJD(1) * qJD(3);
t263 = t237 * t271;
t366 = -t141 * t232 + t235 * t260;
t242 = qJD(3) * t351;
t90 = qJD(1) * t242 + qJD(5) * t260;
t276 = qJD(3) * t234;
t243 = t208 * t276;
t91 = qJD(1) * t243 - qJD(5) * t141;
t28 = qJD(6) * t366 + t232 * t91 + t235 * t90;
t70 = t141 * t235 + t232 * t260;
t29 = -qJD(6) * t70 - t232 * t90 + t235 * t91;
t269 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t263;
t318 = Ifges(7,4) * t70;
t220 = qJD(6) + t227;
t325 = -t220 / 0.2e1;
t340 = -t70 / 0.2e1;
t342 = -t366 / 0.2e1;
t180 = qJD(3) * t254 - qJD(4) * t237 + qJD(2);
t159 = t180 * qJD(1);
t286 = t221 * t237;
t183 = (qJD(4) + t286) * qJD(3);
t100 = t231 * t159 - t183 * t230;
t81 = qJD(3) * t244 + t100;
t101 = t230 * t159 + t231 * t183;
t266 = t234 * t272;
t259 = qJD(1) * t266;
t94 = pkin(8) * t259 + t101;
t16 = -qJD(5) * t49 - t233 * t94 + t236 * t81;
t11 = pkin(5) * t263 - pkin(9) * t90 + t16;
t274 = qJD(5) * t233;
t15 = t233 * t81 + t236 * t94 + t93 * t273 - t274 * t95;
t14 = pkin(9) * t91 + t15;
t293 = t232 * t39;
t9 = t235 * t37 - t293;
t2 = qJD(6) * t9 + t11 * t232 + t14 * t235;
t3 = -qJD(6) * t10 + t11 * t235 - t14 * t232;
t359 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t64 = Ifges(7,4) * t366;
t36 = Ifges(7,1) * t70 + Ifges(7,5) * t220 + t64;
t185 = -qJD(3) * pkin(3) + qJD(4) - t286;
t148 = -pkin(4) * t202 + t185;
t80 = -pkin(5) * t260 + t148;
t373 = t269 + t359 + (Ifges(7,5) * t366 - Ifges(7,6) * t70) * t325 + (t10 * t70 + t366 * t9) * mrSges(7,3) + (-Ifges(7,2) * t70 + t36 + t64) * t342 - t80 * (mrSges(7,1) * t70 + mrSges(7,2) * t366) + (Ifges(7,1) * t366 - t318) * t340;
t372 = -qJD(3) / 0.2e1;
t157 = -t236 * t214 - t215 * t233;
t117 = -pkin(9) * t208 + t157;
t118 = -pkin(9) * t250 + t158;
t55 = t117 * t232 + t118 * t235;
t371 = -qJD(6) * t55 + t374 * t232 + t235 * t375;
t54 = t117 * t235 - t118 * t232;
t370 = qJD(6) * t54 + t232 * t375 - t374 * t235;
t171 = -pkin(4) * t267 + t212;
t367 = pkin(5) * t280 - t171;
t35 = Ifges(7,2) * t366 + Ifges(7,6) * t220 + t318;
t364 = t35 / 0.2e1;
t264 = Ifges(4,6) * t372;
t363 = Ifges(4,5) * t372;
t47 = -t91 * mrSges(6,1) + t90 * mrSges(6,2);
t8 = -t29 * mrSges(7,1) + t28 * mrSges(7,2);
t362 = -t47 - t8;
t360 = -qJD(1) / 0.2e1;
t176 = t208 * t234;
t111 = -t176 * t235 + t232 * t351;
t355 = qJD(6) * t111 + t349 * t232 + t235 * t350;
t113 = -t176 * t232 - t235 * t351;
t354 = -qJD(6) * t113 - t232 * t350 + t349 * t235;
t353 = qJ(2) * (m(3) + m(4));
t200 = t231 * t213;
t262 = -t230 * t238 + pkin(4);
t137 = -pkin(8) * t282 + t234 * t262 + t200;
t281 = t234 * t238;
t161 = t230 * t213 + t231 * t281;
t150 = -pkin(8) * t285 + t161;
t74 = t233 * t137 + t236 * t150;
t348 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t347 = qJD(3) * t237;
t252 = t128 * t231 + t129 * t230;
t303 = Ifges(5,2) * t230;
t306 = Ifges(5,4) * t231;
t255 = t303 - t306;
t307 = Ifges(5,4) * t230;
t256 = -Ifges(5,1) * t231 + t307;
t257 = -mrSges(5,1) * t230 - mrSges(5,2) * t231;
t320 = -t231 / 0.2e1;
t321 = t230 / 0.2e1;
t345 = -t252 * mrSges(5,3) - t185 * t257 - t363 - (Ifges(4,1) * t237 - Ifges(4,4) * t234) * t360 - t202 * t255 / 0.2e1 - t203 * t256 / 0.2e1 - (t203 * Ifges(5,4) + t202 * Ifges(5,2) + Ifges(5,6) * t229) * t321 - (t203 * Ifges(5,1) + t202 * Ifges(5,4) + Ifges(5,5) * t229) * t320;
t344 = t28 / 0.2e1;
t343 = t29 / 0.2e1;
t341 = t366 / 0.2e1;
t339 = t70 / 0.2e1;
t338 = t90 / 0.2e1;
t337 = t91 / 0.2e1;
t112 = -t177 * t235 + t179 * t232;
t334 = t112 / 0.2e1;
t114 = -t177 * t232 - t179 * t235;
t333 = t114 / 0.2e1;
t332 = -t260 / 0.2e1;
t331 = t260 / 0.2e1;
t330 = -t141 / 0.2e1;
t329 = t141 / 0.2e1;
t328 = (Ifges(5,5) * t237 + t234 * t256) * t271 / 0.2e1;
t327 = -t177 / 0.2e1;
t326 = -t179 / 0.2e1;
t324 = t220 / 0.2e1;
t323 = -t227 / 0.2e1;
t322 = t227 / 0.2e1;
t143 = -t208 * t232 - t235 * t250;
t71 = qJD(6) * t143 - t190 * t232 - t235 * t346;
t99 = t167 * t232 + t168 * t235;
t310 = t71 - t99;
t144 = t208 * t235 - t232 * t250;
t72 = -qJD(6) * t144 - t190 * t235 + t232 * t346;
t98 = t167 * t235 - t168 * t232;
t309 = t72 - t98;
t308 = Ifges(4,4) * t237;
t305 = Ifges(6,4) * t141;
t304 = Ifges(5,5) * t231;
t302 = Ifges(5,6) * t230;
t301 = qJ(2) * mrSges(4,1);
t300 = qJ(2) * mrSges(4,2);
t289 = qJD(3) * mrSges(4,2);
t275 = qJD(3) * t238;
t265 = t237 * t275;
t147 = t230 * t180 + t231 * t265;
t278 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t202 - mrSges(5,2) * t203 - mrSges(4,3) * t277;
t268 = Ifges(6,5) * t90 + Ifges(6,6) * t91 + Ifges(6,3) * t263;
t228 = -pkin(4) * t231 - pkin(3);
t209 = t221 * t276;
t73 = t236 * t137 - t150 * t233;
t201 = pkin(4) * t285 - t237 * t238;
t258 = mrSges(4,1) * t234 + mrSges(4,2) * t237;
t52 = pkin(5) * t234 + pkin(9) * t179 + t73;
t53 = -pkin(9) * t177 + t74;
t30 = -t232 * t53 + t235 * t52;
t31 = t232 * t52 + t235 * t53;
t253 = -t100 * t230 + t101 * t231;
t251 = -t128 * t230 + t129 * t231;
t186 = -pkin(4) * t266 + t234 * t275;
t248 = t304 / 0.2e1 - t302 / 0.2e1;
t247 = t257 * qJD(1);
t162 = -pkin(4) * t259 + t209;
t164 = t231 * t180;
t104 = t164 + (t237 * t262 + t270) * qJD(3);
t123 = pkin(8) * t266 + t147;
t32 = t233 * t104 + t236 * t123 + t137 * t273 - t150 * t274;
t33 = -qJD(5) * t74 + t236 * t104 - t123 * t233;
t240 = t128 * mrSges(5,1) + t48 * mrSges(6,1) + t9 * mrSges(7,1) + Ifges(5,3) * t229 / 0.2e1 + t203 * Ifges(5,5) + t202 * Ifges(5,6) + t264 + (-t234 * Ifges(4,2) + t308) * t360 + t220 * Ifges(7,3) + t70 * Ifges(7,5) + t366 * Ifges(7,6) + t227 * Ifges(6,3) + t141 * Ifges(6,5) + t260 * Ifges(6,6) - t10 * mrSges(7,2) - t129 * mrSges(5,2) - t49 * mrSges(6,2);
t217 = -mrSges(4,3) * t229 - t289;
t211 = t258 * qJD(1);
t182 = (mrSges(5,1) * t237 + mrSges(5,3) * t284) * t271;
t181 = (mrSges(5,3) * t230 * t234 - mrSges(5,2) * t237) * t271;
t170 = t247 * t276;
t169 = pkin(5) * t250 + t228;
t166 = mrSges(5,1) * t229 - mrSges(5,3) * t203;
t165 = -mrSges(5,2) * t229 + mrSges(5,3) * t202;
t160 = -t230 * t281 + t200;
t155 = (Ifges(5,6) * t237 + t234 * t255) * t271;
t146 = -t230 * t265 + t164;
t145 = pkin(5) * t177 + t201;
t136 = Ifges(6,4) * t260;
t122 = t237 * t346 + t243;
t120 = -t190 * t237 + t242;
t108 = mrSges(6,1) * t227 - mrSges(6,3) * t141;
t107 = -mrSges(6,2) * t227 + mrSges(6,3) * t260;
t92 = -pkin(5) * t122 + t186;
t77 = -mrSges(6,2) * t263 + mrSges(6,3) * t91;
t76 = mrSges(6,1) * t263 - mrSges(6,3) * t90;
t75 = -mrSges(6,1) * t260 + mrSges(6,2) * t141;
t63 = Ifges(6,1) * t141 + Ifges(6,5) * t227 + t136;
t62 = Ifges(6,2) * t260 + Ifges(6,6) * t227 + t305;
t60 = -pkin(5) * t91 + t162;
t59 = mrSges(7,1) * t220 - mrSges(7,3) * t70;
t58 = -mrSges(7,2) * t220 + mrSges(7,3) * t366;
t46 = t90 * Ifges(6,1) + t91 * Ifges(6,4) + Ifges(6,5) * t263;
t45 = t90 * Ifges(6,4) + t91 * Ifges(6,2) + Ifges(6,6) * t263;
t44 = -qJD(6) * t114 - t120 * t232 + t122 * t235;
t42 = qJD(6) * t112 + t120 * t235 + t122 * t232;
t40 = -mrSges(7,1) * t366 + mrSges(7,2) * t70;
t24 = pkin(9) * t122 + t32;
t23 = -mrSges(7,2) * t263 + mrSges(7,3) * t29;
t22 = mrSges(7,1) * t263 - mrSges(7,3) * t28;
t21 = pkin(5) * t347 - pkin(9) * t120 + t33;
t13 = t235 * t38 - t293;
t12 = -t232 * t38 - t292;
t7 = t28 * Ifges(7,1) + t29 * Ifges(7,4) + Ifges(7,5) * t263;
t6 = t28 * Ifges(7,4) + t29 * Ifges(7,2) + Ifges(7,6) * t263;
t5 = -qJD(6) * t31 + t21 * t235 - t232 * t24;
t4 = qJD(6) * t30 + t21 * t232 + t235 * t24;
t1 = [m(5) * (t100 * t160 + t101 * t161 + t128 * t146 + t129 * t147) + (t10 * t44 + t112 * t2 - t114 * t3 - t42 * t9) * mrSges(7,3) + (-Ifges(6,1) * t179 - Ifges(6,4) * t177) * t338 + (-Ifges(6,4) * t179 - Ifges(6,2) * t177) * t337 + (-t120 * t48 + t122 * t49 - t15 * t177 + t16 * t179) * mrSges(6,3) + t162 * (mrSges(6,1) * t177 - mrSges(6,2) * t179) + (Ifges(7,4) * t114 + Ifges(7,2) * t112) * t343 + (t231 * t328 - t230 * t155 / 0.2e1 - t238 * t170 + (-t100 * t231 - t101 * t230) * mrSges(5,3)) * t237 + (Ifges(7,1) * t114 + Ifges(7,4) * t112) * t344 + (t269 + t268) * t234 / 0.2e1 + ((-m(5) * t238 - t257) * t212 + t240 + t238 * t217 + t264) * t347 + t30 * t22 + t31 * t23 + m(6) * (t148 * t186 + t15 * t74 + t16 * t73 + t162 * t201 + t32 * t49 + t33 * t48) + m(7) * (t10 * t4 + t145 * t60 + t2 * t31 + t3 * t30 + t5 * t9 + t80 * t92) + t42 * t36 / 0.2e1 + t4 * t58 + t5 * t59 + t44 * t364 + t73 * t76 + t74 * t77 + t80 * (-mrSges(7,1) * t44 + mrSges(7,2) * t42) + (((2 * mrSges(3,3)) + t258 + 0.2e1 * t353) * qJD(2) + (0.2e1 * t301 + Ifges(6,5) * t326 + Ifges(6,6) * t327 + Ifges(7,5) * t333 + Ifges(7,6) * t334 + (-0.3e1 / 0.2e1 * Ifges(4,4) + t248) * t237) * t347) * qJD(1) + t92 * t40 + t32 * t107 + t33 * t108 + t60 * (-mrSges(7,1) * t112 + mrSges(7,2) * t114) + t120 * t63 / 0.2e1 + t122 * t62 / 0.2e1 + (t100 * mrSges(5,1) - t101 * mrSges(5,2) + Ifges(6,5) * t338 + Ifges(7,5) * t344 + Ifges(6,6) * t337 + Ifges(7,6) * t343 + (t363 + (m(5) * t185 - t278) * t238 + (-0.2e1 * t300 + (-0.3e1 / 0.2e1 * t304 + 0.3e1 / 0.2e1 * t302 + 0.3e1 / 0.2e1 * Ifges(4,4)) * t234 + (-Ifges(5,1) * t231 ^ 2 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(5,3) + Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1 + (t306 - t303 / 0.2e1) * t230) * t237) * qJD(1) - t345) * qJD(3) + t348 + t359) * t234 + t145 * t8 + t148 * (-mrSges(6,1) * t122 + mrSges(6,2) * t120) + t147 * t165 + t146 * t166 + t161 * t181 + t160 * t182 + t186 * t75 + t201 * t47 + qJD(2) * t211 + (Ifges(6,5) * t120 + Ifges(6,6) * t122) * t322 + (Ifges(7,5) * t42 + Ifges(7,6) * t44) * t324 + t46 * t326 + t45 * t327 + (Ifges(6,1) * t120 + Ifges(6,4) * t122) * t329 + (Ifges(6,4) * t120 + Ifges(6,2) * t122) * t331 + t7 * t333 + t6 * t334 + (Ifges(7,1) * t42 + Ifges(7,4) * t44) * t339 + (Ifges(7,4) * t42 + Ifges(7,2) * t44) * t341; t111 * t22 + t113 * t23 - t176 * t76 - t351 * t77 + t354 * t59 + t355 * t58 + (t231 * t181 - t230 * t182) * t234 + t349 * t108 + t350 * t107 + (-t170 + t362) * t237 + ((t165 * t231 - t166 * t230 + t217) * t237 + (t40 + t75 - t278) * t234) * qJD(3) + m(5) * (t253 * t234 + (t185 * t234 + (t251 - t212) * t237) * qJD(3)) + (t10 * t355 + t111 * t3 + t113 * t2 - t237 * t60 + t276 * t80 + t354 * t9) * m(7) + (t148 * t276 - t15 * t351 - t16 * t176 - t162 * t237 + t349 * t48 + t350 * t49) * m(6) + (-m(5) * t252 - t165 * t230 - t166 * t231 - t211 + (-mrSges(3,3) - t353) * qJD(1)) * qJD(1); (-t167 / 0.2e1 - t190 / 0.2e1) * t62 + (-mrSges(7,1) * t309 + mrSges(7,2) * t310) * t80 + (t10 * t309 + t143 * t2 - t144 * t3 - t310 * t9) * mrSges(7,3) + (t71 / 0.2e1 - t99 / 0.2e1) * t36 + (qJ(4) * t181 + qJD(4) * t165 + t101 * mrSges(5,3) + t155 / 0.2e1) * t231 + (mrSges(6,1) * t280 - mrSges(6,2) * t279) * t148 + (-t100 * mrSges(5,3) - qJ(4) * t182 - qJD(4) * t166 + t328) * t230 + ((-t217 - t289) * t237 + ((-mrSges(5,1) * t231 + mrSges(5,2) * t230 - mrSges(4,1)) * qJD(3) + t278) * t234) * t221 + (-Ifges(6,1) * t346 - Ifges(6,4) * t190) * t329 + (-Ifges(6,4) * t346 - Ifges(6,2) * t190) * t331 + (-t168 / 0.2e1 - t346 / 0.2e1) * t63 + (-Ifges(6,5) * t346 - Ifges(6,6) * t190) * t322 + (-t15 * t250 - t16 * t208 + t279 * t48 - t280 * t49) * mrSges(6,3) + t162 * (mrSges(6,1) * t250 + mrSges(6,2) * t208) + (Ifges(6,4) * t208 - Ifges(6,2) * t250) * t337 + (Ifges(6,1) * t208 - Ifges(6,4) * t250) * t338 + (((t300 + (-Ifges(4,4) / 0.2e1 + t248) * t234) * qJD(1) + (-Ifges(4,5) / 0.2e1 + (Ifges(5,1) * t230 + t306) * t320 + (Ifges(5,2) * t231 + t307) * t321) * qJD(3) + t345) * t234 + (-t240 + (-t301 + t308 / 0.2e1 + (-Ifges(5,3) / 0.2e1 + Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t234) * qJD(1) + t264) * t237 + (Ifges(5,5) * t230 + Ifges(6,5) * t208 + Ifges(7,5) * t144 + Ifges(5,6) * t231 - Ifges(6,6) * t250 + Ifges(7,6) * t143) * t347 / 0.2e1) * qJD(1) - t250 * t45 / 0.2e1 + (t72 / 0.2e1 - t98 / 0.2e1) * t35 + (-pkin(3) * t209 + qJ(4) * t253 + qJD(4) * t251 - t128 * t151 - t129 * t152 - t185 * t212) * m(5) + t370 * t58 + t371 * t59 + (t10 * t370 + t169 * t60 + t2 * t55 + t3 * t54 + t367 * t80 + t371 * t9) * m(7) + t54 * t22 + t55 * t23 + t143 * t6 / 0.2e1 + t144 * t7 / 0.2e1 + t60 * (-mrSges(7,1) * t143 + mrSges(7,2) * t144) + t157 * t76 + t158 * t77 - t152 * t165 - t151 * t166 + t169 * t8 - pkin(3) * t170 - t171 * t75 + t208 * t46 / 0.2e1 + t228 * t47 + t367 * t40 + t368 * t107 + t369 * t108 + (-t148 * t171 + t15 * t158 + t157 * t16 + t162 * t228 + t368 * t49 + t369 * t48) * m(6) + (Ifges(6,5) * t168 + Ifges(6,6) * t167) * t323 + (Ifges(7,5) * t71 + Ifges(7,6) * t72) * t324 + (Ifges(7,5) * t99 + Ifges(7,6) * t98) * t325 + (Ifges(6,1) * t168 + Ifges(6,4) * t167) * t330 + (Ifges(6,4) * t168 + Ifges(6,2) * t167) * t332 + (Ifges(7,1) * t71 + Ifges(7,4) * t72) * t339 + (Ifges(7,1) * t99 + Ifges(7,4) * t98) * t340 + (Ifges(7,4) * t71 + Ifges(7,2) * t72) * t341 + (Ifges(7,4) * t99 + Ifges(7,2) * t98) * t342 + (Ifges(7,4) * t144 + Ifges(7,2) * t143) * t343 + (Ifges(7,1) * t144 + Ifges(7,4) * t143) * t344; -t260 * t107 + t141 * t108 - t202 * t165 + t203 * t166 - t366 * t58 + t70 * t59 + (m(5) * t221 + t247) * t276 - m(5) * (-t128 * t203 + t129 * t202) + (-t10 * t366 + t70 * t9 + t60) * m(7) + (t141 * t48 - t260 * t49 + t162) * m(6) - t362; t268 - m(7) * (t10 * t13 + t12 * t9) + t348 + (-Ifges(6,2) * t141 + t136 + t63) * t332 + (-t141 * t40 + t235 * t22 + t232 * t23 + (-t232 * t59 + t235 * t58) * qJD(6) + (-t141 * t80 + t2 * t232 + t235 * t3 + (t10 * t235 - t232 * t9) * qJD(6)) * m(7)) * pkin(5) + t70 * t364 + (t141 * t49 + t260 * t48) * mrSges(6,3) + (Ifges(6,1) * t260 - t305) * t330 - t148 * (mrSges(6,1) * t141 + mrSges(6,2) * t260) + (Ifges(6,5) * t260 - Ifges(6,6) * t141) * t323 - t13 * t58 - t12 * t59 - t48 * t107 + t49 * t108 + t62 * t329 + t373; t10 * t59 + t35 * t339 - t9 * t58 + t373;];
tauc  = t1(:);
