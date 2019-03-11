% Calculate vector of inverse dynamics joint torques for
% S5RRRRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:56
% EndTime: 2019-03-08 18:37:09
% DurationCPUTime: 7.84s
% Computational Cost: add. (7482->512), mult. (17157->710), div. (0->0), fcn. (13414->14), ass. (0->250)
t185 = qJ(2) + qJ(3);
t181 = qJ(4) + t185;
t168 = sin(t181);
t169 = cos(t181);
t188 = sin(qJ(5));
t320 = mrSges(6,2) * t188;
t391 = t168 * t320 + t169 * (m(6) * pkin(6) + mrSges(6,3));
t190 = sin(qJ(3));
t191 = sin(qJ(2));
t195 = cos(qJ(3));
t196 = cos(qJ(2));
t149 = t190 * t191 - t195 * t196;
t142 = t149 * qJD(1);
t285 = qJD(1) * t196;
t162 = qJD(1) * pkin(1) + pkin(2) * t285;
t113 = -pkin(3) * t142 + t162;
t184 = qJD(2) + qJD(3);
t313 = pkin(2) * qJD(2);
t270 = t195 * t313;
t154 = pkin(3) * t184 + t270;
t194 = cos(qJ(4));
t189 = sin(qJ(4));
t272 = t190 * t313;
t248 = t189 * t272;
t114 = t154 * t194 - t248;
t178 = qJD(4) + t184;
t111 = -pkin(4) * t178 - t114;
t193 = cos(qJ(5));
t241 = mrSges(6,1) * t188 + mrSges(6,2) * t193;
t217 = t111 * t241;
t226 = t190 * t196 + t191 * t195;
t143 = t226 * qJD(1);
t230 = t142 * t189 - t194 * t143;
t83 = t178 * t193 - t188 * t230;
t82 = Ifges(6,4) * t83;
t84 = t178 * t188 + t193 * t230;
t251 = t194 * t142 + t143 * t189;
t96 = qJD(5) - t251;
t33 = t84 * Ifges(6,1) + t96 * Ifges(6,5) + t82;
t300 = t193 * t33;
t307 = t178 * Ifges(5,5);
t312 = t114 * mrSges(5,3);
t327 = t84 * Ifges(6,4);
t32 = t83 * Ifges(6,2) + t96 * Ifges(6,6) + t327;
t347 = t188 / 0.2e1;
t95 = Ifges(5,4) * t251;
t60 = Ifges(5,1) * t230 + t307 + t95;
t390 = -t300 / 0.2e1 + t32 * t347 + t312 - t60 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t230 - t217 - t95 / 0.2e1 - t113 * mrSges(5,2) - t307 / 0.2e1;
t306 = t178 * Ifges(5,6);
t325 = t96 * Ifges(6,3);
t326 = t84 * Ifges(6,5);
t328 = t83 * Ifges(6,6);
t31 = t325 + t326 + t328;
t115 = t154 * t189 + t194 * t272;
t311 = t115 * mrSges(5,3);
t316 = Ifges(5,4) * t230;
t112 = pkin(6) * t178 + t115;
t49 = -pkin(4) * t251 - pkin(6) * t230 + t113;
t40 = t112 * t193 + t188 * t49;
t378 = t40 * mrSges(6,2);
t39 = -t112 * t188 + t193 * t49;
t379 = t39 * mrSges(6,1);
t59 = Ifges(5,2) * t251 + t306 + t316;
t389 = t378 + t311 + t59 / 0.2e1 - t31 / 0.2e1 + t316 / 0.2e1 - t113 * mrSges(5,1) + t306 / 0.2e1 - t379;
t179 = sin(t185);
t180 = cos(t185);
t321 = mrSges(5,2) * t169;
t386 = mrSges(4,1) * t179 + mrSges(5,1) * t168 + mrSges(4,2) * t180 + t321;
t192 = sin(qJ(1));
t197 = cos(qJ(1));
t385 = g(1) * t197 + g(2) * t192;
t164 = t168 * mrSges(5,2);
t384 = mrSges(4,1) * t180 + t169 * mrSges(5,1) - t179 * mrSges(4,2) + t168 * mrSges(6,3) - t164;
t346 = -t196 / 0.2e1;
t377 = -mrSges(5,1) * t178 - mrSges(6,1) * t83 + mrSges(6,2) * t84 + mrSges(5,3) * t230;
t66 = -mrSges(5,1) * t251 + mrSges(5,2) * t230;
t376 = m(5) * t113 + t66;
t186 = qJDD(1) * pkin(1);
t317 = Ifges(3,4) * t196;
t318 = Ifges(3,4) * t191;
t375 = t191 * (-Ifges(3,1) * t196 + t318) + t196 * (Ifges(3,2) * t191 - t317);
t108 = t149 * t189 - t194 * t226;
t279 = qJD(5) * t193;
t213 = t149 * qJD(3);
t109 = qJD(2) * t149 + t213;
t214 = t226 * qJD(3);
t110 = qJD(2) * t226 + t214;
t229 = t194 * t149 + t189 * t226;
t52 = qJD(4) * t229 + t109 * t194 + t110 * t189;
t222 = t108 * t279 + t188 * t52;
t280 = qJD(5) * t188;
t374 = -t39 * t279 - t40 * t280;
t55 = -mrSges(6,2) * t96 + mrSges(6,3) * t83;
t56 = mrSges(6,1) * t96 - mrSges(6,3) * t84;
t373 = -t188 * t55 - t193 * t56;
t183 = qJDD(2) + qJDD(3);
t177 = qJDD(4) + t183;
t278 = qJD(1) * qJD(2);
t152 = -qJDD(1) * t196 + t191 * t278;
t153 = -qJDD(1) * t191 - t196 * t278;
t87 = qJD(1) * t213 + t152 * t190 + t153 * t195;
t88 = qJD(1) * t214 + t152 * t195 - t153 * t190;
t37 = qJD(4) * t251 + t189 * t88 + t194 * t87;
t18 = -qJD(5) * t84 + t177 * t193 - t188 * t37;
t38 = -qJD(4) * t230 - t189 * t87 + t194 * t88;
t36 = qJDD(5) - t38;
t10 = -mrSges(6,2) * t36 + mrSges(6,3) * t18;
t17 = qJD(5) * t83 + t177 * t188 + t193 * t37;
t9 = mrSges(6,1) * t36 - mrSges(6,3) * t17;
t372 = t193 * t10 - t188 * t9;
t89 = -mrSges(5,2) * t178 + mrSges(5,3) * t251;
t370 = -t188 * t56 + t193 * t55 + t89;
t157 = t169 * t320;
t322 = mrSges(6,1) * t193;
t368 = t169 * t322 - t157 + t384;
t367 = mrSges(3,3) + mrSges(2,2) + mrSges(5,3) + mrSges(4,3);
t343 = pkin(2) * t190;
t269 = qJD(3) * t343;
t341 = pkin(2) * t195;
t145 = -qJD(2) * t269 + qJDD(2) * t341;
t122 = pkin(3) * t183 + t145;
t283 = qJD(3) * t195;
t146 = (qJD(2) * t283 + qJDD(2) * t190) * pkin(2);
t73 = -qJD(4) * t115 + t122 * t194 - t146 * t189;
t133 = -pkin(2) * t152 + t186;
t71 = -pkin(3) * t88 + t133;
t11 = -pkin(4) * t38 - pkin(6) * t37 + t71;
t281 = qJD(4) * t194;
t72 = -qJD(4) * t248 + t189 * t122 + t194 * t146 + t154 * t281;
t68 = pkin(6) * t177 + t72;
t2 = qJD(5) * t39 + t11 * t188 + t193 * t68;
t3 = -qJD(5) * t40 + t11 * t193 - t188 * t68;
t366 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t234 = t188 * t40 + t193 * t39;
t365 = m(6) * t234 - t373;
t67 = pkin(4) * t230 - pkin(6) * t251;
t329 = t3 * t188;
t203 = -qJD(5) * t234 - t329;
t330 = t193 * t2;
t364 = m(6) * (t203 + t330) - t56 * t279 - t55 * t280 + t372;
t182 = t196 * pkin(2);
t338 = pkin(3) * t180;
t256 = t182 + t338;
t151 = pkin(1) + t256;
t161 = mrSges(3,1) * t196 - mrSges(3,2) * t191;
t245 = pkin(4) * t169 + pkin(6) * t168;
t173 = t182 + pkin(1);
t345 = m(4) * t173;
t363 = m(3) * pkin(1) + m(5) * t151 + m(6) * (t151 + t245) + mrSges(2,1) + t161 + t345 + t384;
t361 = m(6) * pkin(4);
t360 = t17 / 0.2e1;
t359 = t18 / 0.2e1;
t356 = t36 / 0.2e1;
t353 = -t83 / 0.2e1;
t352 = -t84 / 0.2e1;
t351 = t84 / 0.2e1;
t350 = -t96 / 0.2e1;
t348 = -t143 / 0.2e1;
t342 = pkin(2) * t191;
t340 = pkin(3) * t143;
t339 = pkin(3) * t179;
t337 = pkin(3) * t189;
t336 = pkin(3) * t194;
t335 = pkin(4) * t168;
t332 = g(3) * t168;
t319 = mrSges(6,3) * t193;
t315 = Ifges(6,4) * t188;
t314 = Ifges(6,4) * t193;
t310 = t143 * Ifges(4,4);
t305 = t188 * mrSges(6,3);
t297 = t108 * t188;
t296 = t108 * t193;
t292 = t188 * t197;
t291 = t189 * t190;
t290 = t190 * t194;
t289 = t192 * t188;
t288 = t192 * t193;
t287 = t193 * t197;
t172 = pkin(3) + t341;
t137 = pkin(2) * t290 + t189 * t172;
t286 = qJD(1) * t191;
t284 = qJD(2) * t191;
t282 = qJD(4) * t189;
t276 = m(4) * t342;
t275 = Ifges(6,5) * t17 + Ifges(6,6) * t18 + Ifges(6,3) * t36;
t273 = pkin(2) * t286;
t271 = pkin(2) * t284;
t266 = t168 * t322;
t259 = t300 / 0.2e1;
t257 = -t280 / 0.2e1;
t255 = t162 * t276;
t254 = -t278 / 0.2e1;
t250 = mrSges(4,3) * t272;
t249 = mrSges(4,3) * t270;
t121 = -pkin(3) * t149 + t173;
t247 = t391 * t192;
t246 = t391 * t197;
t244 = mrSges(3,1) * t191 + mrSges(3,2) * t196;
t240 = -t191 * Ifges(3,1) - t317;
t239 = Ifges(6,1) * t193 - t315;
t238 = -t196 * Ifges(3,2) - t318;
t237 = -Ifges(6,2) * t188 + t314;
t236 = -Ifges(3,5) * t196 + Ifges(3,6) * t191;
t235 = Ifges(6,5) * t193 - Ifges(6,6) * t188;
t233 = -t188 * t39 + t193 * t40;
t228 = t189 * t195 + t290;
t227 = t194 * t195 - t291;
t225 = mrSges(5,1) + t322 + t361;
t101 = -pkin(3) * t110 - t271;
t136 = -pkin(2) * t291 + t172 * t194;
t223 = pkin(1) * t244;
t221 = t108 * t280 - t193 * t52;
t220 = t83 * t237;
t219 = t84 * t239;
t218 = t96 * t235;
t54 = -t340 + t67;
t210 = -t245 - t338;
t155 = -t339 - t342;
t204 = m(6) * (t155 - t335) - t266;
t202 = m(6) * (-t335 - t339) - t266;
t201 = t168 * t225 + t321;
t6 = t17 * Ifges(6,4) + t18 * Ifges(6,2) + t36 * Ifges(6,6);
t69 = -pkin(4) * t177 - t73;
t7 = t17 * Ifges(6,1) + t18 * Ifges(6,4) + t36 * Ifges(6,5);
t199 = -t72 * mrSges(5,2) + t2 * t319 + t7 * t347 + (Ifges(6,1) * t188 + t314) * t360 + (Ifges(6,2) * t193 + t315) * t359 + Ifges(5,3) * t177 + (Ifges(6,5) * t188 + Ifges(6,6) * t193) * t356 + t32 * t257 + Ifges(5,6) * t38 + Ifges(5,5) * t37 + t193 * t6 / 0.2e1 + t69 * (t320 - t322) + t73 * mrSges(5,1) + (t259 + t217) * qJD(5) + (t220 + t219 + t218) * qJD(5) / 0.2e1;
t134 = Ifges(4,4) * t142;
t93 = t142 * Ifges(4,2) + t184 * Ifges(4,6) - t310;
t94 = -t143 * Ifges(4,1) + t184 * Ifges(4,5) + t134;
t198 = -t162 * (-mrSges(4,1) * t143 + mrSges(4,2) * t142) - t3 * t305 - t143 * t250 + Ifges(4,3) * t183 - t184 * (Ifges(4,5) * t142 + Ifges(4,6) * t143) / 0.2e1 + t143 * (Ifges(4,1) * t142 + t310) / 0.2e1 + t142 * t249 + t199 + t374 * mrSges(6,3) + t145 * mrSges(4,1) - t146 * mrSges(4,2) + t93 * t348 + Ifges(4,5) * t87 + Ifges(4,6) * t88 - (Ifges(4,2) * t143 + t134 + t94) * t142 / 0.2e1 + (Ifges(6,5) * t352 + Ifges(6,6) * t353 + Ifges(6,3) * t350 + t389) * t230 + (t235 * t350 + t237 * t353 + t239 * t352 + t40 * t305 + t39 * t319 + t390) * t251;
t171 = -pkin(4) - t336;
t141 = Ifges(3,5) * qJD(2) + qJD(1) * t240;
t140 = Ifges(3,6) * qJD(2) + qJD(1) * t238;
t132 = t227 * t313;
t131 = t228 * t313;
t128 = t169 * t287 - t289;
t127 = -t169 * t292 - t288;
t126 = -t169 * t288 - t292;
t125 = t169 * t289 - t287;
t118 = mrSges(4,1) * t184 + mrSges(4,3) * t143;
t117 = -mrSges(4,2) * t184 + mrSges(4,3) * t142;
t105 = -mrSges(4,1) * t142 - mrSges(4,2) * t143;
t53 = qJD(4) * t108 + t109 * t189 - t194 * t110;
t48 = t114 * t193 + t188 * t67;
t47 = -t114 * t188 + t193 * t67;
t44 = t132 * t193 + t188 * t54;
t43 = -t132 * t188 + t193 * t54;
t27 = -mrSges(5,2) * t177 + mrSges(5,3) * t38;
t26 = mrSges(5,1) * t177 - mrSges(5,3) * t37;
t8 = -mrSges(6,1) * t18 + mrSges(6,2) * t17;
t1 = [t53 * t379 + t251 * (Ifges(5,4) * t52 - Ifges(5,2) * t53) / 0.2e1 + t230 * (Ifges(5,1) * t52 - Ifges(5,4) * t53) / 0.2e1 - t109 * t249 + t184 * (Ifges(4,5) * t109 + Ifges(4,6) * t110) / 0.2e1 - t52 * t312 + t111 * (mrSges(6,1) * t222 - mrSges(6,2) * t221) + t83 * (-Ifges(6,4) * t221 - Ifges(6,2) * t222 + Ifges(6,6) * t53) / 0.2e1 + t96 * (-Ifges(6,5) * t221 - Ifges(6,6) * t222 + Ifges(6,3) * t53) / 0.2e1 - t105 * t271 + (m(5) * t71 - mrSges(5,1) * t38 + mrSges(5,2) * t37) * t121 - t53 * t311 - t191 * (Ifges(3,1) * t153 + Ifges(3,4) * t152) / 0.2e1 + t161 * t186 + t140 * t284 / 0.2e1 + t7 * t296 / 0.2e1 - t6 * t297 / 0.2e1 - t223 * t278 + (Ifges(3,4) * t153 + Ifges(3,2) * t152) * t346 + t152 * t238 / 0.2e1 + t153 * t240 / 0.2e1 + t110 * t250 + Ifges(2,3) * qJDD(1) + (t71 * mrSges(5,2) - t73 * mrSges(5,3) + Ifges(5,1) * t37 + Ifges(5,4) * t38 + Ifges(5,5) * t177 + t235 * t356 + t237 * t359 + t239 * t360 + t241 * t69 + t257 * t33) * t108 + (t141 * t346 - t255 + t236 * qJD(2) / 0.2e1) * qJD(2) + t365 * (pkin(4) * t53 - pkin(6) * t52 + t101) + (-t126 * mrSges(6,1) - t125 * mrSges(6,2) + t192 * t363 + t197 * t367) * g(1) + (-t128 * mrSges(6,1) - t127 * mrSges(6,2) + t192 * t367 - t197 * t363) * g(2) + t109 * t94 / 0.2e1 + t110 * t93 / 0.2e1 + t113 * (mrSges(5,1) * t53 + mrSges(5,2) * t52) - t222 * t32 / 0.2e1 + t375 * t254 + (m(3) * t186 - mrSges(3,1) * t152 + mrSges(3,2) * t153) * pkin(1) + t376 * t101 + t142 * (Ifges(4,4) * t109 + Ifges(4,2) * t110) / 0.2e1 + t162 * (-mrSges(4,1) * t110 + mrSges(4,2) * t109) + t173 * (-mrSges(4,1) * t88 + mrSges(4,2) * t87) + t178 * (Ifges(5,5) * t52 - Ifges(5,6) * t53) / 0.2e1 + t52 * t259 + t133 * t345 + (Ifges(4,1) * t109 + Ifges(4,4) * t110) * t348 + (-Ifges(6,1) * t221 - Ifges(6,4) * t222 + Ifges(6,5) * t53) * t351 - t53 * t378 + (-Ifges(3,5) * t191 + 0.2e1 * Ifges(3,6) * t346) * qJDD(2) + (-mrSges(4,1) * t133 + mrSges(4,3) * t146 + Ifges(4,4) * t87 + Ifges(4,2) * t88 + Ifges(4,6) * t183) * t149 - t53 * t59 / 0.2e1 + t52 * t60 / 0.2e1 + t53 * t31 / 0.2e1 + (-t2 * t297 + t221 * t39 - t222 * t40 - t296 * t3) * mrSges(6,3) + (-t56 * t280 + t55 * t279 + m(6) * (qJD(5) * t233 + t188 * t2 + t193 * t3) + t188 * t10 + t193 * t9) * (-pkin(4) * t229 - pkin(6) * t108 + t121) - (-Ifges(5,2) * t38 - Ifges(5,4) * t37 + t71 * mrSges(5,1) - Ifges(5,6) * t177 + Ifges(6,3) * t356 + Ifges(6,6) * t359 + Ifges(6,5) * t360 - t72 * mrSges(5,3) + t275 / 0.2e1 + t366) * t229 - (mrSges(4,2) * t133 - mrSges(4,3) * t145 + Ifges(4,1) * t87 + Ifges(4,4) * t88 + Ifges(4,5) * t183) * t226; -t118 * t269 + t236 * t254 + (m(6) * t69 + t8) * (-pkin(4) - t136) + t198 + (t117 * t283 + m(4) * (t145 * t195 + t146 * t190)) * pkin(2) + m(5) * (t136 * t73 + t137 * t72) + (t255 + (t223 + t375 / 0.2e1) * qJD(1)) * qJD(1) - g(1) * (t197 * t204 + t246) - g(2) * (t192 * t204 + t247) + t141 * t285 / 0.2e1 - t140 * t286 / 0.2e1 + t105 * t273 + t364 * (pkin(6) + t137) + t385 * (-m(5) * t155 + t244 + t276 + t386) - t365 * (t54 - t273) + (-m(6) * (-t182 + t210) + m(5) * t256 + m(4) * t182 + t161 + t368) * g(3) + (m(5) * t115 + m(6) * t233 + t370) * (t172 * t281 + (qJD(3) * t227 - t190 * t282) * pkin(2)) - t376 * (-t273 - t340) + t136 * t26 + t137 * t27 + Ifges(3,6) * t152 + Ifges(3,5) * t153 + (mrSges(4,1) * t183 - mrSges(4,3) * t87) * t341 + (-mrSges(4,2) * t183 + mrSges(4,3) * t88) * t343 + (-m(5) * t114 + m(6) * t111 + t377) * (t172 * t282 + (qJD(3) * t228 + t190 * t281) * pkin(2)) + Ifges(3,3) * qJDD(2); -t117 * t270 + t198 - g(1) * (t197 * t202 + t246) - g(2) * (t192 * t202 + t247) + t118 * t272 - t132 * t89 + t171 * t8 + t26 * t336 + t27 * t337 + t66 * t340 - t44 * t55 - t43 * t56 - t377 * t131 + (-t111 * t131 + t171 * t69 - t39 * t43 - t40 * t44) * m(6) + (t113 * t340 + t114 * t131 - t115 * t132) * m(5) + (m(5) * t338 - m(6) * t210 + t368) * g(3) + t364 * (pkin(6) + t337) + (m(5) * t339 + t386) * t385 + (t377 * t282 + (t111 * t189 + t194 * t233) * qJD(4) * m(6) + (t189 * t72 + t194 * t73 + (-t114 * t189 + t115 * t194) * qJD(4)) * m(5) + t370 * t281) * pkin(3); (t197 * t201 - t246) * g(1) + (t192 * t201 - t247) * g(2) + (t169 * t225 - t157 - t164) * g(3) + (t373 * qJD(5) + (-t329 + t330 + t332 + t374) * m(6) + t372) * pkin(6) - t377 * t115 + (-t325 / 0.2e1 - t328 / 0.2e1 - t326 / 0.2e1 + t389) * t230 - t69 * t361 - m(6) * (t111 * t115 + t39 * t47 + t40 * t48) + (-t218 / 0.2e1 - t220 / 0.2e1 - t219 / 0.2e1 + t234 * mrSges(6,3) + t390) * t251 + (t203 + t332) * mrSges(6,3) - pkin(4) * t8 + t199 - t114 * t89 - t48 * t55 - t47 * t56; -t111 * (mrSges(6,1) * t84 + mrSges(6,2) * t83) + (Ifges(6,1) * t83 - t327) * t352 + t32 * t351 + (Ifges(6,5) * t83 - Ifges(6,6) * t84) * t350 - t39 * t55 + t40 * t56 - g(1) * (mrSges(6,1) * t127 - mrSges(6,2) * t128) - g(2) * (-mrSges(6,1) * t125 + mrSges(6,2) * t126) - t241 * t332 + (t39 * t83 + t40 * t84) * mrSges(6,3) + t275 + (-Ifges(6,2) * t84 + t33 + t82) * t353 + t366;];
tau  = t1;
