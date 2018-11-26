% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2018-11-23 17:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:07:49
% EndTime: 2018-11-23 17:08:01
% DurationCPUTime: 11.66s
% Computational Cost: add. (10693->666), mult. (24689->913), div. (0->0), fcn. (16053->8), ass. (0->309)
t284 = sin(qJ(2));
t274 = t284 * qJD(1);
t269 = pkin(2) * t274;
t287 = cos(qJ(2));
t308 = pkin(8) * t284 - qJ(3) * t287;
t201 = qJD(1) * t308 + t269;
t343 = qJD(1) * t287;
t270 = pkin(7) * t343;
t237 = pkin(3) * t343 + t270;
t283 = sin(qJ(4));
t286 = cos(qJ(4));
t146 = -t201 * t283 + t286 * t237;
t335 = t286 * qJD(5);
t338 = qJD(4) * t283;
t288 = -pkin(2) - pkin(8);
t347 = qJ(5) - t288;
t458 = -(-qJ(5) * t283 * t284 + pkin(4) * t287) * qJD(1) - t146 + t338 * t347 - t335;
t147 = t286 * t201 + t283 * t237;
t242 = t347 * t286;
t326 = t286 * t274;
t457 = qJ(5) * t326 + qJD(4) * t242 + t283 * qJD(5) + t147;
t281 = cos(pkin(10));
t280 = sin(pkin(10));
t350 = t280 * t283;
t183 = -t274 * t350 + t281 * t326;
t337 = qJD(4) * t286;
t208 = t280 * t338 - t281 * t337;
t346 = -t208 + t183;
t300 = t280 * t286 + t281 * t283;
t295 = t300 * t284;
t184 = qJD(1) * t295;
t209 = t300 * qJD(4);
t345 = -t209 - t184;
t282 = sin(qJ(6));
t285 = cos(qJ(6));
t261 = t274 + qJD(4);
t230 = -qJD(2) * t283 - t286 * t343;
t341 = qJD(2) * t286;
t231 = -t283 * t343 + t341;
t158 = t230 * t280 + t231 * t281;
t443 = pkin(9) * t158;
t322 = -qJ(3) * t284 - pkin(1);
t225 = t287 * t288 + t322;
t188 = t225 * qJD(1);
t268 = pkin(7) * t274;
t236 = -pkin(3) * t274 - t268;
t412 = qJD(3) - t236;
t193 = qJD(2) * t288 + t412;
t123 = -t188 * t283 + t286 * t193;
t104 = -qJ(5) * t231 + t123;
t94 = pkin(4) * t261 + t104;
t124 = t188 * t286 + t193 * t283;
t105 = qJ(5) * t230 + t124;
t98 = t280 * t105;
t47 = t281 * t94 - t98;
t34 = pkin(5) * t261 - t443 + t47;
t318 = t281 * t230 - t231 * t280;
t434 = pkin(9) * t318;
t349 = t281 * t105;
t48 = t280 * t94 + t349;
t35 = t48 + t434;
t10 = t282 * t34 + t285 * t35;
t334 = qJD(1) * qJD(2);
t323 = t287 * t334;
t257 = Ifges(7,3) * t323;
t89 = t158 * t285 + t282 * t318;
t374 = Ifges(7,4) * t89;
t253 = qJD(6) + t261;
t380 = -t253 / 0.2e1;
t437 = -t158 * t282 + t285 * t318;
t79 = Ifges(7,4) * t437;
t40 = Ifges(7,1) * t89 + Ifges(7,5) * t253 + t79;
t403 = -t89 / 0.2e1;
t405 = -t437 / 0.2e1;
t9 = -t282 * t35 + t285 * t34;
t279 = qJD(2) * qJ(3);
t215 = t279 + t237;
t160 = -pkin(4) * t230 + qJD(5) + t215;
t97 = -pkin(5) * t318 + t160;
t456 = t257 + (Ifges(7,5) * t437 - Ifges(7,6) * t89) * t380 + (t10 * t89 + t437 * t9) * mrSges(7,3) + (-Ifges(7,2) * t89 + t40 + t79) * t405 - t97 * (mrSges(7,1) * t89 + mrSges(7,2) * t437) + (Ifges(7,1) * t437 - t374) * t403;
t422 = t457 * t280 + t281 * t458;
t421 = t280 * t458 - t457 * t281;
t299 = -t281 * t286 + t350;
t150 = t282 * t299 - t285 * t300;
t333 = qJD(2) * qJD(4);
t336 = qJD(4) * t287;
t342 = qJD(2) * t284;
t173 = -t283 * t333 + (t283 * t342 - t286 * t336) * qJD(1);
t325 = t283 * t336;
t294 = t284 * t341 + t325;
t174 = qJD(1) * t294 - t286 * t333;
t112 = t173 * t281 + t174 * t280;
t324 = t284 * t334;
t260 = pkin(2) * t324;
t339 = qJD(3) * t284;
t293 = qJD(2) * t308 - t339;
t167 = qJD(1) * t293 + t260;
t340 = qJD(2) * t287;
t401 = pkin(3) + pkin(7);
t239 = t401 * t340;
t221 = qJD(1) * t239;
t68 = -qJD(4) * t124 - t167 * t283 + t286 * t221;
t44 = pkin(4) * t323 - qJ(5) * t173 - qJD(5) * t231 + t68;
t67 = t286 * t167 - t188 * t338 + t193 * t337 + t283 * t221;
t46 = qJ(5) * t174 + qJD(5) * t230 + t67;
t15 = -t280 * t46 + t281 * t44;
t11 = pkin(5) * t323 - pkin(9) * t112 + t15;
t111 = -t173 * t280 + t174 * t281;
t16 = t280 * t44 + t281 * t46;
t12 = pkin(9) * t111 + t16;
t2 = qJD(6) * t9 + t11 * t282 + t12 * t285;
t3 = -qJD(6) * t10 + t11 * t285 - t12 * t282;
t292 = qJD(6) * t150 + t208 * t282 - t285 * t209;
t320 = t183 * t282 + t285 * t184;
t353 = t320 - t292;
t415 = -t282 * t300 - t285 * t299;
t118 = t183 * t285 - t184 * t282;
t80 = qJD(6) * t415 - t208 * t285 - t209 * t282;
t446 = t118 + t80;
t455 = -t10 * t446 + t150 * t2 - t415 * t3 + t353 * t9;
t39 = Ifges(7,2) * t437 + Ifges(7,6) * t253 + t374;
t453 = t39 / 0.2e1;
t431 = -qJD(1) / 0.2e1;
t448 = -pkin(5) * t343 - pkin(9) * t345 + t422;
t447 = -pkin(9) * t346 + t421;
t359 = t231 * Ifges(5,4);
t143 = t230 * Ifges(5,2) + t261 * Ifges(5,6) + t359;
t226 = Ifges(5,4) * t230;
t144 = t231 * Ifges(5,1) + t261 * Ifges(5,5) + t226;
t305 = t123 * t283 - t124 * t286;
t369 = Ifges(5,4) * t283;
t310 = Ifges(5,2) * t286 + t369;
t368 = Ifges(5,4) * t286;
t312 = Ifges(5,1) * t283 + t368;
t315 = mrSges(5,1) * t286 - mrSges(5,2) * t283;
t363 = Ifges(5,6) * t286;
t367 = Ifges(5,5) * t283;
t375 = -t286 / 0.2e1;
t376 = -t283 / 0.2e1;
t378 = -t261 / 0.2e1;
t381 = -t231 / 0.2e1;
t382 = -t230 / 0.2e1;
t445 = t305 * mrSges(5,3) + t143 * t375 + t144 * t376 + t215 * t315 + (t363 + t367) * t378 + t310 * t382 + t312 * t381;
t444 = qJD(3) + t268;
t442 = -mrSges(3,1) + mrSges(4,2);
t441 = t158 * Ifges(6,4);
t327 = -pkin(4) * t286 - pkin(3);
t414 = pkin(4) * t337 - t274 * t327 + t444;
t244 = -pkin(2) * t287 + t322;
t216 = t244 * qJD(1);
t247 = -t270 - t279;
t429 = qJD(2) / 0.2e1;
t430 = -qJD(2) / 0.2e1;
t440 = t247 * mrSges(4,1) - t216 * mrSges(4,2) - Ifges(4,5) * t430 - Ifges(3,6) * t429 - t445 + ((Ifges(3,2) + Ifges(4,3)) * t287 + (Ifges(3,4) + Ifges(4,6)) * t284) * t431;
t32 = qJD(6) * t437 + t111 * t282 + t112 * t285;
t33 = -qJD(6) * t89 + t111 * t285 - t112 * t282;
t439 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t32 + Ifges(7,6) * t33;
t438 = t15 * t299 - t16 * t300 - t345 * t47 - t346 * t48;
t77 = Ifges(6,2) * t318 + t261 * Ifges(6,6) + t441;
t435 = t77 / 0.2e1;
t377 = t261 / 0.2e1;
t390 = -t318 / 0.2e1;
t241 = t347 * t283;
t163 = t241 * t280 - t281 * t242;
t130 = pkin(9) * t299 + t163;
t164 = -t281 * t241 - t280 * t242;
t131 = -pkin(9) * t300 + t164;
t70 = t130 * t285 - t131 * t282;
t428 = qJD(6) * t70 + t282 * t448 + t285 * t447;
t71 = t130 * t282 + t131 * t285;
t427 = -qJD(6) * t71 - t282 * t447 + t285 * t448;
t426 = Ifges(6,4) * t318;
t265 = pkin(4) * t281 + pkin(5);
t372 = pkin(4) * t280;
t203 = t265 * t285 - t282 * t372;
t52 = -t104 * t280 - t349;
t36 = t52 - t434;
t53 = t281 * t104 - t98;
t37 = t53 - t443;
t420 = t203 * qJD(6) - t282 * t36 - t285 * t37;
t204 = t265 * t282 + t285 * t372;
t419 = -t204 * qJD(6) + t282 * t37 - t285 * t36;
t416 = pkin(5) * t346 + t414;
t251 = t401 * t284;
t162 = t286 * t225 + t283 * t251;
t413 = t283 * t67 + t286 * t68;
t243 = -qJD(2) * pkin(2) + t444;
t267 = Ifges(3,4) * t343;
t329 = Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t365 = Ifges(4,6) * t287;
t410 = t329 * t261 + t123 * mrSges(5,1) + t243 * mrSges(4,1) + t47 * mrSges(6,1) + t9 * mrSges(7,1) + t231 * Ifges(5,5) + t230 * Ifges(5,6) + Ifges(3,1) * t274 / 0.2e1 + Ifges(3,5) * t429 + t267 / 0.2e1 + Ifges(4,4) * t430 + (-t284 * Ifges(4,2) - t365) * t431 + t253 * Ifges(7,3) + t89 * Ifges(7,5) + t437 * Ifges(7,6) + t158 * Ifges(6,5) + t318 * Ifges(6,6) - t10 * mrSges(7,2) - t124 * mrSges(5,2) - t216 * mrSges(4,3) - t48 * mrSges(6,2) + (Ifges(5,3) + Ifges(6,3)) * t377;
t409 = t68 * mrSges(5,1) + t15 * mrSges(6,1) - t67 * mrSges(5,2) - t16 * mrSges(6,2) + Ifges(5,5) * t173 + Ifges(6,5) * t112 + Ifges(5,6) * t174 + Ifges(6,6) * t111 + t439;
t407 = t32 / 0.2e1;
t406 = t33 / 0.2e1;
t404 = t437 / 0.2e1;
t402 = t89 / 0.2e1;
t400 = pkin(1) * mrSges(3,1);
t399 = pkin(1) * mrSges(3,2);
t397 = t111 / 0.2e1;
t396 = t112 / 0.2e1;
t196 = t299 * t287;
t197 = t300 * t287;
t134 = t196 * t285 + t197 * t282;
t395 = t134 / 0.2e1;
t135 = t196 * t282 - t197 * t285;
t394 = t135 / 0.2e1;
t393 = t143 / 0.2e1;
t392 = t150 / 0.2e1;
t391 = t415 / 0.2e1;
t389 = t318 / 0.2e1;
t388 = -t158 / 0.2e1;
t387 = t158 / 0.2e1;
t386 = t196 / 0.2e1;
t385 = -t197 / 0.2e1;
t384 = -t300 / 0.2e1;
t383 = -t299 / 0.2e1;
t379 = t253 / 0.2e1;
t373 = pkin(4) * t231;
t273 = pkin(2) * t342;
t178 = t273 + t293;
t220 = t286 * t239;
t321 = qJ(5) * t287 - t225;
t62 = pkin(4) * t340 + t220 + t321 * t337 + (-qJ(5) * t342 - qJD(4) * t251 + qJD(5) * t287 - t178) * t283;
t84 = t286 * t178 - t225 * t338 + t283 * t239 + t251 * t337;
t69 = qJ(5) * t294 - t287 * t335 + t84;
t26 = t280 * t62 + t281 * t69;
t366 = Ifges(5,5) * t286;
t364 = Ifges(5,6) * t283;
t352 = qJD(2) * mrSges(3,2);
t348 = t286 * t287;
t266 = t283 * pkin(4) + qJ(3);
t233 = t286 * t251;
t138 = pkin(4) * t284 + t283 * t321 + t233;
t145 = -qJ(5) * t348 + t162;
t73 = t280 * t138 + t281 * t145;
t165 = -mrSges(5,1) * t230 + mrSges(5,2) * t231;
t249 = -mrSges(4,1) * t343 - qJD(2) * mrSges(4,3);
t344 = -t249 + t165;
t252 = t401 * t287;
t332 = -Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t331 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t330 = -0.3e1 / 0.2e1 * Ifges(4,6) - 0.3e1 / 0.2e1 * Ifges(3,4);
t328 = m(4) * pkin(7) + mrSges(4,1);
t210 = pkin(4) * t348 + t252;
t8 = -t33 * mrSges(7,1) + t32 * mrSges(7,2);
t25 = -t280 * t69 + t281 * t62;
t57 = -t111 * mrSges(6,1) + t112 * mrSges(6,2);
t72 = t281 * t138 - t145 * t280;
t238 = t401 * t342;
t317 = m(4) * t243 + (mrSges(4,1) + mrSges(3,3)) * t274 + t442 * qJD(2);
t316 = m(4) * t247 - mrSges(3,3) * t343 + t249 + t352;
t314 = mrSges(5,1) * t283 + mrSges(5,2) * t286;
t313 = Ifges(5,1) * t286 - t369;
t311 = -Ifges(5,2) * t283 + t368;
t58 = pkin(5) * t284 + pkin(9) * t197 + t72;
t59 = pkin(9) * t196 + t73;
t23 = -t282 * t59 + t285 * t58;
t24 = t282 * t58 + t285 * t59;
t153 = mrSges(5,1) * t323 - mrSges(5,3) * t173;
t154 = -mrSges(5,2) * t323 + mrSges(5,3) * t174;
t304 = t286 * t153 + t283 * t154;
t176 = -mrSges(5,2) * t261 + mrSges(5,3) * t230;
t177 = mrSges(5,1) * t261 - mrSges(5,3) * t231;
t303 = t286 * t176 - t283 * t177;
t297 = -qJ(3) * t340 - t339;
t278 = qJD(2) * qJD(3);
t198 = -qJD(1) * t238 + t278;
t136 = -pkin(4) * t174 + t198;
t166 = -pkin(4) * t325 + (-pkin(7) + t327) * t342;
t259 = Ifges(5,3) * t323;
t258 = Ifges(6,3) * t323;
t240 = pkin(7) * t324 - t278;
t234 = (mrSges(4,2) * t287 - mrSges(4,3) * t284) * qJD(1);
t205 = t273 + t297;
t187 = pkin(5) * t300 + t266;
t185 = qJD(1) * t297 + t260;
t161 = -t225 * t283 + t233;
t148 = -pkin(5) * t196 + t210;
t140 = qJD(2) * t295 + qJD(4) * t196;
t139 = -t299 * t342 + t300 * t336;
t133 = mrSges(6,1) * t261 - mrSges(6,3) * t158;
t132 = -mrSges(6,2) * t261 + mrSges(6,3) * t318;
t117 = pkin(5) * t158 + t373;
t115 = -mrSges(5,1) * t174 + mrSges(5,2) * t173;
t107 = t173 * Ifges(5,1) + t174 * Ifges(5,4) + Ifges(5,5) * t323;
t106 = t173 * Ifges(5,4) + t174 * Ifges(5,2) + Ifges(5,6) * t323;
t96 = mrSges(6,1) * t323 - mrSges(6,3) * t112;
t95 = -mrSges(6,2) * t323 + mrSges(6,3) * t111;
t93 = -mrSges(6,1) * t318 + mrSges(6,2) * t158;
t92 = -pkin(5) * t139 + t166;
t85 = -qJD(4) * t162 - t178 * t283 + t220;
t78 = t158 * Ifges(6,1) + t261 * Ifges(6,5) + t426;
t75 = mrSges(7,1) * t253 - mrSges(7,3) * t89;
t74 = -mrSges(7,2) * t253 + mrSges(7,3) * t437;
t66 = -pkin(5) * t111 + t136;
t56 = t112 * Ifges(6,1) + t111 * Ifges(6,4) + Ifges(6,5) * t323;
t55 = t112 * Ifges(6,4) + t111 * Ifges(6,2) + Ifges(6,6) * t323;
t50 = -qJD(6) * t135 + t139 * t285 - t140 * t282;
t49 = qJD(6) * t134 + t139 * t282 + t140 * t285;
t41 = -mrSges(7,1) * t437 + mrSges(7,2) * t89;
t28 = -mrSges(7,2) * t323 + mrSges(7,3) * t33;
t27 = mrSges(7,1) * t323 - mrSges(7,3) * t32;
t20 = pkin(9) * t139 + t26;
t19 = pkin(5) * t340 - pkin(9) * t140 + t25;
t7 = t32 * Ifges(7,1) + t33 * Ifges(7,4) + Ifges(7,5) * t323;
t6 = t32 * Ifges(7,4) + t33 * Ifges(7,2) + Ifges(7,6) * t323;
t5 = -qJD(6) * t24 + t19 * t285 - t20 * t282;
t4 = qJD(6) * t23 + t19 * t282 + t20 * t285;
t1 = [m(4) * (t185 * t244 + t205 * t216) + (t257 / 0.2e1 + t258 / 0.2e1 + t259 / 0.2e1 - t185 * mrSges(4,3) + t409) * t284 + m(5) * (t123 * t85 + t124 * t84 + t161 * t68 + t162 * t67 + t198 * t252 - t215 * t238) + t136 * (-mrSges(6,1) * t196 - mrSges(6,2) * t197) + (-Ifges(6,1) * t197 + Ifges(6,4) * t196) * t396 + (-Ifges(6,4) * t197 + Ifges(6,2) * t196) * t397 + (t139 * t48 - t140 * t47 + t15 * t197 + t16 * t196) * mrSges(6,3) + m(7) * (t10 * t4 + t148 * t66 + t2 * t24 + t23 * t3 + t5 * t9 + t92 * t97) + m(6) * (t136 * t210 + t15 * t72 + t16 * t73 + t160 * t166 + t25 * t47 + t26 * t48) + t252 * t115 + t205 * t234 - t238 * t165 + t210 * t57 + t84 * t176 + t85 * t177 + t166 * t93 + t160 * (-mrSges(6,1) * t139 + mrSges(6,2) * t140) + t161 * t153 + t162 * t154 + t140 * t78 / 0.2e1 + t148 * t8 + t26 * t132 + t25 * t133 + t66 * (-mrSges(7,1) * t134 + mrSges(7,2) * t135) + t73 * t95 + t72 * t96 + t97 * (-mrSges(7,1) * t50 + mrSges(7,2) * t49) + t92 * t41 + t4 * t74 + t5 * t75 + t49 * t40 / 0.2e1 + t24 * t28 + t23 * t27 + t50 * t453 + (-t173 * t312 / 0.2e1 - t174 * t310 / 0.2e1 + t198 * t315 + t185 * mrSges(4,2) + t107 * t376 + t106 * t375 - t328 * t240 + (t283 * t68 - t286 * t67) * mrSges(5,3) + ((t364 - t366) * t377 + t311 * t382 + t313 * t381 - t215 * t314 + t144 * t375 + t283 * t393 + (t123 * t286 + t124 * t283) * mrSges(5,3)) * qJD(4) + (t332 * qJD(2) + t317 * pkin(7) + (Ifges(7,5) * t394 + Ifges(7,6) * t395 + Ifges(6,5) * t385 + Ifges(6,6) * t386 + (-t367 / 0.2e1 - t363 / 0.2e1 - t330) * t287 - t244 * mrSges(4,3) - 0.2e1 * t399 + (Ifges(7,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + t328 * pkin(7) + t329) * t284) * qJD(1) + t410) * qJD(2)) * t287 + (Ifges(6,5) * t140 + Ifges(6,6) * t139) * t377 + (Ifges(7,5) * t49 + Ifges(7,6) * t50) * t379 + t56 * t385 + t55 * t386 + (Ifges(6,1) * t140 + Ifges(6,4) * t139) * t387 + (Ifges(6,4) * t140 + Ifges(6,2) * t139) * t389 + t7 * t394 + t6 * t395 + t139 * t435 + (Ifges(7,1) * t49 + Ifges(7,4) * t50) * t402 + (Ifges(7,4) * t49 + Ifges(7,2) * t50) * t404 + (Ifges(7,4) * t135 + Ifges(7,2) * t134) * t406 + (Ifges(7,1) * t135 + Ifges(7,4) * t134) * t407 + (t331 * qJD(2) + t316 * pkin(7) + (-t244 * mrSges(4,2) + t284 * t330 - 0.2e1 * t400) * qJD(1) + t440) * t342 + (t10 * t50 + t134 * t2 - t135 * t3 - t49 * t9) * mrSges(7,3); (t208 / 0.2e1 - t183 / 0.2e1) * t77 + t136 * (mrSges(6,1) * t300 - mrSges(6,2) * t299) + (-Ifges(6,1) * t299 - Ifges(6,4) * t300) * t396 + (-Ifges(6,4) * t299 - Ifges(6,2) * t300) * t397 + t66 * (-mrSges(7,1) * t150 + mrSges(7,2) * t415) + (Ifges(7,4) * t415 + Ifges(7,2) * t150) * t406 + (Ifges(7,1) * t415 + Ifges(7,4) * t150) * t407 + (-t209 / 0.2e1 - t184 / 0.2e1) * t78 + (-Ifges(6,5) * t209 + Ifges(6,6) * t208) * t377 + (-Ifges(6,1) * t209 + Ifges(6,4) * t208) * t387 + (-Ifges(6,4) * t209 + Ifges(6,2) * t208) * t389 + (Ifges(7,5) * t320 + Ifges(7,6) * t118) * t380 + (Ifges(7,1) * t320 + Ifges(7,4) * t118) * t403 + (Ifges(7,4) * t320 + Ifges(7,2) * t118) * t405 + (-m(4) * t216 - t234) * (-qJ(3) * t343 + t269) + (-t320 / 0.2e1 + t292 / 0.2e1) * t40 + t286 * t107 / 0.2e1 + t445 * qJD(4) + t266 * t57 - t240 * mrSges(4,3) - t236 * t165 + (mrSges(7,1) * t446 - mrSges(7,2) * t353) * t97 + t187 * t8 - t147 * t176 - t146 * t177 + t163 * t96 + t164 * t95 + t455 * mrSges(7,3) + qJ(3) * t115 + t70 * t27 + t71 * t28 + (t198 * qJ(3) - t123 * t146 - t124 * t147 + t215 * t412) * m(5) + (t304 + m(5) * t413 + (-m(5) * t305 + t303) * qJD(4)) * t288 + t414 * t93 + t416 * t41 + t174 * t311 / 0.2e1 + t173 * t313 / 0.2e1 + t198 * t314 + t421 * t132 + t422 * t133 + (t136 * t266 + t15 * t163 + t16 * t164 + t160 * t414 + t421 * t48 + t422 * t47) * m(6) + (-t118 / 0.2e1 - t80 / 0.2e1) * t39 + (Ifges(7,5) * t292 - Ifges(7,6) * t80) * t379 + (Ifges(7,1) * t292 - Ifges(7,4) * t80) * t402 + (Ifges(7,4) * t292 - Ifges(7,2) * t80) * t404 - t413 * mrSges(5,3) + t106 * t376 + (Ifges(6,5) * t184 + Ifges(6,6) * t183) * t378 + t56 * t383 + t55 * t384 + (Ifges(6,1) * t184 + Ifges(6,4) * t183) * t388 + (Ifges(6,4) * t184 + Ifges(6,2) * t183) * t390 + t7 * t391 + t6 * t392 + t427 * t75 + t428 * t74 + (t10 * t428 + t187 * t66 + t2 * t71 + t3 * t70 + t416 * t97 + t427 * t9) * m(7) + m(4) * (-qJ(3) * t240 - qJD(3) * t247) + t438 * mrSges(6,3) + (((-qJ(3) * mrSges(4,1) + t331) * qJD(2) + (t400 + (Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1) * t284) * qJD(1) + (-t316 + t352) * pkin(7) - t440) * t284 + ((Ifges(6,5) * t383 + Ifges(6,6) * t384 + Ifges(7,5) * t391 + Ifges(7,6) * t392 + t366 / 0.2e1 - t364 / 0.2e1 - pkin(2) * mrSges(4,1) + t332) * qJD(2) - t267 / 0.2e1 + ((-m(4) * pkin(2) + t442) * qJD(2) - t317) * pkin(7) + (-Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t274 + (t399 - t365 / 0.2e1) * qJD(1) - t410) * t287) * qJD(1) + t344 * qJD(3) + (mrSges(6,1) * t346 + mrSges(6,2) * t345) * t160; t415 * t27 - t150 * t28 + t300 * t95 - t299 * t96 - t353 * t75 + t446 * t74 + t345 * t133 + t346 * t132 + t303 * qJD(4) + (t234 + t303) * t274 + (t328 * t343 - t344 - t41 - t93) * qJD(2) - m(4) * (-qJD(2) * t247 - t216 * t274) + t304 + (-qJD(2) * t97 - t455) * m(7) + (-qJD(2) * t160 - t438) * m(6) + (-qJD(2) * t215 - t261 * t305 + t413) * m(5); t456 + ((t15 * t281 + t16 * t280) * pkin(4) - t160 * t373 - t47 * t52 - t48 * t53) * m(6) + t258 + t259 + (-Ifges(6,2) * t158 + t426 + t78) * t390 + (Ifges(5,5) * t230 + Ifges(6,5) * t318 - Ifges(5,6) * t231 - Ifges(6,6) * t158) * t378 + (t158 * t48 + t318 * t47) * mrSges(6,3) - t160 * (mrSges(6,1) * t158 + mrSges(6,2) * t318) + t158 * t435 + (-t231 * t93 + t280 * t95 + t281 * t96) * pkin(4) + t409 + (-Ifges(5,2) * t231 + t144 + t226) * t382 - t215 * (mrSges(5,1) * t231 + mrSges(5,2) * t230) + t89 * t453 + t203 * t27 + t204 * t28 + (t123 * t230 + t124 * t231) * mrSges(5,3) - t123 * t176 + t124 * t177 - t53 * t132 - t52 * t133 - t117 * t41 + t419 * t75 + t420 * t74 + (t10 * t420 - t117 * t97 + t2 * t204 + t203 * t3 + t419 * t9) * m(7) + (Ifges(5,1) * t230 - t359) * t381 + t231 * t393 + (Ifges(6,1) * t318 - t441) * t388; -t318 * t132 + t158 * t133 - t437 * t74 + t89 * t75 + t57 + t8 + (-t10 * t437 + t89 * t9 + t66) * m(7) + (t158 * t47 - t318 * t48 + t136) * m(6); t10 * t75 + t39 * t402 - t9 * t74 + t439 + t456;];
tauc  = t1(:);
