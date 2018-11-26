% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:15:02
% EndTime: 2018-11-23 15:15:15
% DurationCPUTime: 13.35s
% Computational Cost: add. (9450->607), mult. (24183->848), div. (0->0), fcn. (18338->12), ass. (0->299)
t222 = sin(qJ(2));
t217 = sin(pkin(6));
t287 = qJD(1) * t217;
t267 = t222 * t287;
t221 = sin(qJ(3));
t282 = qJD(3) * t221;
t395 = pkin(3) * t282 - t267;
t321 = -qJ(4) - pkin(8);
t260 = qJD(3) * t321;
t225 = cos(qJ(3));
t280 = qJD(4) * t225;
t176 = t221 * t260 + t280;
t177 = -qJD(4) * t221 + t225 * t260;
t216 = sin(pkin(12));
t297 = cos(pkin(12));
t234 = -t216 * t221 + t225 * t297;
t226 = cos(qJ(2));
t266 = t226 * t287;
t380 = t176 * t297 + t216 * t177 - t234 * t266;
t191 = t216 * t225 + t221 * t297;
t181 = t191 * qJD(3);
t182 = t234 * qJD(3);
t394 = pkin(4) * t181 - pkin(9) * t182 + t395;
t214 = -pkin(3) * t225 - pkin(2);
t138 = -pkin(4) * t234 - pkin(9) * t191 + t214;
t202 = t321 * t221;
t203 = t321 * t225;
t155 = t216 * t202 - t203 * t297;
t220 = sin(qJ(5));
t224 = cos(qJ(5));
t278 = qJD(5) * t224;
t279 = qJD(5) * t220;
t383 = t138 * t278 - t155 * t279 + t394 * t220 + t380 * t224;
t393 = -t380 * t220 + t394 * t224;
t144 = t224 * t155;
t325 = pkin(10) * t224;
t392 = -t182 * t325 + pkin(5) * t181 + (-t144 + (pkin(10) * t191 - t138) * t220) * qJD(5) + t393;
t236 = t182 * t220 + t191 * t278;
t391 = pkin(10) * t236 - t383;
t257 = qJD(2) * t297;
t285 = qJD(2) * t221;
t179 = -t216 * t285 + t225 * t257;
t331 = -t179 / 0.2e1;
t390 = Ifges(5,2) * t331;
t326 = pkin(3) * t216;
t212 = pkin(9) + t326;
t322 = pkin(10) + t212;
t259 = qJD(5) * t322;
t293 = t179 * t220;
t283 = qJD(2) * t225;
t180 = -t216 * t283 - t221 * t257;
t275 = pkin(3) * t285;
t119 = -pkin(4) * t180 - pkin(9) * t179 + t275;
t197 = qJD(2) * pkin(8) + t267;
t254 = qJ(4) * qJD(2) + t197;
t218 = cos(pkin(6));
t286 = qJD(1) * t221;
t265 = t218 * t286;
t150 = t225 * t254 + t265;
t139 = t216 * t150;
t288 = t218 * t225;
t209 = qJD(1) * t288;
t149 = -t221 * t254 + t209;
t85 = t149 * t297 - t139;
t52 = t220 * t119 + t224 * t85;
t389 = -pkin(10) * t293 + t220 * t259 + t52;
t292 = t179 * t224;
t51 = t224 * t119 - t220 * t85;
t388 = pkin(5) * t180 + pkin(10) * t292 - t224 * t259 - t51;
t387 = t279 - t293;
t223 = cos(qJ(6));
t173 = qJD(5) - t179;
t158 = qJD(3) * t220 - t180 * t224;
t143 = qJD(3) * pkin(3) + t149;
t258 = t297 * t150;
t81 = t216 * t143 + t258;
t76 = qJD(3) * pkin(9) + t81;
t171 = qJD(2) * t214 + qJD(4) - t266;
t98 = -pkin(4) * t179 + pkin(9) * t180 + t171;
t44 = -t220 * t76 + t224 * t98;
t35 = -pkin(10) * t158 + t44;
t27 = pkin(5) * t173 + t35;
t219 = sin(qJ(6));
t157 = qJD(3) * t224 + t180 * t220;
t45 = t220 * t98 + t224 * t76;
t36 = pkin(10) * t157 + t45;
t300 = t219 * t36;
t12 = t223 * t27 - t300;
t299 = t223 * t36;
t13 = t219 * t27 + t299;
t169 = qJD(2) * t181;
t165 = Ifges(7,3) * t169;
t255 = t223 * t157 - t158 * t219;
t94 = t157 * t219 + t158 * t223;
t327 = Ifges(7,4) * t94;
t168 = qJD(6) + t173;
t335 = -t168 / 0.2e1;
t344 = -t94 / 0.2e1;
t80 = t143 * t297 - t139;
t75 = -qJD(3) * pkin(4) - t80;
t60 = -t157 * pkin(5) + t75;
t386 = t165 + (Ifges(7,5) * t255 - Ifges(7,6) * t94) * t335 + (t12 * t255 + t13 * t94) * mrSges(7,3) - t60 * (mrSges(7,1) * t94 + mrSges(7,2) * t255) + (Ifges(7,1) * t255 - t327) * t344;
t77 = t224 * t138 - t155 * t220;
t59 = -pkin(5) * t234 - t191 * t325 + t77;
t291 = t191 * t220;
t78 = t220 * t138 + t144;
t63 = -pkin(10) * t291 + t78;
t23 = -t219 * t63 + t223 * t59;
t385 = qJD(6) * t23 + t392 * t219 - t391 * t223;
t24 = t219 * t59 + t223 * t63;
t384 = -qJD(6) * t24 + t391 * t219 + t392 * t223;
t382 = -qJD(5) * t78 + t393;
t116 = t176 * t216 - t297 * t177;
t152 = t191 * t266;
t381 = t116 - t152;
t170 = qJD(2) * t182;
t104 = qJD(5) * t157 + t170 * t224;
t289 = t217 * t226;
t263 = qJD(2) * t289;
t253 = t225 * t263;
t132 = qJD(1) * t253 + qJD(3) * t209 - t197 * t282;
t103 = (-qJ(4) * t282 + t280) * qJD(2) + t132;
t160 = t197 * t225 + t265;
t281 = qJD(3) * t225;
t364 = -t160 * qJD(3) + (-qJ(4) * t281 + (-qJD(4) - t266) * t221) * qJD(2);
t57 = t103 * t297 + t216 * t364;
t277 = qJD(2) * qJD(3);
t262 = t221 * t277;
t284 = qJD(2) * t222;
t264 = t217 * t284;
t178 = pkin(3) * t262 + qJD(1) * t264;
t90 = pkin(4) * t169 - pkin(9) * t170 + t178;
t18 = -qJD(5) * t45 - t220 * t57 + t224 * t90;
t8 = pkin(5) * t169 - pkin(10) * t104 + t18;
t105 = -qJD(5) * t158 - t170 * t220;
t17 = t220 * t90 + t224 * t57 + t98 * t278 - t279 * t76;
t9 = pkin(10) * t105 + t17;
t2 = qJD(6) * t12 + t219 * t8 + t223 * t9;
t3 = -qJD(6) * t13 - t219 * t9 + t223 * t8;
t30 = qJD(6) * t255 + t104 * t223 + t105 * t219;
t31 = -qJD(6) * t94 - t104 * t219 + t105 * t223;
t379 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t30 + Ifges(7,6) * t31;
t88 = Ifges(7,4) * t255;
t378 = -Ifges(7,2) * t94 + t88;
t376 = -Ifges(5,6) / 0.2e1;
t352 = t30 / 0.2e1;
t351 = t31 / 0.2e1;
t333 = t169 / 0.2e1;
t188 = t322 * t220;
t189 = t322 * t224;
t136 = -t188 * t219 + t189 * t223;
t375 = -qJD(6) * t136 + t219 * t389 + t223 * t388;
t135 = -t188 * t223 - t189 * t219;
t374 = qJD(6) * t135 + t219 * t388 - t223 * t389;
t369 = t171 * mrSges(5,2);
t83 = t149 * t216 + t258;
t368 = pkin(5) * t387 - t83;
t302 = t180 * mrSges(5,3);
t298 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t157 + mrSges(6,2) * t158 - t302;
t367 = Ifges(5,5) * qJD(3);
t239 = t219 * t220 - t223 * t224;
t129 = t239 * t191;
t114 = t239 * t179;
t362 = qJD(5) + qJD(6);
t145 = t362 * t239;
t366 = t114 - t145;
t194 = t219 * t224 + t220 * t223;
t113 = t194 * t179;
t146 = t362 * t194;
t365 = t113 - t146;
t363 = t17 * t224 - t18 * t220;
t360 = qJD(3) * t376;
t318 = Ifges(5,4) * t180;
t359 = t390 + t318 / 0.2e1;
t358 = t18 * mrSges(6,1) - t17 * mrSges(6,2) + Ifges(6,5) * t104 + Ifges(6,6) * t105 + t379;
t251 = mrSges(6,1) * t220 + mrSges(6,2) * t224;
t237 = t75 * t251;
t246 = Ifges(6,5) * t224 - Ifges(6,6) * t220;
t316 = Ifges(6,4) * t224;
t248 = -Ifges(6,2) * t220 + t316;
t317 = Ifges(6,4) * t220;
t250 = Ifges(6,1) * t224 - t317;
t328 = t224 / 0.2e1;
t329 = -t220 / 0.2e1;
t336 = t158 / 0.2e1;
t311 = t158 * Ifges(6,4);
t71 = t157 * Ifges(6,2) + t173 * Ifges(6,6) + t311;
t156 = Ifges(6,4) * t157;
t72 = t158 * Ifges(6,1) + t173 * Ifges(6,5) + t156;
t357 = t173 * t246 / 0.2e1 + t250 * t336 + t157 * t248 / 0.2e1 + t237 + t72 * t328 + t71 * t329;
t356 = t171 * mrSges(5,1) + t44 * mrSges(6,1) + t12 * mrSges(7,1) - t45 * mrSges(6,2) - t13 * mrSges(7,2) + t359 + t360;
t227 = qJD(2) ^ 2;
t354 = Ifges(7,4) * t352 + Ifges(7,2) * t351 + Ifges(7,6) * t333;
t353 = Ifges(7,1) * t352 + Ifges(7,4) * t351 + Ifges(7,5) * t333;
t39 = Ifges(7,2) * t255 + Ifges(7,6) * t168 + t327;
t350 = -t39 / 0.2e1;
t349 = t39 / 0.2e1;
t40 = Ifges(7,1) * t94 + Ifges(7,5) * t168 + t88;
t348 = -t40 / 0.2e1;
t347 = t40 / 0.2e1;
t346 = -t255 / 0.2e1;
t345 = t255 / 0.2e1;
t343 = t94 / 0.2e1;
t341 = t104 / 0.2e1;
t340 = t105 / 0.2e1;
t338 = -t157 / 0.2e1;
t337 = -t158 / 0.2e1;
t334 = t168 / 0.2e1;
t332 = -t173 / 0.2e1;
t330 = t180 / 0.2e1;
t324 = t255 * Ifges(7,6);
t323 = t94 * Ifges(7,5);
t320 = mrSges(5,3) * t169;
t319 = Ifges(4,4) * t221;
t172 = Ifges(5,4) * t179;
t290 = t217 * t222;
t183 = -t221 * t290 + t288;
t184 = t218 * t221 + t225 * t290;
t124 = -t183 * t297 + t184 * t216;
t56 = t103 * t216 - t297 * t364;
t314 = t124 * t56;
t154 = -t297 * t202 - t203 * t216;
t313 = t154 * t56;
t312 = t157 * Ifges(6,6);
t310 = t158 * Ifges(6,5);
t309 = t168 * Ifges(7,3);
t307 = t170 * mrSges(5,3);
t306 = t173 * Ifges(6,3);
t305 = t179 * mrSges(5,3);
t301 = t180 * Ifges(5,1);
t296 = Ifges(4,5) * qJD(3);
t295 = Ifges(4,6) * qJD(3);
t50 = -mrSges(7,1) * t255 + mrSges(7,2) * t94;
t276 = -t50 - t298;
t272 = mrSges(4,3) * t285;
t271 = mrSges(4,3) * t283;
t268 = t297 * pkin(3);
t118 = t169 * mrSges(5,1) + t170 * mrSges(5,2);
t213 = -t268 - pkin(4);
t252 = mrSges(6,1) * t224 - mrSges(6,2) * t220;
t249 = Ifges(6,1) * t220 + t316;
t247 = Ifges(6,2) * t224 + t317;
t245 = Ifges(6,5) * t220 + Ifges(6,6) * t224;
t244 = -t17 * t220 - t18 * t224;
t243 = -t45 * t220 - t44 * t224;
t242 = t44 * t220 - t45 * t224;
t125 = t216 * t183 + t184 * t297;
t238 = -t125 * t224 + t220 * t289;
t99 = -t125 * t220 - t224 * t289;
t55 = t219 * t99 - t223 * t238;
t54 = t219 * t238 + t223 * t99;
t106 = -mrSges(6,2) * t173 + mrSges(6,3) * t157;
t107 = mrSges(6,1) * t173 - mrSges(6,3) * t158;
t241 = t106 * t224 - t107 * t220;
t133 = -t197 * t281 + (-qJD(3) * t218 - t263) * t286;
t240 = t132 * t225 - t133 * t221;
t235 = (mrSges(4,1) * t221 + mrSges(4,2) * t225) * qJD(2);
t215 = Ifges(4,4) * t283;
t201 = -qJD(3) * mrSges(4,2) + t271;
t200 = qJD(3) * mrSges(4,1) - t272;
t199 = -t224 * pkin(5) + t213;
t198 = -qJD(2) * pkin(2) - t266;
t187 = qJD(3) * t235;
t186 = Ifges(4,1) * t285 + t215 + t296;
t185 = t295 + (t225 * Ifges(4,2) + t319) * qJD(2);
t166 = Ifges(6,3) * t169;
t161 = -qJD(3) * mrSges(5,2) + t305;
t159 = -t197 * t221 + t209;
t148 = -qJD(3) * t184 - t221 * t263;
t147 = qJD(3) * t183 + t253;
t134 = -mrSges(5,1) * t179 - mrSges(5,2) * t180;
t128 = t194 * t191;
t127 = t172 - t301 + t367;
t115 = pkin(5) * t291 + t154;
t84 = t147 * t297 + t216 * t148;
t82 = t147 * t216 - t148 * t297;
t74 = -mrSges(6,2) * t169 + mrSges(6,3) * t105;
t73 = mrSges(6,1) * t169 - mrSges(6,3) * t104;
t70 = t306 + t310 + t312;
t67 = mrSges(7,1) * t168 - mrSges(7,3) * t94;
t66 = -mrSges(7,2) * t168 + mrSges(7,3) * t255;
t65 = pkin(5) * t236 + t116;
t58 = -mrSges(6,1) * t105 + mrSges(6,2) * t104;
t49 = t104 * Ifges(6,1) + t105 * Ifges(6,4) + t169 * Ifges(6,5);
t48 = t104 * Ifges(6,4) + t105 * Ifges(6,2) + t169 * Ifges(6,6);
t47 = t129 * t362 - t194 * t182;
t46 = -t146 * t191 - t182 * t239;
t42 = qJD(5) * t238 - t220 * t84 + t224 * t264;
t41 = qJD(5) * t99 + t220 * t264 + t224 * t84;
t38 = t309 + t323 + t324;
t34 = -pkin(5) * t105 + t56;
t26 = -mrSges(7,2) * t169 + mrSges(7,3) * t31;
t25 = mrSges(7,1) * t169 - mrSges(7,3) * t30;
t16 = t223 * t35 - t300;
t15 = -t219 * t35 - t299;
t14 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t7 = -qJD(6) * t55 - t219 * t41 + t223 * t42;
t6 = qJD(6) * t54 + t219 * t42 + t223 * t41;
t1 = [-t125 * t320 - t238 * t74 + t41 * t106 + t42 * t107 + t147 * t201 + t148 * t200 + t84 * t161 + t54 * t25 + t55 * t26 + t6 * t66 + t7 * t67 + t99 * t73 + (-t183 * t225 - t184 * t221) * mrSges(4,3) * t277 - t276 * t82 + (t14 + t58 + t307) * t124 + ((-mrSges(3,2) * t227 - t118 - t187) * t226 + (-mrSges(3,1) * t227 + (t134 + qJD(2) * (-mrSges(4,1) * t225 + mrSges(4,2) * t221)) * qJD(2)) * t222) * t217 + m(6) * (-t17 * t238 + t18 * t99 + t41 * t45 + t42 * t44 + t75 * t82 + t314) + m(7) * (t12 * t7 + t124 * t34 + t13 * t6 + t2 * t55 + t3 * t54 + t60 * t82) + m(5) * (t314 + t125 * t57 - t80 * t82 + t81 * t84 + (t171 * t284 - t178 * t226) * t217) + m(4) * (t132 * t184 + t133 * t183 + t147 * t160 + t148 * t159 + (t198 - t266) * t264); (Ifges(7,5) * t46 + Ifges(7,6) * t47) * t334 + (-Ifges(7,1) * t129 - Ifges(7,4) * t128) * t352 + t34 * (mrSges(7,1) * t128 - mrSges(7,2) * t129) + (-Ifges(7,5) * t129 - Ifges(7,6) * t128) * t333 + (-t12 * t46 - t128 * t2 + t129 * t3 + t13 * t47) * mrSges(7,3) + (-Ifges(7,4) * t129 - Ifges(7,2) * t128) * t351 + (t312 / 0.2e1 + t310 / 0.2e1 + t306 / 0.2e1 + t38 / 0.2e1 + t70 / 0.2e1 + t324 / 0.2e1 + t323 / 0.2e1 + t309 / 0.2e1 + t356 + t359) * t181 + (Ifges(5,5) * t182 / 0.2e1 + t181 * t376 + (t198 * mrSges(4,2) + t186 / 0.2e1 - t159 * mrSges(4,3) - pkin(8) * t200 + t296 / 0.2e1) * t225 + (t198 * mrSges(4,1) - t185 / 0.2e1 - t160 * mrSges(4,3) + pkin(3) * t134 - pkin(8) * t201 - t295 / 0.2e1 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t283) * t221) * qJD(3) + (t154 * t170 - t155 * t169 - t181 * t81 - t182 * t80) * mrSges(5,3) + t380 * t161 + (t250 * t341 + t248 * t340 + t246 * t333 + Ifges(5,1) * t170 + t178 * mrSges(5,2) - Ifges(5,4) * t169 + t48 * t329 + t49 * t328 + (mrSges(5,3) + t251) * t56 + t244 * mrSges(6,3) + (t247 * t338 + t249 * t337 + t75 * t252 + t245 * t332 + t72 * t329 - t224 * t71 / 0.2e1 + t242 * mrSges(6,3)) * qJD(5)) * t191 + t298 * t116 - (-Ifges(5,4) * t170 + t178 * mrSges(5,1) + t165 / 0.2e1 + t166 / 0.2e1 - t57 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) + Ifges(7,3) / 0.2e1) * t169 + t358) * t234 + t276 * t152 + (-(pkin(2) * t284 + t198 * t222 + (-t159 * t221 + t160 * t225) * t226) * t287 + ((-t159 * t225 - t160 * t221) * qJD(3) + t240) * pkin(8)) * m(4) + t240 * mrSges(4,3) + (t369 - t301 / 0.2e1 + t172 / 0.2e1 + t127 / 0.2e1 + t243 * mrSges(6,3) + t357) * t182 + (t155 * t57 + t171 * t395 + t178 * t214 + t380 * t81 - t381 * t80 + t313) * m(5) + t214 * t118 - pkin(2) * t187 + (Ifges(7,1) * t46 + Ifges(7,4) * t47) * t343 + (Ifges(7,4) * t46 + Ifges(7,2) * t47) * t345 + t382 * t107 + t383 * t106 + (t17 * t78 + t18 * t77 + t381 * t75 + t382 * t44 + t383 * t45 + t313) * m(6) + t384 * t67 + t385 * t66 + (t115 * t34 + t2 * t24 + t23 * t3 + (-t152 + t65) * t60 + t385 * t13 + t384 * t12) * m(7) + (0.3e1 / 0.2e1 * t225 ^ 2 - 0.3e1 / 0.2e1 * t221 ^ 2) * Ifges(4,4) * t277 + t46 * t347 + t47 * t349 - t129 * t353 - t128 * t354 + t23 * t25 + t24 * t26 + t60 * (-mrSges(7,1) * t47 + mrSges(7,2) * t46) + t65 * t50 + t77 * t73 + t78 * t74 + t115 * t14 + t154 * t58 + ((t221 * t200 - t225 * t201) * t226 - t222 * t134) * t287; t247 * t340 + t48 * t328 + t374 * t66 + (t12 * t375 + t13 * t374 + t135 * t3 + t136 * t2 + t199 * t34 + t368 * t60) * m(7) + t375 * t67 + (m(6) * t213 - mrSges(5,1) - t252) * t56 + (-Ifges(7,4) * t114 - Ifges(7,2) * t113) * t346 + (-Ifges(7,5) * t114 - Ifges(7,6) * t113) * t335 + (-Ifges(7,1) * t114 - Ifges(7,4) * t113) * t344 + (m(6) * t363 - t106 * t279 - t107 * t278 - t220 * t73 + t224 * t74) * t212 + (-mrSges(7,1) * t365 + mrSges(7,2) * t366) * t60 + (m(6) * t212 * t243 + t357) * qJD(5) + (t271 - t201) * t159 - t320 * t326 - t268 * t307 - t81 * t302 - t72 * t292 / 0.2e1 + t71 * t293 / 0.2e1 + (Ifges(7,5) * t194 - Ifges(7,6) * t239 + t245) * t333 + t34 * (mrSges(7,1) * t239 + mrSges(7,2) * t194) + (Ifges(7,4) * t194 - Ifges(7,2) * t239) * t351 + (Ifges(7,1) * t194 - Ifges(7,4) * t239) * t352 + (-t12 * t366 + t13 * t365 - t194 * t3 - t2 * t239) * mrSges(7,3) - t239 * t354 + (-Ifges(7,5) * t145 - Ifges(7,6) * t146) * t334 + (-Ifges(7,1) * t145 - Ifges(7,4) * t146) * t343 + (-Ifges(7,4) * t145 - Ifges(7,2) * t146) * t345 + t185 * t285 / 0.2e1 - t134 * t275 + (-m(6) * t75 - t298) * t83 + t368 * t50 + (-t171 * t275 + t80 * t83 - t81 * t85 + (t216 * t57 - t297 * t56) * pkin(3)) * m(5) - t221 * t227 * (Ifges(4,1) * t225 - t319) / 0.2e1 + (t272 + t200) * t160 + t80 * t305 - (-Ifges(4,2) * t285 + t186 + t215) * t283 / 0.2e1 + t277 * Ifges(4,5) * t225 / 0.2e1 - Ifges(4,6) * t262 / 0.2e1 + (t127 + t172) * t331 - m(6) * (t44 * t51 + t45 * t52) + (t318 + t70 + t38) * t330 + (t250 * t337 + t248 * t338 + Ifges(5,1) * t330 + t246 * t332 - t237 - t367 / 0.2e1 - t369) * t179 + (-Ifges(6,5) * t337 - Ifges(7,5) * t344 - Ifges(6,6) * t338 - Ifges(7,6) * t346 - Ifges(6,3) * t332 - Ifges(7,3) * t335 + t356 + t360 + t390) * t180 + (-t387 * t45 + (-t278 + t292) * t44 + t363) * mrSges(6,3) + t220 * t49 / 0.2e1 + t213 * t58 + t199 * t14 + t249 * t341 - t145 * t347 - t114 * t348 - t146 * t349 - t113 * t350 + t194 * t353 - t198 * t235 - t57 * mrSges(5,2) - t52 * t106 - t51 * t107 - t132 * mrSges(4,2) + t133 * mrSges(4,1) + t135 * t25 + t136 * t26 - t85 * t161 - Ifges(5,6) * t169 + Ifges(5,5) * t170; -t239 * t25 + t194 * t26 + t220 * t74 + t224 * t73 + t365 * t67 + t366 * t66 + t241 * qJD(5) - t276 * t180 + (-t161 - t241) * t179 + t118 + (t12 * t365 + t13 * t366 + t180 * t60 + t194 * t2 - t239 * t3) * m(7) + (-t173 * t242 + t180 * t75 - t244) * m(6) + (-t179 * t81 - t180 * t80 + t178) * m(5); t71 * t336 + (Ifges(6,1) * t157 - t311) * t337 + (Ifges(6,5) * t157 - Ifges(6,6) * t158) * t332 + (-t158 * t50 + t219 * t26 + t223 * t25 + (-t219 * t67 + t223 * t66) * qJD(6) + (-t158 * t60 + t2 * t219 + t223 * t3 + (-t12 * t219 + t13 * t223) * qJD(6)) * m(7)) * pkin(5) - t94 * t350 - m(7) * (t12 * t15 + t13 * t16) + (-Ifges(6,2) * t158 + t156 + t72) * t338 + t255 * t348 + (t157 * t44 + t158 * t45) * mrSges(6,3) + t358 + t166 + t378 * t346 - t16 * t66 - t15 * t67 - t44 * t106 + t45 * t107 - t75 * (mrSges(6,1) * t158 + mrSges(6,2) * t157) + t386; t39 * t343 - t12 * t66 + t13 * t67 + (t378 + t40) * t346 + t379 + t386;];
tauc  = t1(:);
