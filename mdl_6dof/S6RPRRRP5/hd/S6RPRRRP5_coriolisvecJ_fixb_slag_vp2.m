% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:26
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:26:32
% EndTime: 2018-11-23 16:26:39
% DurationCPUTime: 6.60s
% Computational Cost: add. (12849->523), mult. (34420->682), div. (0->0), fcn. (26706->8), ass. (0->241)
t359 = Ifges(6,1) + Ifges(7,1);
t357 = Ifges(7,4) + Ifges(6,5);
t358 = Ifges(6,6) - Ifges(7,6);
t213 = sin(qJ(5));
t216 = cos(qJ(5));
t211 = sin(pkin(10));
t215 = sin(qJ(3));
t212 = cos(pkin(10));
t218 = cos(qJ(3));
t289 = t212 * t218;
t195 = -t211 * t215 + t289;
t186 = t195 * qJD(1);
t196 = t211 * t218 + t212 * t215;
t187 = t196 * qJD(1);
t214 = sin(qJ(4));
t217 = cos(qJ(4));
t238 = t186 * t214 + t217 * t187;
t278 = qJD(3) + qJD(4);
t146 = t213 * t238 - t216 * t278;
t331 = -t146 / 0.2e1;
t188 = t195 * qJD(3);
t179 = qJD(1) * t188;
t189 = t196 * qJD(3);
t180 = qJD(1) * t189;
t120 = qJD(4) * t238 + t179 * t214 + t217 * t180;
t269 = t217 * t186 - t187 * t214;
t119 = qJD(4) * t269 + t179 * t217 - t180 * t214;
t82 = -qJD(5) * t146 + t216 * t119;
t147 = t213 * t278 + t216 * t238;
t83 = qJD(5) * t147 + t213 * t119;
t356 = (-Ifges(6,4) + Ifges(7,5)) * t83 + t359 * t82 + t357 * t120;
t143 = Ifges(6,4) * t146;
t156 = qJD(5) - t269;
t298 = Ifges(7,5) * t146;
t347 = t147 * t359 + t357 * t156 - t143 + t298;
t24 = mrSges(6,1) * t83 + mrSges(6,2) * t82;
t308 = pkin(7) + qJ(2);
t201 = t308 * t211;
t197 = qJD(1) * t201;
t202 = t308 * t212;
t198 = qJD(1) * t202;
t204 = qJD(2) * t289;
t279 = qJD(1) * qJD(2);
t282 = qJD(3) * t218;
t137 = -t197 * t282 + qJD(1) * t204 + (-qJD(3) * t198 - t211 * t279) * t215;
t127 = -pkin(8) * t180 + t137;
t166 = -t197 * t215 + t198 * t218;
t228 = t196 * qJD(2);
t138 = -qJD(1) * t228 - qJD(3) * t166;
t226 = -pkin(8) * t179 + t138;
t145 = pkin(8) * t186 + t166;
t140 = t217 * t145;
t165 = -t218 * t197 - t198 * t215;
t144 = -pkin(8) * t187 + t165;
t141 = qJD(3) * pkin(3) + t144;
t91 = t214 * t141 + t140;
t32 = qJD(4) * t91 + t127 * t214 - t217 * t226;
t355 = m(6) * t32 + t24;
t354 = t213 * t357 + t216 * t358;
t296 = Ifges(7,5) * t216;
t299 = Ifges(6,4) * t216;
t353 = t213 * t359 - t296 + t299;
t280 = qJD(5) * t216;
t281 = qJD(5) * t213;
t139 = t214 * t145;
t90 = t217 * t141 - t139;
t31 = qJD(4) * t90 + t217 * t127 + t214 * t226;
t317 = pkin(3) * t180;
t56 = pkin(4) * t120 - pkin(9) * t119 + t317;
t86 = pkin(9) * t278 + t91;
t272 = -pkin(2) * t212 - pkin(1);
t199 = qJD(1) * t272 + qJD(2);
t167 = -pkin(3) * t186 + t199;
t92 = -pkin(4) * t269 - pkin(9) * t238 + t167;
t7 = t213 * t56 + t216 * t31 + t92 * t280 - t281 * t86;
t38 = t213 * t92 + t216 * t86;
t8 = -qJD(5) * t38 - t213 * t31 + t216 * t56;
t352 = -t8 * t213 + t216 * t7;
t2 = qJ(6) * t120 + qJD(6) * t156 + t7;
t37 = -t213 * t86 + t216 * t92;
t343 = qJD(6) - t37;
t27 = -pkin(5) * t156 + t343;
t28 = qJ(6) * t156 + t38;
t4 = -pkin(5) * t120 - t8;
t351 = t2 * t216 + t213 * t4 + t27 * t280 - t28 * t281;
t255 = -Ifges(6,2) * t213 + t299;
t260 = mrSges(7,1) * t213 - mrSges(7,3) * t216;
t262 = mrSges(6,1) * t213 + mrSges(6,2) * t216;
t85 = -pkin(4) * t278 - t90;
t49 = t146 * pkin(5) - t147 * qJ(6) + t85;
t350 = t255 * t331 + t85 * t262 + t49 * t260 - (t213 * t38 + t216 * t37) * mrSges(6,3);
t330 = t146 / 0.2e1;
t329 = -t147 / 0.2e1;
t326 = -t156 / 0.2e1;
t295 = Ifges(5,2) * t269;
t349 = t295 / 0.2e1;
t348 = pkin(5) * t238;
t183 = pkin(5) * t281 - qJ(6) * t280 - qJD(6) * t213;
t246 = pkin(5) * t213 - qJ(6) * t216;
t344 = t246 * t269;
t346 = -t91 - t344 + t183;
t345 = qJ(6) * t238;
t168 = -t218 * t201 - t202 * t215;
t153 = -pkin(8) * t196 + t168;
t169 = -t215 * t201 + t218 * t202;
t154 = pkin(8) * t195 + t169;
t110 = t153 * t214 + t154 * t217;
t164 = t195 * t214 + t196 * t217;
t174 = -pkin(3) * t195 + t272;
t237 = t217 * t195 - t196 * t214;
t111 = -pkin(4) * t237 - pkin(9) * t164 + t174;
t342 = t216 * t110 + t213 * t111;
t341 = t217 * t153 - t154 * t214;
t155 = Ifges(5,4) * t269;
t251 = Ifges(6,5) * t216 - Ifges(6,6) * t213;
t229 = t156 * t251;
t253 = Ifges(7,4) * t216 + Ifges(7,6) * t213;
t230 = t156 * t253;
t297 = Ifges(7,5) * t213;
t257 = Ifges(7,1) * t216 + t297;
t231 = t147 * t257;
t300 = Ifges(6,4) * t213;
t259 = Ifges(6,1) * t216 - t300;
t232 = t147 * t259;
t249 = Ifges(7,3) * t213 + t296;
t233 = t146 * t249;
t319 = -t216 / 0.2e1;
t320 = t213 / 0.2e1;
t321 = -t213 / 0.2e1;
t303 = Ifges(5,1) * t238;
t339 = t155 / 0.2e1 + t303 / 0.2e1;
t142 = Ifges(7,5) * t147;
t71 = Ifges(7,6) * t156 + Ifges(7,3) * t146 + t142;
t301 = Ifges(6,4) * t147;
t74 = -Ifges(6,2) * t146 + Ifges(6,6) * t156 + t301;
t220 = (t213 * t28 - t216 * t27) * mrSges(7,2) + t90 * mrSges(5,3) - Ifges(5,5) * t278 - t233 / 0.2e1 - t232 / 0.2e1 - t231 / 0.2e1 - t230 / 0.2e1 - t229 / 0.2e1 - t167 * mrSges(5,2) + t71 * t321 + t74 * t320 + t347 * t319 - t339 - t350;
t340 = -t155 / 0.2e1 + t220;
t338 = (m(3) * qJ(2) + mrSges(3,3)) * (t211 ^ 2 + t212 ^ 2);
t150 = -t201 * t282 + t204 + (-qJD(2) * t211 - qJD(3) * t202) * t215;
t133 = -pkin(8) * t189 + t150;
t151 = -qJD(3) * t169 - t228;
t134 = -pkin(8) * t188 + t151;
t45 = qJD(4) * t341 + t133 * t217 + t134 * t214;
t125 = qJD(4) * t237 + t188 * t217 - t189 * t214;
t126 = qJD(4) * t164 + t188 * t214 + t217 * t189;
t316 = pkin(3) * t189;
t61 = pkin(4) * t126 - pkin(9) * t125 + t316;
t13 = -qJD(5) * t342 - t213 * t45 + t216 * t61;
t123 = pkin(4) * t238 - pkin(9) * t269;
t275 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t276 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t277 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t337 = t275 * t146 + t277 * t147 + t276 * t156 + t167 * mrSges(5,1) + t37 * mrSges(6,1) - t27 * mrSges(7,1) - t38 * mrSges(6,2) - t91 * mrSges(5,3) + t28 * mrSges(7,3) - Ifges(5,4) * t238 - Ifges(5,6) * t278 - Ifges(6,6) * t330 - Ifges(7,6) * t331 - t349 - t357 * t329 - (Ifges(7,2) + Ifges(6,3)) * t326;
t336 = t82 / 0.2e1;
t335 = -t83 / 0.2e1;
t334 = t83 / 0.2e1;
t332 = t120 / 0.2e1;
t328 = t147 / 0.2e1;
t324 = -t187 / 0.2e1;
t323 = t188 / 0.2e1;
t322 = -t189 / 0.2e1;
t318 = t216 / 0.2e1;
t315 = pkin(3) * t214;
t314 = pkin(3) * t217;
t40 = -mrSges(7,2) * t83 + mrSges(7,3) * t120;
t43 = -mrSges(6,2) * t120 - mrSges(6,3) * t83;
t307 = t40 + t43;
t41 = mrSges(6,1) * t120 - mrSges(6,3) * t82;
t42 = -t120 * mrSges(7,1) + t82 * mrSges(7,2);
t306 = -t41 + t42;
t94 = t144 * t217 - t139;
t99 = pkin(3) * t187 + t123;
t51 = t213 * t99 + t216 * t94;
t305 = mrSges(6,3) * t146;
t304 = mrSges(6,3) * t147;
t294 = t341 * t32;
t293 = t187 * Ifges(4,4);
t53 = t213 * t123 + t216 * t90;
t292 = -mrSges(5,1) * t278 + mrSges(6,1) * t146 + mrSges(6,2) * t147 + mrSges(5,3) * t238;
t288 = t213 * t217;
t287 = t216 * t217;
t104 = -mrSges(7,2) * t146 + mrSges(7,3) * t156;
t105 = -mrSges(6,2) * t156 - t305;
t285 = t104 + t105;
t106 = mrSges(6,1) * t156 - t304;
t107 = -mrSges(7,1) * t156 + mrSges(7,2) * t147;
t284 = t106 - t107;
t270 = t120 * mrSges(5,1) + t119 * mrSges(5,2);
t93 = t144 * t214 + t140;
t265 = -t2 * t213 + t216 * t4;
t264 = -t213 * t7 - t216 * t8;
t263 = mrSges(6,1) * t216 - mrSges(6,2) * t213;
t261 = mrSges(7,1) * t216 + mrSges(7,3) * t213;
t254 = Ifges(6,2) * t216 + t300;
t248 = -Ifges(7,3) * t216 + t297;
t247 = pkin(5) * t216 + qJ(6) * t213;
t244 = t213 * t27 + t216 * t28;
t241 = t213 * t37 - t216 * t38;
t50 = -t213 * t94 + t216 * t99;
t52 = t123 * t216 - t213 * t90;
t58 = -t110 * t213 + t111 * t216;
t200 = -pkin(4) - t247;
t12 = -t110 * t281 + t111 * t280 + t213 * t61 + t216 * t45;
t227 = t8 * mrSges(6,1) - t4 * mrSges(7,1) - t7 * mrSges(6,2) + t2 * mrSges(7,3);
t46 = qJD(4) * t110 + t133 * t214 - t217 * t134;
t223 = t307 * t216 + t306 * t213 + (-t213 * t285 - t216 * t284) * qJD(5) + m(6) * (-t280 * t37 - t281 * t38 + t352) + m(7) * t351;
t11 = pkin(5) * t83 - qJ(6) * t82 - qJD(6) * t147 + t32;
t19 = t82 * Ifges(7,5) + t120 * Ifges(7,6) + t83 * Ifges(7,3);
t20 = t82 * Ifges(6,4) - t83 * Ifges(6,2) + t120 * Ifges(6,6);
t222 = -t31 * mrSges(5,2) + Ifges(5,5) * t119 - Ifges(5,6) * t120 - t11 * t261 + t19 * t319 + t20 * t318 + t248 * t334 + t254 * t335 + t353 * t336 + t354 * t332 + t356 * t320 + (-t263 - mrSges(5,1)) * t32 + (-t74 / 0.2e1 + t71 / 0.2e1) * t281 + t347 * t280 / 0.2e1 + t352 * mrSges(6,3) + t350 * qJD(5) + t351 * mrSges(7,2) + (t231 + t230 + t229 + t233 + t232) * qJD(5) / 0.2e1;
t191 = t200 - t314;
t181 = Ifges(4,4) * t186;
t175 = t179 * mrSges(4,2);
t173 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t187;
t172 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t186;
t170 = qJD(4) * t315 + t183;
t158 = t187 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t181;
t157 = t186 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t293;
t148 = -mrSges(5,2) * t278 + mrSges(5,3) * t269;
t122 = -mrSges(5,1) * t269 + mrSges(5,2) * t238;
t117 = Ifges(7,2) * t120;
t115 = Ifges(6,3) * t120;
t96 = mrSges(7,1) * t146 - mrSges(7,3) * t147;
t95 = pkin(5) * t147 + qJ(6) * t146;
t81 = Ifges(7,4) * t82;
t80 = Ifges(6,5) * t82;
t79 = Ifges(6,6) * t83;
t78 = Ifges(7,6) * t83;
t62 = t164 * t246 - t341;
t57 = t93 + t344;
t48 = pkin(5) * t237 - t58;
t47 = -qJ(6) * t237 + t342;
t36 = -t52 - t348;
t35 = t53 + t345;
t34 = -t50 - t348;
t33 = t51 + t345;
t23 = mrSges(7,1) * t83 - mrSges(7,3) * t82;
t14 = t246 * t125 + (qJD(5) * t247 - qJD(6) * t216) * t164 + t46;
t9 = -pkin(5) * t126 - t13;
t6 = qJ(6) * t126 - qJD(6) * t237 + t12;
t1 = [t342 * t43 + (t137 * t195 - t138 * t196 - t165 * t188 - t166 * t189 - t168 * t179 - t169 * t180) * mrSges(4,3) + t272 * (t180 * mrSges(4,1) + t175) + (-t195 * t180 + t186 * t322) * Ifges(4,2) + (t195 * t179 - t196 * t180 + t186 * t323 + t187 * t322) * Ifges(4,4) + m(6) * (t12 * t38 + t13 * t37 + t342 * t7 + t46 * t85 + t58 * t8 - t294) + m(5) * (-t294 + t110 * t31 + t45 * t91 - t46 * t90 + (t167 * t189 + t174 * t180) * pkin(3)) - (mrSges(5,1) * t317 + t80 / 0.2e1 - t79 / 0.2e1 + t115 / 0.2e1 + t81 / 0.2e1 + t117 / 0.2e1 + t78 / 0.2e1 - Ifges(5,4) * t119 - t31 * mrSges(5,3) + t275 * t83 + t277 * t82 + (Ifges(5,2) + t276) * t120 + t227) * t237 + (-t110 * t120 - t119 * t341) * mrSges(5,3) - t341 * t24 + t62 * t23 + t58 * t41 + 0.2e1 * t338 * t279 + t47 * t40 + t48 * t42 + m(4) * (t137 * t169 + t138 * t168 + t150 * t166 + t151 * t165) + m(7) * (t11 * t62 + t14 * t49 + t2 * t47 + t27 * t9 + t28 * t6 + t4 * t48) + qJD(3) * (Ifges(4,5) * t188 - Ifges(4,6) * t189) / 0.2e1 + t199 * (mrSges(4,1) * t189 + mrSges(4,2) * t188) + (-t295 / 0.2e1 + t337) * t126 + (t255 * t335 + t19 * t320 + t11 * t260 + t249 * t334 + mrSges(5,2) * t317 + Ifges(5,1) * t119 - Ifges(5,4) * t120 + (mrSges(5,3) + t262) * t32 + t264 * mrSges(6,3) + t265 * mrSges(7,2) + (-mrSges(7,2) * t244 + mrSges(6,3) * t241 + t248 * t331 + t254 * t330 + t261 * t49 + t263 * t85 + t319 * t74 + t326 * t354 + t329 * t353) * qJD(5) + (t257 + t259) * t336 + (t253 + t251) * t332 + (qJD(5) * t347 + t20) * t321 + (qJD(5) * t71 + t356) * t318) * t164 + t157 * t322 + t158 * t323 + t14 * t96 + t174 * t270 + t6 * t104 + t12 * t105 + (-t220 + t339) * t125 + t13 * t106 + t9 * t107 + t45 * t148 + t292 * t46 + t150 * t172 + t151 * t173 + t122 * t316 + (t196 * t179 + t187 * t323) * Ifges(4,1); -t269 * t148 - t186 * t172 + t187 * t173 + t175 - (-m(5) * pkin(3) - mrSges(4,1)) * t180 + (-t96 - t292) * t238 + (t156 * t285 - t306) * t216 + (-t156 * t284 + t307) * t213 - m(4) * (-t165 * t187 + t166 * t186) - m(5) * (-t238 * t90 + t269 * t91) + t270 - t338 * qJD(1) ^ 2 + (t156 * t244 - t238 * t49 - t265) * m(7) + (-t156 * t241 - t238 * t85 - t264) * m(6); -(-Ifges(4,2) * t187 + t158 + t181) * t186 / 0.2e1 + (-t337 + t349) * t238 + t222 - m(5) * (-t90 * t93 + t91 * t94) + (t165 * t186 + t166 * t187) * mrSges(4,3) - m(7) * (t27 * t34 + t28 * t33 + t49 * t57) - m(6) * (t37 * t50 + t38 * t51 + t85 * t93) + t355 * (-pkin(4) - t314) + (-t57 + t170) * t96 + (-t187 * t122 + (-t119 * t217 - t120 * t214) * mrSges(5,3) + (t292 * t214 + (-t213 * t284 + t216 * t285 + t148) * t217 + m(6) * (t214 * t85 + t287 * t38 - t288 * t37) + m(7) * (t27 * t288 + t28 * t287)) * qJD(4) + (t214 * t31 - t217 * t32 + 0.2e1 * t167 * t324 + (-t214 * t90 + t217 * t91) * qJD(4)) * m(5)) * pkin(3) + (Ifges(4,1) * t186 - t293) * t324 + m(7) * (t11 * t191 + t170 * t49) - t33 * t104 - t51 * t105 - t50 * t106 - t34 * t107 + (-t303 / 0.2e1 + t340) * t269 - t137 * mrSges(4,2) + t138 * mrSges(4,1) - t94 * t148 - t292 * t93 - t165 * t172 + t166 * t173 + Ifges(4,5) * t179 - Ifges(4,6) * t180 + t187 * t157 / 0.2e1 - qJD(3) * (Ifges(4,5) * t186 - Ifges(4,6) * t187) / 0.2e1 + t223 * (pkin(9) + t315) + t191 * t23 - t199 * (mrSges(4,1) * t187 + mrSges(4,2) * t186); t223 * pkin(9) - t292 * t91 - m(6) * (t37 * t52 + t38 * t53 + t85 * t91) + t346 * t96 + t222 - t337 * t238 + ((Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t238 + t340) * t269 - t35 * t104 - t53 * t105 - t52 * t106 - t36 * t107 - t90 * t148 + t200 * t23 - t355 * pkin(4) + (t11 * t200 - t27 * t36 - t28 * t35 + t346 * t49) * m(7); (t146 * t27 + t147 * t28) * mrSges(7,2) + t117 + t115 + qJ(6) * t40 - pkin(5) * t42 + t227 + t81 + t80 - t79 + t78 + (Ifges(7,3) * t147 - t298) * t331 + t74 * t328 - t95 * t96 + qJD(6) * t104 - t49 * (mrSges(7,1) * t147 + mrSges(7,3) * t146) - t85 * (mrSges(6,1) * t147 - mrSges(6,2) * t146) + (t284 + t304) * t38 + (-t285 - t305) * t37 + (-t357 * t146 - t147 * t358) * t326 + (-pkin(5) * t4 + qJ(6) * t2 - t27 * t38 + t28 * t343 - t49 * t95) * m(7) + (-Ifges(6,2) * t147 - t143 + t347) * t330 + (-t146 * t359 + t142 - t301 + t71) * t329; -t156 * t104 + t147 * t96 + 0.2e1 * (t4 / 0.2e1 + t49 * t328 + t28 * t326) * m(7) + t42;];
tauc  = t1(:);
