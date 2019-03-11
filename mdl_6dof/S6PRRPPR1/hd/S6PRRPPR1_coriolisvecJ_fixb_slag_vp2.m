% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:58:33
% EndTime: 2019-03-08 20:58:56
% DurationCPUTime: 11.84s
% Computational Cost: add. (6660->547), mult. (17720->778), div. (0->0), fcn. (13455->12), ass. (0->263)
t216 = sin(qJ(2));
t211 = sin(pkin(6));
t262 = qJD(1) * t211;
t245 = t216 * t262;
t215 = sin(qJ(3));
t257 = qJD(3) * t215;
t349 = pkin(3) * t257 - t245;
t210 = sin(pkin(11));
t218 = cos(qJ(3));
t276 = cos(pkin(11));
t186 = t210 * t218 + t215 * t276;
t172 = t186 * qJD(3);
t236 = t276 * t218;
t222 = -t210 * t215 + t236;
t174 = t222 * qJD(3);
t348 = pkin(4) * t172 - qJ(5) * t174 - qJD(5) * t186 + t349;
t296 = -qJ(4) - pkin(8);
t238 = qJD(3) * t296;
t255 = qJD(4) * t218;
t168 = t215 * t238 + t255;
t169 = -qJD(4) * t215 + t218 * t238;
t219 = cos(qJ(2));
t244 = t219 * t262;
t339 = t168 * t276 + t210 * t169 - t222 * t244;
t209 = sin(pkin(12));
t212 = cos(pkin(12));
t342 = -t339 * t209 + t348 * t212;
t341 = t348 * t209 + t339 * t212;
t300 = pkin(9) * t212;
t347 = pkin(5) * t172 - t174 * t300 + t342;
t267 = t174 * t209;
t346 = pkin(9) * t267 - t341;
t260 = qJD(2) * t215;
t171 = -qJD(2) * t236 + t210 * t260;
t307 = -t171 / 0.2e1;
t345 = Ifges(5,6) * qJD(3) / 0.2e1;
t214 = sin(qJ(6));
t217 = cos(qJ(6));
t206 = -pkin(3) * t218 - pkin(2);
t133 = -pkin(4) * t222 - qJ(5) * t186 + t206;
t196 = t296 * t215;
t197 = t296 * t218;
t148 = t210 * t196 - t197 * t276;
t71 = t212 * t133 - t148 * t209;
t50 = -pkin(5) * t222 - t186 * t300 + t71;
t266 = t186 * t209;
t72 = t209 * t133 + t212 * t148;
t54 = -pkin(9) * t266 + t72;
t17 = t214 * t50 + t217 * t54;
t344 = -qJD(6) * t17 + t346 * t214 + t217 * t347;
t16 = -t214 * t54 + t217 * t50;
t343 = qJD(6) * t16 + t214 * t347 - t346 * t217;
t117 = t168 * t210 - t276 * t169;
t145 = t186 * t244;
t340 = t117 - t145;
t173 = t186 * qJD(2);
t151 = qJD(3) * t209 + t173 * t212;
t235 = t212 * qJD(3) - t173 * t209;
t338 = -t151 * t214 + t217 * t235;
t88 = t151 * t217 + t214 * t235;
t165 = qJD(2) * t206 + qJD(4) - t244;
t282 = t151 * Ifges(6,5);
t283 = t235 * Ifges(6,6);
t290 = Ifges(5,4) * t173;
t213 = cos(pkin(6));
t263 = t213 * t218;
t201 = qJD(1) * t263;
t191 = qJD(2) * pkin(8) + t245;
t234 = qJ(4) * qJD(2) + t191;
t143 = -t215 * t234 + t201;
t139 = qJD(3) * pkin(3) + t143;
t261 = qJD(1) * t215;
t243 = t213 * t261;
t144 = t218 * t234 + t243;
t237 = t276 * t144;
t75 = t210 * t139 + t237;
t70 = qJD(3) * qJ(5) + t75;
t99 = pkin(4) * t171 - qJ(5) * t173 + t165;
t30 = -t209 * t70 + t212 * t99;
t31 = t209 * t99 + t212 * t70;
t20 = pkin(5) * t171 - pkin(9) * t151 + t30;
t21 = pkin(9) * t235 + t31;
t5 = t20 * t217 - t21 * t214;
t6 = t20 * t214 + t21 * t217;
t337 = t31 * mrSges(6,2) + t6 * mrSges(7,2) + Ifges(5,2) * t307 + t290 / 0.2e1 + t345 - t165 * mrSges(5,1) - t30 * mrSges(6,1) - t5 * mrSges(7,1) - t282 / 0.2e1 - t283 / 0.2e1;
t164 = qJD(2) * t174;
t224 = t209 * t214 - t212 * t217;
t44 = qJD(6) * t338 - t164 * t224;
t316 = t44 / 0.2e1;
t187 = t209 * t217 + t212 * t214;
t45 = -qJD(6) * t88 - t164 * t187;
t315 = t45 / 0.2e1;
t163 = qJD(2) * t172;
t310 = t163 / 0.2e1;
t301 = pkin(3) * t210;
t203 = qJ(5) + t301;
t297 = pkin(9) + t203;
t181 = t297 * t209;
t182 = t297 * t212;
t128 = -t181 * t217 - t182 * t214;
t268 = t171 * t212;
t252 = pkin(3) * t260;
t116 = pkin(4) * t173 + qJ(5) * t171 + t252;
t135 = t210 * t144;
t81 = t143 * t276 - t135;
t36 = t212 * t116 - t209 * t81;
t22 = pkin(5) * t173 + pkin(9) * t268 + t36;
t269 = t171 * t209;
t37 = t209 * t116 + t212 * t81;
t27 = pkin(9) * t269 + t37;
t336 = -qJD(5) * t224 + qJD(6) * t128 - t214 * t22 - t217 * t27;
t129 = -t181 * t214 + t182 * t217;
t335 = -qJD(5) * t187 - qJD(6) * t129 + t214 * t27 - t217 * t22;
t270 = t164 * t212;
t271 = t164 * t209;
t103 = mrSges(6,1) * t271 + mrSges(6,2) * t270;
t14 = -t45 * mrSges(7,1) + t44 * mrSges(7,2);
t334 = t103 + t14;
t292 = mrSges(5,3) * t173;
t333 = qJD(3) * mrSges(5,1) + mrSges(6,1) * t235 - mrSges(6,2) * t151 - t292;
t332 = Ifges(5,5) * qJD(3);
t254 = qJD(2) * qJD(3);
t240 = t215 * t254;
t110 = t224 * t171;
t175 = t224 * qJD(6);
t328 = -t175 - t110;
t109 = t187 * t171;
t176 = t187 * qJD(6);
t327 = -t176 - t109;
t104 = -mrSges(6,2) * t171 + mrSges(6,3) * t235;
t105 = mrSges(6,1) * t171 - mrSges(6,3) * t151;
t326 = t104 * t212 - t105 * t209;
t264 = t211 * t219;
t241 = qJD(2) * t264;
t233 = t218 * t241;
t130 = qJD(1) * t233 + qJD(3) * t201 - t191 * t257;
t102 = (-qJ(4) * t257 + t255) * qJD(2) + t130;
t156 = t191 * t218 + t243;
t256 = qJD(3) * t218;
t221 = -t156 * qJD(3) + (-qJ(4) * t256 + (-qJD(4) - t244) * t215) * qJD(2);
t47 = t276 * t102 + t210 * t221;
t40 = qJD(3) * qJD(5) + t47;
t259 = qJD(2) * t216;
t242 = t211 * t259;
t170 = pkin(3) * t240 + qJD(1) * t242;
t64 = pkin(4) * t163 - qJ(5) * t164 - qJD(5) * t173 + t170;
t18 = -t209 * t40 + t212 * t64;
t19 = t209 * t64 + t212 * t40;
t325 = -t18 * t209 + t19 * t212;
t13 = pkin(5) * t163 - pkin(9) * t270 + t18;
t15 = -pkin(9) * t271 + t19;
t1 = qJD(6) * t5 + t13 * t214 + t15 * t217;
t2 = -qJD(6) * t6 + t13 * t217 - t15 * t214;
t324 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t44 + Ifges(7,6) * t45;
t286 = Ifges(6,2) * t209;
t288 = Ifges(6,4) * t212;
t229 = -t286 + t288;
t289 = Ifges(6,4) * t209;
t230 = Ifges(6,1) * t212 - t289;
t231 = mrSges(6,1) * t209 + mrSges(6,2) * t212;
t303 = t212 / 0.2e1;
t304 = -t209 / 0.2e1;
t74 = t139 * t276 - t135;
t66 = -qJD(3) * pkin(4) + qJD(5) - t74;
t322 = (t151 * Ifges(6,1) + Ifges(6,4) * t235 + t171 * Ifges(6,5)) * t303 + (t151 * Ifges(6,4) + Ifges(6,2) * t235 + Ifges(6,6) * t171) * t304 + t165 * mrSges(5,2) + t66 * t231 + t151 * t230 / 0.2e1 + t235 * t229 / 0.2e1;
t220 = qJD(2) ^ 2;
t320 = Ifges(7,4) * t316 + Ifges(7,2) * t315 + Ifges(7,6) * t310;
t319 = Ifges(7,1) * t316 + Ifges(7,4) * t315 + Ifges(7,5) * t310;
t167 = qJD(6) + t171;
t302 = Ifges(7,4) * t88;
t25 = Ifges(7,2) * t338 + t167 * Ifges(7,6) + t302;
t318 = t25 / 0.2e1;
t85 = Ifges(7,4) * t338;
t26 = t88 * Ifges(7,1) + t167 * Ifges(7,5) + t85;
t317 = t26 / 0.2e1;
t314 = -t338 / 0.2e1;
t313 = t338 / 0.2e1;
t312 = -t88 / 0.2e1;
t311 = t88 / 0.2e1;
t309 = -t167 / 0.2e1;
t308 = t167 / 0.2e1;
t306 = t171 / 0.2e1;
t305 = -t173 / 0.2e1;
t299 = t338 * Ifges(7,6);
t298 = t88 * Ifges(7,5);
t295 = mrSges(5,3) * t163;
t294 = mrSges(5,3) * t164;
t293 = mrSges(5,3) * t171;
t291 = Ifges(4,4) * t215;
t265 = t211 * t216;
t177 = -t215 * t265 + t263;
t178 = t213 * t215 + t218 * t265;
t122 = -t177 * t276 + t178 * t210;
t46 = t102 * t210 - t276 * t221;
t285 = t122 * t46;
t147 = -t276 * t196 - t197 * t210;
t284 = t147 * t46;
t281 = t167 * Ifges(7,3);
t280 = t173 * Ifges(5,1);
t275 = Ifges(4,5) * qJD(3);
t274 = Ifges(4,6) * qJD(3);
t258 = qJD(2) * t218;
t34 = -mrSges(7,1) * t338 + mrSges(7,2) * t88;
t253 = -t34 + t333;
t250 = mrSges(4,3) * t260;
t249 = mrSges(4,3) * t258;
t246 = t276 * pkin(3);
t119 = t163 * mrSges(5,1) + t164 * mrSges(5,2);
t79 = t143 * t210 + t237;
t205 = -t246 - pkin(4);
t228 = Ifges(6,5) * t212 - Ifges(6,6) * t209;
t227 = t18 * t212 + t19 * t209;
t226 = t209 * t30 - t212 * t31;
t123 = t210 * t177 + t178 * t276;
t100 = -t123 * t209 - t212 * t264;
t101 = t123 * t212 - t209 * t264;
t38 = t100 * t217 - t101 * t214;
t39 = t100 * t214 + t101 * t217;
t131 = -t191 * t256 + (-qJD(3) * t213 - t241) * t261;
t225 = t130 * t218 - t131 * t215;
t223 = (mrSges(4,1) * t215 + mrSges(4,2) * t218) * qJD(2);
t207 = Ifges(4,4) * t258;
t195 = -qJD(3) * mrSges(4,2) + t249;
t194 = qJD(3) * mrSges(4,1) - t250;
t193 = -t212 * pkin(5) + t205;
t192 = -qJD(2) * pkin(2) - t244;
t183 = qJD(3) * t223;
t180 = Ifges(4,1) * t260 + t207 + t275;
t179 = t274 + (t218 * Ifges(4,2) + t291) * qJD(2);
t166 = Ifges(5,4) * t171;
t161 = Ifges(7,3) * t163;
t157 = -qJD(3) * mrSges(5,2) - t293;
t155 = -t191 * t215 + t201;
t142 = qJD(3) * t177 + t233;
t141 = -qJD(3) * t178 - t215 * t241;
t132 = mrSges(5,1) * t171 + mrSges(5,2) * t173;
t125 = -t166 + t280 + t332;
t121 = t224 * t186;
t120 = t187 * t186;
t113 = pkin(5) * t266 + t147;
t112 = mrSges(6,1) * t163 - mrSges(6,3) * t270;
t111 = -mrSges(6,2) * t163 - mrSges(6,3) * t271;
t82 = pkin(5) * t267 + t117;
t80 = t210 * t141 + t142 * t276;
t78 = -t141 * t276 + t142 * t210;
t77 = t163 * Ifges(6,5) + t164 * t230;
t76 = t163 * Ifges(6,6) + t164 * t229;
t67 = t171 * Ifges(6,3) + t282 + t283;
t61 = mrSges(7,1) * t167 - mrSges(7,3) * t88;
t60 = -mrSges(7,2) * t167 + mrSges(7,3) * t338;
t59 = t209 * t242 + t212 * t80;
t58 = -t209 * t80 + t212 * t242;
t57 = -pkin(5) * t269 + t79;
t56 = -t174 * t187 + t175 * t186;
t55 = -t174 * t224 - t176 * t186;
t51 = -pkin(5) * t235 + t66;
t32 = pkin(5) * t271 + t46;
t29 = -mrSges(7,2) * t163 + mrSges(7,3) * t45;
t28 = mrSges(7,1) * t163 - mrSges(7,3) * t44;
t24 = t281 + t298 + t299;
t10 = -qJD(6) * t39 - t214 * t59 + t217 * t58;
t9 = qJD(6) * t38 + t214 * t58 + t217 * t59;
t3 = [-t123 * t295 + t10 * t61 + t100 * t112 + t101 * t111 + t59 * t104 + t58 * t105 + t141 * t194 + t142 * t195 + t80 * t157 + t38 * t28 + t39 * t29 + t9 * t60 + (-t177 * t218 - t178 * t215) * mrSges(4,3) * t254 - t253 * t78 + (t294 + t334) * t122 + ((-mrSges(3,2) * t220 - t119 - t183) * t219 + (-mrSges(3,1) * t220 + (t132 + qJD(2) * (-mrSges(4,1) * t218 + mrSges(4,2) * t215)) * qJD(2)) * t216) * t211 + m(5) * (t285 + t123 * t47 - t74 * t78 + t75 * t80 + (t165 * t259 - t170 * t219) * t211) + m(4) * (t130 * t178 + t131 * t177 + t141 * t155 + t142 * t156 + (t192 - t244) * t242) + m(7) * (t1 * t39 + t10 * t5 + t122 * t32 + t2 * t38 + t51 * t78 + t6 * t9) + m(6) * (t100 * t18 + t101 * t19 + t30 * t58 + t31 * t59 + t66 * t78 + t285); (t280 / 0.2e1 + t228 * t306 + t125 / 0.2e1 - t74 * mrSges(5,3) + (-t209 * t31 - t212 * t30) * mrSges(6,3) + t322) * t174 + (t24 / 0.2e1 + t67 / 0.2e1 + t281 / 0.2e1 + t299 / 0.2e1 + t298 / 0.2e1 - t75 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t171 - t337) * t172 - (t18 * mrSges(6,1) - t19 * mrSges(6,2) + t170 * mrSges(5,1) + t161 / 0.2e1 - t47 * mrSges(5,3) + t228 * t164 + (Ifges(5,2) + Ifges(6,3) + Ifges(7,3) / 0.2e1) * t163 + t324) * t222 + (-t163 * t186 + t164 * t222 + t172 * t305 + t174 * t307) * Ifges(5,4) + (((-t155 * t218 - t156 * t215) * qJD(3) + t225) * pkin(8) - (pkin(2) * t259 + t192 * t216 + (-t155 * t215 + t156 * t218) * t219) * t262) * m(4) + (t148 * t47 + t165 * t349 + t170 * t206 + t339 * t75 - t340 * t74 + t284) * m(5) + t339 * t157 + t225 * mrSges(4,3) + (t147 * t164 - t148 * t163) * mrSges(5,3) + (-0.3e1 / 0.2e1 * t215 ^ 2 + 0.3e1 / 0.2e1 * t218 ^ 2) * Ifges(4,4) * t254 + (-Ifges(7,5) * t121 - Ifges(7,6) * t120) * t310 + (-Ifges(7,4) * t121 - Ifges(7,2) * t120) * t315 + (-Ifges(7,1) * t121 - Ifges(7,4) * t120) * t316 + (-t1 * t120 + t121 * t2 - t5 * t55 + t56 * t6) * mrSges(7,3) + t32 * (mrSges(7,1) * t120 - mrSges(7,2) * t121) + t206 * t119 - pkin(2) * t183 + t147 * t103 + t72 * t111 + t71 * t112 + t113 * t14 + t82 * t34 + t51 * (-mrSges(7,1) * t56 + mrSges(7,2) * t55) + t17 * t29 + t16 * t28 - t333 * t117 + ((t194 * t215 - t195 * t218) * t219 - t216 * t132) * t262 + (Ifges(7,5) * t55 + Ifges(7,6) * t56) * t308 + (Ifges(7,1) * t55 + Ifges(7,4) * t56) * t311 + (Ifges(7,4) * t55 + Ifges(7,2) * t56) * t313 + t341 * t104 + t342 * t105 + (t18 * t71 + t19 * t72 + t30 * t342 + t31 * t341 + t340 * t66 + t284) * m(6) + t343 * t60 + t344 * t61 + (t1 * t17 + t113 * t32 + t16 * t2 + t343 * t6 + (-t145 + t82) * t51 + t344 * t5) * m(7) + t55 * t317 + t56 * t318 - t121 * t319 - t120 * t320 + t253 * t145 + (Ifges(5,5) * t174 / 0.2e1 - Ifges(5,6) * t172 / 0.2e1 + (t180 / 0.2e1 + t192 * mrSges(4,2) - t155 * mrSges(4,3) - pkin(8) * t194 + t275 / 0.2e1) * t218 + (-t179 / 0.2e1 + t192 * mrSges(4,1) + pkin(3) * t132 - t156 * mrSges(4,3) - pkin(8) * t195 - t274 / 0.2e1 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t258) * t215) * qJD(3) + (t77 * t303 + t76 * t304 + t228 * t310 + t170 * mrSges(5,2) + (mrSges(5,3) + t231) * t46 - t227 * mrSges(6,3) + (Ifges(6,1) * t212 ^ 2 / 0.2e1 + Ifges(5,1) + (-t288 + t286 / 0.2e1) * t209) * t164) * t186; (t249 - t195) * t155 + ((t210 * t47 - t276 * t46) * pkin(3) - t165 * t252 + t74 * t79 - t75 * t81) * m(5) + t254 * Ifges(4,5) * t218 / 0.2e1 + (-mrSges(6,1) * t212 + mrSges(6,2) * t209 - mrSges(5,1)) * t46 + (t125 - t166) * t306 + (t212 * t111 - t209 * t112) * t203 + (-t1 * t224 - t187 * t2 + t327 * t6 - t328 * t5) * mrSges(7,3) + (Ifges(6,5) * t209 + Ifges(7,5) * t187 + Ifges(6,6) * t212 - Ifges(7,6) * t224) * t310 + t32 * (mrSges(7,1) * t224 + mrSges(7,2) * t187) + (Ifges(7,4) * t187 - Ifges(7,2) * t224) * t315 + (Ifges(7,1) * t187 - Ifges(7,4) * t224) * t316 - t224 * t320 + t75 * t292 + t76 * t303 + (Ifges(7,5) * t110 + Ifges(7,6) * t109) * t309 + (t332 / 0.2e1 - t228 * t307 + t322) * t171 + t333 * t79 - t192 * t223 + t326 * qJD(5) + (-mrSges(7,1) * t327 + mrSges(7,2) * t328) * t51 + (t250 + t194) * t156 + (Ifges(7,4) * t110 + Ifges(7,2) * t109) * t314 + (Ifges(7,1) * t110 + Ifges(7,4) * t109) * t312 + (-Ifges(5,1) * t171 + t24 - t290 + t67) * t305 - Ifges(4,6) * t240 / 0.2e1 + (-t268 * t30 - t269 * t31 + t325) * mrSges(6,3) + (-t226 * qJD(5) + t203 * t325 + t205 * t46 - t30 * t36 - t31 * t37 - t66 * t79) * m(6) + (-Ifges(7,5) * t175 - Ifges(7,6) * t176) * t308 + (-Ifges(7,1) * t175 - Ifges(7,4) * t176) * t311 + (-Ifges(7,4) * t175 - Ifges(7,2) * t176) * t313 - (-Ifges(4,2) * t260 + t180 + t207) * t258 / 0.2e1 + t209 * t77 / 0.2e1 + t205 * t103 + t193 * t14 + Ifges(5,5) * t164 - Ifges(5,6) * t163 - t81 * t157 + t128 * t28 + t129 * t29 - t130 * mrSges(4,2) + t131 * mrSges(4,1) - t109 * t25 / 0.2e1 - t110 * t26 / 0.2e1 - t37 * t104 - t36 * t105 - t57 * t34 - t47 * mrSges(5,2) - t215 * t220 * (Ifges(4,1) * t218 - t291) / 0.2e1 + t335 * t61 + t336 * t60 + (t1 * t129 + t128 * t2 + t193 * t32 + t335 * t5 + t336 * t6 - t51 * t57) * m(7) + (Ifges(7,5) * t312 - Ifges(5,2) * t306 + Ifges(7,6) * t314 + Ifges(6,3) * t307 + Ifges(7,3) * t309 + t337 + t345) * t173 - t175 * t317 - t176 * t318 + t187 * t319 - t132 * t252 + t179 * t260 / 0.2e1 + (Ifges(6,1) * t209 + t288) * t270 / 0.2e1 - (Ifges(6,2) * t212 + t289) * t271 / 0.2e1 - t74 * t293 - t246 * t294 - t295 * t301; t209 * t111 + t212 * t112 - t224 * t28 + t187 * t29 + t327 * t61 + t328 * t60 + t253 * t173 + (t157 + t326) * t171 + t119 + (t1 * t187 - t173 * t51 - t2 * t224 + t327 * t5 + t328 * t6) * m(7) + (-t171 * t226 - t173 * t66 + t227) * m(6) + (t171 * t75 + t173 * t74 + t170) * m(5); -t235 * t104 + t151 * t105 - t338 * t60 + t88 * t61 + (-t338 * t6 + t5 * t88 + t32) * m(7) + (t151 * t30 - t235 * t31 + t46) * m(6) + t334; t161 - t51 * (mrSges(7,1) * t88 + mrSges(7,2) * t338) + (Ifges(7,1) * t338 - t302) * t312 + t25 * t311 + (Ifges(7,5) * t338 - Ifges(7,6) * t88) * t309 - t5 * t60 + t6 * t61 + (t338 * t5 + t6 * t88) * mrSges(7,3) + (-Ifges(7,2) * t88 + t26 + t85) * t314 + t324;];
tauc  = t3(:);
