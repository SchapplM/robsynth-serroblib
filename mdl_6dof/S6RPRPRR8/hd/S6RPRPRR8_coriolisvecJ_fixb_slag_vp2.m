% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPRR8
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:57:56
% EndTime: 2019-03-09 03:58:14
% DurationCPUTime: 10.44s
% Computational Cost: add. (8888->562), mult. (19849->777), div. (0->0), fcn. (13677->8), ass. (0->269)
t198 = sin(qJ(5));
t196 = sin(pkin(10));
t294 = pkin(3) * t196;
t187 = pkin(8) + t294;
t292 = pkin(9) + t187;
t238 = qJD(5) * t292;
t199 = sin(qJ(3));
t271 = cos(pkin(10));
t236 = t271 * t199;
t296 = cos(qJ(3));
t245 = qJD(1) * t296;
t155 = -qJD(1) * t236 - t196 * t245;
t269 = t155 * t198;
t202 = -pkin(1) - pkin(7);
t179 = qJD(1) * t202 + qJD(2);
t262 = qJD(1) * t199;
t148 = -qJ(4) * t262 + t179 * t199;
t136 = t196 * t148;
t170 = t296 * t179;
t149 = -qJ(4) * t245 + t170;
t105 = t149 * t271 - t136;
t228 = t271 * t296;
t217 = qJD(1) * t228;
t154 = t196 * t262 - t217;
t232 = pkin(3) * t245;
t107 = -t154 * pkin(4) - t155 * pkin(8) + t232;
t201 = cos(qJ(5));
t52 = t201 * t105 + t198 * t107;
t361 = -pkin(9) * t269 + t198 * t238 + t52;
t268 = t155 * t201;
t51 = -t105 * t198 + t201 * t107;
t360 = pkin(5) * t154 + pkin(9) * t268 - t201 * t238 - t51;
t260 = qJD(5) * t198;
t359 = t260 - t269;
t200 = cos(qJ(6));
t197 = sin(qJ(6));
t130 = qJD(3) * t201 + t154 * t198;
t139 = qJD(3) * pkin(3) + t149;
t237 = t271 * t148;
t92 = t196 * t139 + t237;
t84 = qJD(3) * pkin(8) + t92;
t172 = pkin(3) * t262 + qJD(1) * qJ(2) + qJD(4);
t95 = -pkin(4) * t155 + pkin(8) * t154 + t172;
t50 = t198 * t95 + t201 * t84;
t37 = pkin(9) * t130 + t50;
t275 = t197 * t37;
t151 = qJD(5) - t155;
t131 = qJD(3) * t198 - t154 * t201;
t49 = -t198 * t84 + t201 * t95;
t36 = -pkin(9) * t131 + t49;
t31 = pkin(5) * t151 + t36;
t10 = t200 * t31 - t275;
t273 = t200 * t37;
t11 = t197 * t31 + t273;
t234 = t200 * t130 - t131 * t197;
t257 = qJD(1) * qJD(3);
t242 = t199 * t257;
t144 = -qJD(3) * t217 + t196 * t242;
t207 = -t196 * t296 - t236;
t145 = t207 * t257;
t85 = qJD(5) * t130 + t145 * t201;
t86 = -qJD(5) * t131 - t145 * t198;
t29 = qJD(6) * t234 + t197 * t86 + t200 * t85;
t75 = t130 * t197 + t131 * t200;
t30 = -qJD(6) * t75 - t197 * t85 + t200 * t86;
t255 = Ifges(7,5) * t29 + Ifges(7,6) * t30 - Ifges(7,3) * t144;
t295 = Ifges(7,4) * t75;
t143 = qJD(6) + t151;
t306 = -t143 / 0.2e1;
t313 = -t75 / 0.2e1;
t259 = qJD(5) * t201;
t244 = qJD(3) * t296;
t231 = t179 * t244;
t258 = t199 * qJD(4);
t128 = t231 + (-qJ(4) * t244 - t258) * qJD(1);
t243 = t296 * qJD(4);
t261 = qJD(3) * t199;
t247 = t179 * t261;
t204 = -t247 + (qJ(4) * t261 - t243) * qJD(1);
t70 = t128 * t271 + t196 * t204;
t194 = qJD(1) * qJD(2);
t230 = qJD(1) * t244;
t171 = pkin(3) * t230 + t194;
t79 = -pkin(4) * t144 - pkin(8) * t145 + t171;
t18 = t198 * t79 + t201 * t70 + t95 * t259 - t260 * t84;
t12 = pkin(9) * t86 + t18;
t19 = -qJD(5) * t50 - t198 * t70 + t201 * t79;
t9 = -pkin(5) * t144 - pkin(9) * t85 + t19;
t2 = qJD(6) * t10 + t12 * t200 + t197 * t9;
t3 = -qJD(6) * t11 - t12 * t197 + t200 * t9;
t328 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t91 = t139 * t271 - t136;
t83 = -qJD(3) * pkin(4) - t91;
t59 = -t130 * pkin(5) + t83;
t358 = t255 + t328 + (Ifges(7,5) * t234 - Ifges(7,6) * t75) * t306 + (t10 * t234 + t11 * t75) * mrSges(7,3) - t59 * (mrSges(7,1) * t75 + mrSges(7,2) * t234) + (Ifges(7,1) * t234 - t295) * t313;
t71 = Ifges(7,4) * t234;
t357 = -Ifges(7,2) * t75 + t71;
t356 = qJ(2) * (m(3) + m(4)) + mrSges(3,3);
t321 = t29 / 0.2e1;
t320 = t30 / 0.2e1;
t304 = -t144 / 0.2e1;
t156 = -qJD(3) * t228 + t196 * t261;
t353 = t156 / 0.2e1;
t163 = t292 * t198;
t164 = t292 * t201;
t119 = -t163 * t197 + t164 * t200;
t352 = -qJD(6) * t119 + t361 * t197 + t360 * t200;
t118 = -t163 * t200 - t164 * t197;
t351 = qJD(6) * t118 + t360 * t197 - t361 * t200;
t344 = t151 * Ifges(6,3);
t345 = t130 * Ifges(6,6);
t350 = t131 * Ifges(6,5) + t75 * Ifges(7,5) + Ifges(7,6) * t234 + t143 * Ifges(7,3) + t344 + t345;
t104 = t149 * t196 + t237;
t343 = t359 * pkin(5) - t104;
t291 = mrSges(5,3) * t145;
t48 = -mrSges(6,1) * t86 + mrSges(6,2) * t85;
t342 = t291 + t48;
t168 = t197 * t201 + t198 * t200;
t331 = qJD(5) + qJD(6);
t124 = t331 * t168;
t219 = t197 * t198 - t200 * t201;
t341 = -t168 * qJD(1) + t124 * t207 + t156 * t219;
t336 = t219 * t207;
t340 = t219 * qJD(1) + t168 * t156 - t331 * t336;
t284 = t154 * mrSges(5,3);
t339 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t130 + mrSges(6,2) * t131 - t284;
t338 = Ifges(5,5) * qJD(3);
t337 = Ifges(5,6) * qJD(3);
t165 = t196 * t199 - t228;
t112 = t168 * t165;
t263 = qJ(4) - t202;
t189 = t199 * pkin(3) + qJ(2);
t120 = -pkin(4) * t207 + pkin(8) * t165 + t189;
t173 = t263 * t199;
t174 = t263 * t296;
t126 = -t173 * t271 - t196 * t174;
t121 = t201 * t126;
t61 = t198 * t120 + t121;
t103 = t219 * t155;
t123 = t331 * t219;
t335 = t103 - t123;
t102 = t168 * t155;
t334 = t102 - t124;
t65 = -mrSges(6,1) * t144 - mrSges(6,3) * t85;
t66 = mrSges(6,2) * t144 + mrSges(6,3) * t86;
t333 = -t198 * t65 + t201 * t66;
t332 = t18 * t201 - t19 * t198;
t157 = t207 * qJD(3);
t69 = t128 * t196 - t271 * t204;
t279 = t165 * t69;
t330 = t156 * t92 - t157 * t91 + t207 * t70 - t279;
t329 = -m(6) * t83 - t339;
t327 = t19 * mrSges(6,1) - t18 * mrSges(6,2);
t297 = t198 / 0.2e1;
t289 = Ifges(6,4) * t131;
t63 = t130 * Ifges(6,2) + t151 * Ifges(6,6) + t289;
t129 = Ifges(6,4) * t130;
t64 = t131 * Ifges(6,1) + t151 * Ifges(6,5) + t129;
t326 = -t201 * t64 / 0.2e1 + t63 * t297 - t338 / 0.2e1 - t172 * mrSges(5,2);
t325 = -t172 * mrSges(5,1) - t49 * mrSges(6,1) - t10 * mrSges(7,1) + t50 * mrSges(6,2) + t11 * mrSges(7,2);
t324 = qJD(1) ^ 2;
t323 = Ifges(7,4) * t321 + Ifges(7,2) * t320 + Ifges(7,6) * t304;
t322 = Ifges(7,1) * t321 + Ifges(7,4) * t320 + Ifges(7,5) * t304;
t33 = Ifges(7,2) * t234 + Ifges(7,6) * t143 + t295;
t319 = -t33 / 0.2e1;
t318 = t33 / 0.2e1;
t34 = Ifges(7,1) * t75 + Ifges(7,5) * t143 + t71;
t317 = -t34 / 0.2e1;
t316 = t34 / 0.2e1;
t315 = -t234 / 0.2e1;
t314 = t234 / 0.2e1;
t312 = t75 / 0.2e1;
t311 = t85 / 0.2e1;
t310 = t86 / 0.2e1;
t309 = -t130 / 0.2e1;
t308 = -t131 / 0.2e1;
t307 = t131 / 0.2e1;
t305 = t143 / 0.2e1;
t303 = -t151 / 0.2e1;
t302 = -t154 / 0.2e1;
t301 = t154 / 0.2e1;
t300 = -t155 / 0.2e1;
t290 = Ifges(4,4) * t199;
t288 = Ifges(6,4) * t198;
t287 = Ifges(6,4) * t201;
t125 = -t173 * t196 + t271 * t174;
t286 = t125 * t69;
t285 = t144 * mrSges(5,3);
t283 = t154 * Ifges(5,4);
t282 = t155 * mrSges(5,3);
t267 = t157 * t201;
t266 = t165 * t198;
t265 = t165 * t201;
t177 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t245;
t264 = t199 * t177;
t180 = pkin(3) * t244 + qJD(2);
t44 = -mrSges(7,1) * t234 + mrSges(7,2) * t75;
t256 = t44 + t339;
t254 = Ifges(6,5) * t85 + Ifges(6,6) * t86 - Ifges(6,3) * t144;
t251 = Ifges(4,4) * t296;
t248 = t271 * pkin(3);
t240 = t259 / 0.2e1;
t106 = -pkin(4) * t156 - pkin(8) * t157 + t180;
t146 = t261 * t263 - t243;
t147 = -qJD(3) * t174 - t258;
t99 = t196 * t146 + t147 * t271;
t239 = t201 * t106 - t198 * t99;
t235 = -t144 * mrSges(5,1) + t145 * mrSges(5,2);
t60 = t201 * t120 - t126 * t198;
t98 = -t271 * t146 + t147 * t196;
t188 = -t248 - pkin(4);
t227 = mrSges(6,1) * t198 + mrSges(6,2) * t201;
t226 = Ifges(6,1) * t201 - t288;
t225 = -Ifges(6,2) * t198 + t287;
t224 = Ifges(6,5) * t201 - Ifges(6,6) * t198;
t47 = -pkin(5) * t207 + pkin(9) * t265 + t60;
t53 = pkin(9) * t266 + t61;
t21 = -t197 * t53 + t200 * t47;
t22 = t197 * t47 + t200 * t53;
t223 = t50 * t198 + t49 * t201;
t222 = t198 * t49 - t201 * t50;
t93 = -mrSges(6,2) * t151 + mrSges(6,3) * t130;
t94 = mrSges(6,1) * t151 - mrSges(6,3) * t131;
t221 = -t198 * t94 + t201 * t93;
t220 = -t198 * t93 - t201 * t94;
t134 = -qJD(3) * mrSges(5,2) + t282;
t216 = -t134 - t221;
t215 = -Ifges(4,5) * t199 - Ifges(4,6) * t296;
t214 = t83 * t227;
t213 = -t157 * t198 + t165 * t259;
t212 = t165 * t260 + t267;
t25 = t198 * t106 + t120 * t259 - t126 * t260 + t201 * t99;
t211 = qJ(2) * (mrSges(4,1) * t296 - mrSges(4,2) * t199);
t210 = t199 * (-Ifges(4,2) * t296 - t290);
t209 = (Ifges(4,1) * t296 - t290) * qJD(1);
t208 = (-Ifges(4,2) * t199 + t251) * qJD(1);
t169 = (t199 * mrSges(4,1) + mrSges(4,2) * t296) * qJD(1);
t206 = (-Ifges(4,1) * t199 - t251) * t296;
t205 = -qJD(5) * t223 + t332;
t176 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t262;
t175 = -t201 * pkin(5) + t188;
t159 = Ifges(4,5) * qJD(3) + t209;
t158 = Ifges(4,6) * qJD(3) + t208;
t150 = Ifges(5,4) * t155;
t117 = -mrSges(5,1) * t155 - mrSges(5,2) * t154;
t114 = t219 * t165;
t111 = t168 * t207;
t110 = -t154 * Ifges(5,1) + t150 + t338;
t109 = t155 * Ifges(5,2) - t283 + t337;
t89 = -pkin(5) * t266 + t125;
t56 = mrSges(7,1) * t143 - mrSges(7,3) * t75;
t55 = -mrSges(7,2) * t143 + mrSges(7,3) * t234;
t54 = -pkin(5) * t213 + t98;
t45 = -pkin(5) * t86 + t69;
t43 = t85 * Ifges(6,1) + t86 * Ifges(6,4) - t144 * Ifges(6,5);
t42 = t85 * Ifges(6,4) + t86 * Ifges(6,2) - t144 * Ifges(6,6);
t41 = -t123 * t165 - t157 * t168;
t39 = t112 * t331 - t219 * t157;
t26 = -qJD(5) * t61 + t239;
t24 = mrSges(7,2) * t144 + mrSges(7,3) * t30;
t23 = -mrSges(7,1) * t144 - mrSges(7,3) * t29;
t20 = pkin(9) * t213 + t25;
t17 = -pkin(9) * t267 - pkin(5) * t156 + (-t121 + (-pkin(9) * t165 - t120) * t198) * qJD(5) + t239;
t14 = t200 * t36 - t275;
t13 = -t197 * t36 - t273;
t8 = -mrSges(7,1) * t30 + mrSges(7,2) * t29;
t5 = -qJD(6) * t22 + t17 * t200 - t197 * t20;
t4 = qJD(6) * t21 + t17 * t197 + t20 * t200;
t1 = [t130 * (Ifges(6,4) * t212 + Ifges(6,2) * t213) / 0.2e1 + (t176 * t244 - t177 * t261) * t202 - (t255 + t254) * t207 / 0.2e1 + (t155 * t157 / 0.2e1 - t165 * t144 + t207 * t145) * Ifges(5,4) + (t144 * t207 + t155 * t353) * Ifges(5,2) - (t171 * mrSges(5,1) + Ifges(6,5) * t311 + Ifges(7,5) * t321 + Ifges(6,6) * t310 + Ifges(7,6) * t320 + (Ifges(6,3) + Ifges(7,3)) * t304 + t327 + t328) * t207 + m(5) * (t126 * t70 + t171 * t189 + t172 * t180 + t92 * t99 + t286) - (-t110 / 0.2e1 - Ifges(5,1) * t302 + t326) * t157 - (t208 + t158) * t244 / 0.2e1 - (t209 + t159) * t261 / 0.2e1 + (qJD(5) * t64 + t42) * t266 / 0.2e1 + t330 * mrSges(5,3) + (Ifges(7,5) * t114 + Ifges(7,6) * t112) * t304 + (Ifges(7,5) * t39 + Ifges(7,6) * t41) * t305 + t83 * (-mrSges(6,1) * t213 + mrSges(6,2) * t212) + t151 * (Ifges(6,5) * t212 + Ifges(6,6) * t213) / 0.2e1 + (-t350 / 0.2e1 + t302 * Ifges(5,4) - t344 / 0.2e1 - t345 / 0.2e1 - Ifges(7,3) * t305 - Ifges(6,5) * t307 - Ifges(7,5) * t312 - Ifges(7,6) * t314 + t109 / 0.2e1 + t325) * t156 + m(7) * (t10 * t5 + t11 * t4 + t2 * t22 + t21 * t3 + t45 * t89 + t54 * t59) + (-t10 * t39 + t11 * t41 + t112 * t2 - t114 * t3) * mrSges(7,3) + (t18 * t266 + t19 * t265 - t212 * t49 + t213 * t50) * mrSges(6,3) + t342 * t125 + 0.2e1 * t356 * t194 + m(6) * (t18 * t61 + t19 * t60 + t25 * t50 + t26 * t49 + t286) + (0.2e1 * t211 - t210 + t206) * t257 + (Ifges(7,1) * t39 + Ifges(7,4) * t41) * t312 + (Ifges(7,1) * t114 + Ifges(7,4) * t112) * t321 + (Ifges(7,4) * t39 + Ifges(7,2) * t41) * t314 + (Ifges(7,4) * t114 + Ifges(7,2) * t112) * t320 + t22 * t24 + t21 * t23 + (Ifges(5,6) * t353 + t215 * qJD(3) / 0.2e1) * qJD(3) + (-m(5) * t91 - t329) * t98 - t227 * t279 - t43 * t265 / 0.2e1 + (Ifges(6,1) * t212 + Ifges(6,4) * t213) * t307 + t54 * t44 + t4 * t55 + t5 * t56 + t59 * (-mrSges(7,1) * t41 + mrSges(7,2) * t39) + t60 * t65 + t61 * t66 + t89 * t8 + t25 * t93 + t26 * t94 + t45 * (-mrSges(7,1) * t112 + mrSges(7,2) * t114) + t126 * t285 + t39 * t316 + t41 * t318 + t114 * t322 + t112 * t323 + (-t171 * mrSges(5,2) - Ifges(5,1) * t145 - t224 * t304 - t225 * t310 - t226 * t311 + t240 * t63) * t165 + t99 * t134 + 0.2e1 * qJD(2) * t169 + t180 * t117 + t189 * t235; t111 * t23 + t336 * t24 + t340 * t56 + t341 * t55 + (t176 * t296 - t264) * qJD(3) - t356 * t324 + (t8 + t342) * t165 - t256 * t157 + t216 * t156 - (qJD(5) * t220 + t285 + t333) * t207 + m(6) * (t156 * t222 - t157 * t83 - t205 * t207 + t279) - m(5) * t330 + (-m(5) * t172 - m(6) * t223 - t117 - t169 + t220) * qJD(1) + (t10 * t340 + t11 * t341 + t111 * t3 - t157 * t59 + t165 * t45 + t2 * t336) * m(7); (t150 + t110) * t300 + (t283 + t350) * t301 + t351 * t55 + t352 * t56 + (t10 * t352 + t11 * t351 + t118 * t3 + t119 * t2 + t175 * t45 + t343 * t59) * m(7) + (m(6) * t188 - mrSges(6,1) * t201 + mrSges(6,2) * t198 - mrSges(5,1)) * t69 + (-Ifges(7,4) * t103 - Ifges(7,2) * t102) * t315 + (-Ifges(7,5) * t103 - Ifges(7,6) * t102) * t306 + ((t196 * t70 - t271 * t69) * pkin(3) + t91 * t104 - t92 * t105 - t172 * t232) * m(5) + (Ifges(6,5) * t198 + Ifges(7,5) * t168 + Ifges(6,6) * t201 - Ifges(7,6) * t219) * t304 + (Ifges(7,4) * t168 - Ifges(7,2) * t219) * t320 + (Ifges(7,1) * t168 - Ifges(7,4) * t219) * t321 + t45 * (mrSges(7,1) * t219 + mrSges(7,2) * t168) + (-t10 * t335 + t11 * t334 - t168 * t3 - t2 * t219) * mrSges(7,3) - t219 * t323 + (-Ifges(7,5) * t123 - Ifges(7,6) * t124) * t305 + (-Ifges(7,1) * t123 - Ifges(7,4) * t124) * t312 + (-Ifges(7,4) * t123 - Ifges(7,2) * t124) * t314 + (Ifges(5,1) * t301 + t224 * t303 + t225 * t309 + t226 * t308 - t214 + t326) * t155 + (-Ifges(7,1) * t103 - Ifges(7,4) * t102) * t313 + (t130 * t225 + t131 * t226 + t151 * t224) * qJD(5) / 0.2e1 + (m(6) * t205 - t259 * t94 - t260 * t93 + t333) * t187 + (-mrSges(7,1) * t334 + mrSges(7,2) * t335) * t59 + (-t359 * t50 + (-t259 + t268) * t49 + t332) * mrSges(6,3) - m(6) * (t49 * t51 + t50 * t52) + t343 * t44 + (Ifges(5,2) * t300 - Ifges(6,3) * t303 - Ifges(7,3) * t306 - Ifges(6,5) * t308 - Ifges(6,6) * t309 - Ifges(7,5) * t313 - Ifges(7,6) * t315 - t337 / 0.2e1 - t325) * t154 + (t210 / 0.2e1 - t211 - t206 / 0.2e1) * t324 + t329 * t104 + t64 * t240 - t248 * t291 - t92 * t284 - t63 * t260 / 0.2e1 + t159 * t262 / 0.2e1 - t215 * t257 / 0.2e1 - t176 * t170 - mrSges(4,1) * t247 + t158 * t245 / 0.2e1 + qJD(5) * t214 - t70 * mrSges(5,2) - t52 * t93 - t51 * t94 + t179 * t264 + t91 * t282 + t285 * t294 + t43 * t297 + t109 * t302 + (Ifges(6,2) * t201 + t288) * t310 + (Ifges(6,1) * t198 + t287) * t311 - t123 * t316 - t103 * t317 - t124 * t318 - t102 * t319 + t168 * t322 + t118 * t23 + t119 * t24 - t105 * t134 + Ifges(5,6) * t144 + Ifges(5,5) * t145 - Ifges(4,6) * t230 - mrSges(4,2) * t231 + t175 * t8 + t188 * t48 + t201 * t42 / 0.2e1 - t117 * t232 - Ifges(4,5) * t242; -t219 * t23 + t168 * t24 + t198 * t66 + t201 * t65 + t334 * t56 + t335 * t55 + t221 * qJD(5) + t216 * t155 + t256 * t154 + t235 + (t10 * t334 + t11 * t335 + t154 * t59 + t168 * t2 - t219 * t3) * m(7) + (-t151 * t222 + t154 * t83 + t18 * t198 + t19 * t201) * m(6) + (-t154 * t91 - t155 * t92 + t171) * m(5); t327 - t75 * t319 - m(7) * (t10 * t13 + t11 * t14) + t234 * t317 + (t130 * t49 + t131 * t50) * mrSges(6,3) + t254 + (-t131 * t44 + t197 * t24 + t200 * t23 + (-t197 * t56 + t200 * t55) * qJD(6) + (-t131 * t59 + t197 * t2 + t200 * t3 + (-t10 * t197 + t11 * t200) * qJD(6)) * m(7)) * pkin(5) + (-Ifges(6,2) * t131 + t129 + t64) * t309 - t14 * t55 - t13 * t56 - t49 * t93 + t50 * t94 + (Ifges(6,5) * t130 - Ifges(6,6) * t131) * t303 + t63 * t307 + (Ifges(6,1) * t130 - t289) * t308 + t357 * t315 - t83 * (mrSges(6,1) * t131 + mrSges(6,2) * t130) + t358; t33 * t312 - t10 * t55 + t11 * t56 + (t34 + t357) * t315 + t358;];
tauc  = t1(:);
