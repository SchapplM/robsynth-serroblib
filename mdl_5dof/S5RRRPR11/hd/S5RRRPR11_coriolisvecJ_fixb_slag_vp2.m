% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:26
% EndTime: 2019-12-31 21:32:45
% DurationCPUTime: 8.15s
% Computational Cost: add. (4153->533), mult. (10105->711), div. (0->0), fcn. (6102->6), ass. (0->272)
t354 = qJD(2) / 0.2e1;
t233 = Ifges(3,5) * t354;
t350 = Ifges(4,1) + Ifges(5,1);
t343 = Ifges(5,4) + Ifges(4,5);
t353 = Ifges(4,6) - Ifges(5,6);
t179 = sin(qJ(2));
t182 = cos(qJ(2));
t225 = pkin(2) * t179 - pkin(7) * t182;
t146 = t225 * qJD(1);
t178 = sin(qJ(3));
t129 = t178 * t146;
t251 = qJD(1) * t179;
t170 = qJ(4) * t251;
t247 = qJD(3) * t178;
t181 = cos(qJ(3));
t256 = t179 * t181;
t257 = t178 * t182;
t304 = pkin(7) - pkin(8);
t352 = t304 * t247 + t129 + t170 + (-pkin(6) * t256 + pkin(8) * t257) * qJD(1);
t159 = t304 * t181;
t285 = pkin(6) * t178;
t237 = -pkin(3) - t285;
t197 = (-pkin(4) + t237) * t179;
t255 = t181 * t182;
t189 = -pkin(8) * t255 + t197;
t258 = t146 * t181;
t351 = -qJD(1) * t189 + qJD(3) * t159 + t258;
t153 = -pkin(2) * t182 - pkin(7) * t179 - pkin(1);
t131 = t153 * qJD(1);
t250 = qJD(1) * t182;
t174 = pkin(6) * t250;
t157 = qJD(2) * pkin(7) + t174;
t82 = t181 * t131 - t178 * t157;
t320 = qJD(4) - t82;
t167 = qJD(3) - t250;
t236 = t178 * t251;
t243 = t181 * qJD(2);
t140 = t236 - t243;
t235 = t181 * t251;
t141 = qJD(2) * t178 + t235;
t177 = sin(qJ(5));
t180 = cos(qJ(5));
t200 = t140 * t177 + t141 * t180;
t78 = t140 * t180 - t141 * t177;
t72 = Ifges(6,4) * t78;
t349 = Ifges(6,2) * t200 - t72;
t348 = -pkin(8) * t141 + t320;
t347 = qJD(4) * t178 + t174;
t172 = Ifges(3,4) * t250;
t156 = -qJD(2) * pkin(2) + pkin(6) * t251;
t83 = t178 * t131 + t181 * t157;
t201 = t178 * t83 + t181 * t82;
t61 = -pkin(3) * t167 + t320;
t161 = t167 * qJ(4);
t62 = t161 + t83;
t202 = t178 * t62 - t181 * t61;
t270 = Ifges(5,5) * t181;
t209 = Ifges(5,3) * t178 + t270;
t275 = Ifges(4,4) * t181;
t213 = -Ifges(4,2) * t178 + t275;
t218 = mrSges(5,1) * t178 - mrSges(5,3) * t181;
t220 = mrSges(4,1) * t178 + mrSges(4,2) * t181;
t268 = Ifges(5,6) * t178;
t269 = Ifges(4,6) * t178;
t273 = Ifges(4,5) * t181;
t274 = Ifges(5,4) * t181;
t287 = t181 / 0.2e1;
t289 = t178 / 0.2e1;
t290 = -t178 / 0.2e1;
t295 = t141 / 0.2e1;
t297 = t140 / 0.2e1;
t298 = -t140 / 0.2e1;
t271 = Ifges(5,5) * t178;
t276 = Ifges(4,4) * t178;
t319 = t181 * t350 + t271 - t276;
t137 = Ifges(4,4) * t140;
t272 = Ifges(5,5) * t140;
t326 = t141 * t350 + t343 * t167 - t137 + t272;
t330 = t167 / 0.2e1;
t194 = qJ(4) * t141 - t156;
t64 = pkin(3) * t140 - t194;
t136 = Ifges(5,5) * t141;
t65 = Ifges(5,6) * t167 + Ifges(5,3) * t140 + t136;
t277 = Ifges(4,4) * t141;
t68 = -Ifges(4,2) * t140 + Ifges(4,6) * t167 + t277;
t315 = t202 * mrSges(5,2) + t201 * mrSges(4,3) - t156 * t220 - t209 * t297 - t213 * t298 - t64 * t218 - t65 * t289 - t68 * t290 - t319 * t295 - (t268 + t274 - t269 + t273) * t330 - t326 * t287;
t346 = Ifges(3,1) * t251 / 0.2e1 + t172 / 0.2e1 + t233 - t315;
t242 = qJD(1) * qJD(2);
t231 = t179 * t242;
t162 = qJD(5) - t167;
t294 = -t162 / 0.2e1;
t234 = t179 * t247;
t241 = qJD(2) * qJD(3);
t102 = t181 * t241 + (t182 * t243 - t234) * qJD(1);
t246 = qJD(3) * t181;
t248 = qJD(2) * t182;
t103 = t178 * t241 + (t178 * t248 + t179 * t246) * qJD(1);
t19 = qJD(5) * t78 + t102 * t180 + t103 * t177;
t20 = -qJD(5) * t200 - t102 * t177 + t103 * t180;
t327 = Ifges(6,5) * t19 + Ifges(6,6) * t20;
t305 = pkin(3) + pkin(4);
t51 = -t140 * t305 + t194;
t34 = -t167 * t305 + t348;
t56 = pkin(8) * t140 + t83;
t47 = t161 + t56;
t8 = -t177 * t47 + t180 * t34;
t9 = t177 * t34 + t180 * t47;
t345 = (t200 * t9 + t78 * t8) * mrSges(6,3) - Ifges(6,3) * t231 + (Ifges(6,5) * t78 - Ifges(6,6) * t200) * t294 - t51 * (mrSges(6,1) * t200 + mrSges(6,2) * t78) + t327;
t344 = -qJD(2) / 0.2e1;
t342 = t343 * t231 + (-Ifges(4,4) + Ifges(5,5)) * t103 + t350 * t102;
t158 = t304 * t178;
t98 = t158 * t180 - t159 * t177;
t341 = qJD(5) * t98 + t351 * t177 - t352 * t180;
t99 = t158 * t177 + t159 * t180;
t340 = -qJD(5) * t99 + t352 * t177 + t351 * t180;
t260 = qJ(4) * t181;
t193 = -t178 * t305 + t260;
t339 = t167 * t193 + t347;
t206 = pkin(3) * t178 - t260;
t338 = t167 * t206 - t347;
t337 = t343 * t178 + t353 * t181;
t336 = t350 * t178 - t270 + t275;
t286 = Ifges(6,4) * t200;
t333 = Ifges(6,1) * t78 - t286;
t311 = t19 / 0.2e1;
t310 = t20 / 0.2e1;
t25 = Ifges(6,1) * t200 + Ifges(6,5) * t162 + t72;
t331 = t25 / 0.2e1;
t307 = -t200 / 0.2e1;
t329 = -t231 / 0.2e1;
t232 = Ifges(3,6) * t344;
t151 = t180 * qJ(4) - t177 * t305;
t322 = -qJD(5) * t151 - t177 * t348 - t180 * t56;
t150 = -t177 * qJ(4) - t180 * t305;
t321 = qJD(5) * t150 - t177 * t56 + t180 * t348;
t199 = t177 * t181 - t178 * t180;
t124 = t199 * t179;
t318 = qJD(3) - qJD(5);
t261 = qJ(4) * t178;
t316 = -t181 * t305 - t261;
t238 = Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1;
t239 = Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1;
t240 = Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t278 = Ifges(3,4) * t179;
t314 = t238 * t140 + t240 * t141 + t239 * t167 + t62 * mrSges(5,3) + t82 * mrSges(4,1) + t9 * mrSges(6,2) + t232 - (t182 * Ifges(3,2) + t278) * qJD(1) / 0.2e1 - t162 * Ifges(6,3) - t200 * Ifges(6,5) - t78 * Ifges(6,6) + Ifges(4,6) * t298 + Ifges(5,6) * t297 - t61 * mrSges(5,1) - t8 * mrSges(6,1) - t83 * mrSges(4,2) + (Ifges(4,3) + Ifges(5,2)) * t330 + t343 * t295;
t313 = Ifges(6,4) * t311 + Ifges(6,2) * t310 + Ifges(6,6) * t329;
t312 = Ifges(6,1) * t311 + Ifges(6,4) * t310 + Ifges(6,5) * t329;
t309 = -t78 / 0.2e1;
t308 = t78 / 0.2e1;
t306 = t200 / 0.2e1;
t303 = pkin(1) * mrSges(3,1);
t302 = pkin(1) * mrSges(3,2);
t301 = t102 / 0.2e1;
t300 = -t103 / 0.2e1;
t299 = t103 / 0.2e1;
t296 = -t141 / 0.2e1;
t293 = t162 / 0.2e1;
t292 = -t167 / 0.2e1;
t288 = -t181 / 0.2e1;
t284 = pkin(8) * t179;
t280 = mrSges(4,3) * t140;
t279 = mrSges(4,3) * t141;
t192 = t182 * t199;
t114 = qJD(1) * t192;
t86 = t318 * t199;
t266 = -t114 + t86;
t198 = t177 * t178 + t180 * t181;
t191 = t182 * t198;
t115 = qJD(1) * t191;
t85 = t318 * t198;
t265 = t115 - t85;
t262 = qJ(4) * t140;
t259 = qJD(2) * mrSges(3,2);
t106 = -mrSges(4,2) * t167 - t280;
t109 = -mrSges(5,2) * t140 + mrSges(5,3) * t167;
t254 = -t106 - t109;
t107 = mrSges(4,1) * t167 - t279;
t108 = -mrSges(5,1) * t167 + mrSges(5,2) * t141;
t253 = -t107 + t108;
t149 = t225 * qJD(2);
t252 = t178 * t149 + t153 * t246;
t169 = pkin(6) * t255;
t113 = t178 * t153 + t169;
t249 = qJD(2) * t179;
t244 = qJD(4) * t181;
t168 = pkin(6) * t257;
t112 = t153 * t181 - t168;
t230 = pkin(6) * t231;
t228 = m(4) * t156 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t140 + mrSges(4,2) * t141 + mrSges(3,3) * t251;
t132 = qJD(1) * t149;
t36 = t131 * t246 + t178 * t132 - t157 * t247 - t181 * t230;
t26 = qJ(4) * t231 + t167 * qJD(4) + t36;
t10 = pkin(8) * t103 + t26;
t223 = t131 * t247 - t132 * t181 + t157 * t246;
t11 = -pkin(8) * t102 + t197 * t242 + t223;
t1 = qJD(5) * t8 + t10 * t180 + t11 * t177;
t2 = -qJD(5) * t9 - t10 * t177 + t11 * t180;
t227 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t226 = t237 * t179;
t100 = -qJ(4) * t182 + t113;
t222 = qJD(3) * t169 - t149 * t181 + t153 * t247;
t221 = mrSges(4,1) * t181 - mrSges(4,2) * t178;
t219 = mrSges(5,1) * t181 + mrSges(5,3) * t178;
t212 = Ifges(4,2) * t181 + t276;
t208 = -Ifges(5,3) * t181 + t271;
t207 = pkin(3) * t181 + t261;
t57 = -mrSges(6,2) * t162 + mrSges(6,3) * t78;
t58 = mrSges(6,1) * t162 - mrSges(6,3) * t200;
t205 = -t177 * t58 + t180 * t57;
t176 = t182 * pkin(3);
t71 = pkin(4) * t182 + t168 + t176 + (-t153 - t284) * t181;
t81 = t178 * t284 + t100;
t27 = -t177 * t81 + t180 * t71;
t28 = t177 * t71 + t180 * t81;
t195 = qJD(2) * t226;
t33 = qJD(1) * t195 + t223;
t204 = t178 * t33 + t181 * t26;
t37 = t178 * t230 - t223;
t203 = -t178 * t37 + t181 * t36;
t105 = -pkin(6) * t235 + t129;
t75 = -mrSges(5,1) * t231 + t102 * mrSges(5,2);
t196 = pkin(6) + t206;
t190 = -pkin(6) + t193;
t188 = qJ(4) * t102 - qJD(2) * t174 + qJD(4) * t141;
t53 = (-t179 * t243 - t182 * t247) * pkin(6) + t252;
t187 = -t37 * mrSges(4,1) + t33 * mrSges(5,1) + t36 * mrSges(4,2) - t26 * mrSges(5,3) + t227;
t171 = qJ(4) * t249;
t166 = Ifges(5,2) * t231;
t165 = Ifges(4,3) * t231;
t155 = mrSges(3,3) * t250 - t259;
t152 = -pkin(2) - t207;
t135 = pkin(2) - t316;
t125 = t198 * t179;
t118 = t196 * t179;
t104 = pkin(6) * t236 + t258;
t101 = -t112 + t176;
t97 = t190 * t179;
t96 = Ifges(5,4) * t102;
t95 = Ifges(4,5) * t102;
t94 = Ifges(4,6) * t103;
t93 = Ifges(5,6) * t103;
t90 = mrSges(5,1) * t140 - mrSges(5,3) * t141;
t89 = pkin(3) * t141 + t262;
t88 = qJD(1) * t226 - t258;
t87 = t105 + t170;
t76 = -mrSges(4,2) * t231 - mrSges(4,3) * t103;
t74 = mrSges(4,1) * t231 - mrSges(4,3) * t102;
t73 = -mrSges(5,2) * t103 + mrSges(5,3) * t231;
t59 = -t141 * t305 - t262;
t54 = t249 * t285 - t222;
t52 = (qJD(3) * t207 - t244) * t179 + t196 * t248;
t50 = t195 + t222;
t49 = mrSges(4,1) * t103 + mrSges(4,2) * t102;
t48 = mrSges(5,1) * t103 - mrSges(5,3) * t102;
t46 = -qJD(4) * t182 + t171 + t53;
t41 = t102 * Ifges(4,4) - t103 * Ifges(4,2) + Ifges(4,6) * t231;
t40 = t102 * Ifges(5,5) + Ifges(5,6) * t231 + t103 * Ifges(5,3);
t39 = qJD(2) * t191 + t124 * t318;
t38 = -qJD(2) * t192 + t179 * t85;
t35 = (qJD(3) * t316 + t244) * t179 + t190 * t248;
t32 = -mrSges(6,1) * t78 + mrSges(6,2) * t200;
t31 = t171 + (-pkin(6) * qJD(2) + pkin(8) * qJD(3)) * t256 + (-qJD(4) + (-pkin(6) * qJD(3) + pkin(8) * qJD(2)) * t178) * t182 + t252;
t30 = pkin(3) * t103 - t188;
t29 = pkin(8) * t234 + qJD(2) * t189 + t222;
t24 = Ifges(6,2) * t78 + Ifges(6,6) * t162 + t286;
t14 = mrSges(6,2) * t231 + mrSges(6,3) * t20;
t13 = -mrSges(6,1) * t231 - mrSges(6,3) * t19;
t12 = -t103 * t305 + t188;
t7 = -mrSges(6,1) * t20 + mrSges(6,2) * t19;
t4 = -qJD(5) * t28 - t177 * t31 + t180 * t29;
t3 = qJD(5) * t27 + t177 * t29 + t180 * t31;
t5 = [t52 * t90 + t51 * (-mrSges(6,1) * t38 + mrSges(6,2) * t39) + t3 * t57 + t4 * t58 + t35 * t32 + t38 * t24 / 0.2e1 + t27 * t13 + t28 * t14 + m(4) * (t112 * t37 + t113 * t36 + t53 * t83 + t54 * t82) + (Ifges(6,4) * t125 - Ifges(6,2) * t124) * t310 + (Ifges(6,1) * t125 - Ifges(6,4) * t124) * t311 + t12 * (mrSges(6,1) * t124 + mrSges(6,2) * t125) + (-t1 * t124 - t125 * t2 + t38 * t9 - t39 * t8) * mrSges(6,3) + m(5) * (t100 * t26 + t101 * t33 + t118 * t30 + t46 * t62 + t50 * t61 + t52 * t64) + m(6) * (t1 * t28 + t12 * t97 + t2 * t27 + t3 * t9 + t35 * t51 + t4 * t8) + (Ifges(6,1) * t39 + Ifges(6,4) * t38) * t306 + (Ifges(6,4) * t39 + Ifges(6,2) * t38) * t308 + t125 * t312 - t124 * t313 + (Ifges(6,5) * t39 + Ifges(6,6) * t38) * t293 + (((-0.2e1 * t302 + 0.3e1 / 0.2e1 * Ifges(3,4) * t182) * qJD(1) + t228 * pkin(6) + t233 + t346) * qJD(2) + t187 - t238 * t103 - t240 * t102 - t165 / 0.2e1 - t166 / 0.2e1 - t93 / 0.2e1 + t94 / 0.2e1 - t95 / 0.2e1 - t96 / 0.2e1 + t327) * t182 + t97 * t7 + t100 * t73 + t101 * t75 + t53 * t106 + t54 * t107 + t50 * t108 + t46 * t109 + t112 * t74 + t113 * t76 + t118 * t48 + (t209 * t299 + t213 * t300 + pkin(6) * t49 + t40 * t289 + t30 * t218 + (-t178 * t36 - t181 * t37) * mrSges(4,3) + (-t178 * t26 + t181 * t33) * mrSges(5,2) + (t64 * t219 + t208 * t298 + t212 * t297 + t156 * t221 + t68 * t288 + (t178 * t82 - t181 * t83) * mrSges(4,3) + (-t178 * t61 - t181 * t62) * mrSges(5,2) + t336 * t296 + t337 * t292) * qJD(3) + ((-0.2e1 * t303 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t274 / 0.2e1 + t268 / 0.2e1 + t273 / 0.2e1 - t269 / 0.2e1) * t179 - Ifges(6,5) * t125 / 0.2e1 + Ifges(6,6) * t124 / 0.2e1 + (-0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(6,3) + 0.3e1 / 0.2e1 * Ifges(3,1) + (m(4) * pkin(6) + t220) * pkin(6) - t239) * t182) * qJD(1) - pkin(6) * t155 + t232 + t314) * qJD(2) + t319 * t301 + (qJD(3) * t326 + t41) * t290 + (qJD(3) * t65 + t342) * t287) * t179 + t39 * t331; -pkin(2) * t49 + (t85 / 0.2e1 - t115 / 0.2e1) * t25 - t315 * qJD(3) - t30 * t219 + (mrSges(6,1) * t266 - mrSges(6,2) * t265) * t51 - m(4) * (t104 * t82 + t105 * t83) + t208 * t299 + t212 * t300 + t41 * t287 + t40 * t288 - t198 * t313 - t199 * t312 + (-t1 * t198 + t199 * t2 + t265 * t8 - t266 * t9) * mrSges(6,3) + (-Ifges(6,4) * t199 - Ifges(6,2) * t198) * t310 + (-Ifges(6,1) * t199 - Ifges(6,4) * t198) * t311 + t12 * (mrSges(6,1) * t198 - mrSges(6,2) * t199) + (Ifges(6,4) * t85 - Ifges(6,2) * t86) * t308 + (Ifges(6,5) * t85 - Ifges(6,6) * t86) * t293 + (Ifges(6,1) * t85 - Ifges(6,4) * t86) * t306 + (Ifges(6,1) * t115 - Ifges(6,4) * t114) * t307 + (Ifges(6,4) * t115 - Ifges(6,2) * t114) * t309 + (Ifges(6,5) * t115 - Ifges(6,6) * t114) * t294 + (-t86 / 0.2e1 + t114 / 0.2e1) * t24 + t203 * mrSges(4,3) + t204 * mrSges(5,2) + t98 * t13 + t99 * t14 - t105 * t106 - t104 * t107 - t88 * t108 - t87 * t109 + t336 * t301 + t338 * t90 + (t152 * t30 + t338 * t64 - t61 * t88 - t62 * t87) * m(5) + t135 * t7 + t339 * t32 + t152 * t48 + t340 * t58 + t341 * t57 + (t1 * t99 + t12 * t135 + t2 * t98 + t339 * t51 + t340 * t8 + t341 * t9) * m(6) + t342 * t289 + (((-Ifges(6,5) * t199 - Ifges(6,6) * t198) * t344 + (t278 / 0.2e1 + t303) * qJD(1) + (t155 + t259) * pkin(6) + t232 + t337 * t354 - t314) * t179 + ((t302 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t179) * qJD(1) + ((-m(4) * pkin(2) - mrSges(3,1) - t221) * qJD(2) - t228) * pkin(6) - t172 / 0.2e1 + t233 - t346) * t182) * qJD(1) + ((-m(4) * t201 - m(5) * t202 + t178 * t254 + t181 * t253) * qJD(3) + (t73 + t76) * t181 + (-t74 + t75) * t178 + m(4) * t203 + m(5) * t204) * pkin(7); -t89 * t90 + qJ(4) * t73 - pkin(3) * t75 - t59 * t32 + t321 * t57 + (t1 * t151 + t150 * t2 + t321 * t9 + t322 * t8 - t51 * t59) * m(6) + t322 * t58 + (-Ifges(4,2) * t141 - t137 + t326) * t297 + t68 * t295 + (Ifges(5,3) * t141 - t272) * t298 + (-t343 * t140 - t353 * t141) * t292 + (-t140 * t350 + t136 - t277 + t65) * t296 - t187 + (-t253 + t279) * t83 + (t254 - t280) * t82 - t345 + (-pkin(3) * t33 + qJ(4) * t26 + t320 * t62 - t61 * t83 - t64 * t89) * m(5) + t78 * t331 + (t140 * t61 + t141 * t62) * mrSges(5,2) + t349 * t309 + t165 + t166 + (t24 - t333) * t307 + qJD(4) * t109 - t64 * (mrSges(5,1) * t141 + mrSges(5,3) * t140) + t150 * t13 + t151 * t14 - t156 * (mrSges(4,1) * t141 - mrSges(4,2) * t140) + t93 - t94 + t95 + t96; t180 * t13 + t177 * t14 + (-t32 + t90) * t141 + t205 * qJD(5) + (-t109 - t205) * t167 + t75 + (t1 * t177 - t141 * t51 + t180 * t2 + t162 * (-t177 * t8 + t180 * t9)) * m(6) + (t141 * t64 - t167 * t62 + t33) * m(5); t333 * t307 + t24 * t306 - t8 * t57 + t9 * t58 + t227 + (t25 - t349) * t309 + t345;];
tauc = t5(:);
