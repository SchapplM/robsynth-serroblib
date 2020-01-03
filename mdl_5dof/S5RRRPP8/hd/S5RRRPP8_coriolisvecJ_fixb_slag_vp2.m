% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:31
% EndTime: 2019-12-31 21:07:45
% DurationCPUTime: 6.45s
% Computational Cost: add. (2402->470), mult. (6157->581), div. (0->0), fcn. (3343->4), ass. (0->230)
t304 = -Ifges(4,4) + Ifges(6,6);
t303 = Ifges(4,1) + Ifges(6,3);
t299 = Ifges(6,4) + Ifges(5,5);
t298 = Ifges(4,5) + Ifges(6,5);
t302 = Ifges(6,2) + Ifges(5,3);
t212 = Ifges(3,5) * qJD(2) / 0.2e1;
t150 = cos(qJ(2));
t226 = qJD(1) * t150;
t140 = pkin(6) * t226;
t147 = sin(qJ(3));
t223 = qJD(3) * t147;
t301 = -qJD(4) * t147 - t140 + (-t147 * t226 + t223) * pkin(3);
t300 = t304 * t147;
t148 = sin(qJ(2));
t220 = qJD(1) * qJD(2);
t210 = t148 * t220;
t149 = cos(qJ(3));
t209 = t150 * t220;
t214 = t148 * t223;
t219 = qJD(2) * qJD(3);
t74 = qJD(1) * t214 + (-t209 - t219) * t149;
t222 = qJD(3) * t149;
t224 = qJD(2) * t150;
t288 = t147 * t224 + t148 * t222;
t75 = t288 * qJD(1) + t147 * t219;
t297 = t302 * t75 + (Ifges(5,6) - Ifges(6,6)) * t74 + t299 * t210;
t296 = t298 * t210 - t303 * t74 + t304 * t75;
t236 = qJ(4) * t149;
t167 = qJ(5) * t147 - t236;
t159 = t167 * t150;
t295 = -qJD(1) * t159 + qJD(3) * t167 - qJD(5) * t149 + t301;
t294 = -qJ(4) * t222 + t226 * t236 + t301;
t207 = pkin(6) * t149 - qJ(4);
t193 = t207 * t148;
t233 = t147 * t150;
t156 = -pkin(4) * t233 - t193;
t271 = pkin(4) + pkin(7);
t197 = pkin(2) * t148 - pkin(7) * t150;
t109 = t197 * qJD(1);
t95 = t147 * t109;
t293 = -qJD(1) * t156 - t271 * t223 - t95;
t120 = t271 * t149;
t257 = pkin(6) * t147;
t218 = -pkin(3) - t257;
t165 = (-qJ(5) + t218) * t148;
t230 = t149 * t150;
t155 = pkin(4) * t230 + t165;
t231 = t149 * t109;
t292 = -qJD(1) * t155 + qJD(3) * t120 + t231;
t291 = t303 * t149 + t300;
t243 = Ifges(6,6) * t149;
t246 = Ifges(5,6) * t149;
t290 = t302 * t147 + t243 - t246;
t248 = Ifges(4,4) * t149;
t289 = t246 + t248 + (Ifges(4,1) + Ifges(5,2)) * t147;
t227 = qJD(1) * t148;
t216 = t149 * t227;
t108 = qJD(2) * t147 + t216;
t118 = qJD(2) * pkin(7) + t140;
t114 = -pkin(2) * t150 - t148 * pkin(7) - pkin(1);
t98 = t114 * qJD(1);
t55 = t118 * t147 - t149 * t98;
t163 = pkin(4) * t108 + t55;
t287 = qJD(4) + t163;
t217 = t147 * t227;
t221 = t149 * qJD(2);
t107 = t217 - t221;
t286 = -pkin(4) * t107 + qJD(5);
t138 = Ifges(3,4) * t226;
t133 = -qJD(3) + t226;
t256 = pkin(3) + qJ(5);
t10 = t133 * t256 + t287;
t117 = -qJD(2) * pkin(2) + pkin(6) * t227;
t56 = t118 * t149 + t147 * t98;
t168 = t56 * t147 - t55 * t149;
t280 = -qJD(4) - t55;
t33 = pkin(3) * t133 - t280;
t34 = qJ(4) * t133 - t56;
t169 = t34 * t147 + t33 * t149;
t247 = Ifges(5,6) * t147;
t179 = -Ifges(5,2) * t149 + t247;
t184 = -Ifges(4,2) * t147 + t248;
t188 = -mrSges(6,2) * t149 + mrSges(6,3) * t147;
t189 = -mrSges(5,2) * t147 - mrSges(5,3) * t149;
t191 = mrSges(4,1) * t147 + mrSges(4,2) * t149;
t20 = -t34 + t286;
t161 = -qJ(4) * t108 + t117;
t21 = t107 * t256 + t161;
t259 = t149 / 0.2e1;
t260 = -t149 / 0.2e1;
t261 = t147 / 0.2e1;
t262 = -t147 / 0.2e1;
t264 = -t133 / 0.2e1;
t265 = t108 / 0.2e1;
t266 = -t108 / 0.2e1;
t267 = t107 / 0.2e1;
t268 = -t107 / 0.2e1;
t277 = -Ifges(4,6) + t299;
t278 = Ifges(5,4) - t298;
t102 = Ifges(6,6) * t108;
t241 = t108 * Ifges(5,6);
t281 = t302 * t107 - t299 * t133 + t102 - t241;
t104 = Ifges(4,4) * t107;
t245 = Ifges(6,6) * t107;
t282 = t303 * t108 - t298 * t133 - t104 + t245;
t36 = pkin(3) * t107 + t161;
t103 = Ifges(5,6) * t107;
t41 = -Ifges(5,4) * t133 - Ifges(5,2) * t108 + t103;
t242 = t108 * Ifges(4,4);
t42 = -Ifges(4,2) * t107 - Ifges(4,6) * t133 + t242;
t151 = t169 * mrSges(5,1) + (t10 * t149 - t20 * t147) * mrSges(6,1) - t168 * mrSges(4,3) + t117 * t191 + t179 * t266 + t184 * t268 + t21 * t188 + t36 * t189 + t41 * t260 + t42 * t262 + t290 * t267 + t291 * t265 + t281 * t261 + t282 * t259 + (t277 * t147 - t278 * t149) * t264;
t285 = t151 + Ifges(3,1) * t227 / 0.2e1 + t138 / 0.2e1 + t212;
t284 = t247 - t300 + (Ifges(4,2) + t302) * t149;
t211 = -Ifges(3,6) * qJD(2) / 0.2e1;
t283 = pkin(6) * (t148 * t221 + t150 * t223);
t279 = -t56 - t286;
t200 = Ifges(4,6) / 0.2e1 - Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1;
t201 = -Ifges(6,5) / 0.2e1 - Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1;
t202 = Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1;
t250 = Ifges(3,4) * t148;
t276 = t200 * t107 + t201 * t108 + t202 * t133 - t20 * mrSges(6,2) - t33 * mrSges(5,2) - Ifges(4,6) * t268 - Ifges(5,4) * t266 - t211 + (t150 * Ifges(3,2) + t250) * qJD(1) / 0.2e1 + t10 * mrSges(6,3) + t34 * mrSges(5,3) + t55 * mrSges(4,1) + t56 * mrSges(4,2) - t299 * t267 - t298 * t265 - (Ifges(4,3) + Ifges(6,1) + Ifges(5,1)) * t264;
t275 = -t74 / 0.2e1;
t274 = t74 / 0.2e1;
t273 = -t75 / 0.2e1;
t270 = pkin(1) * mrSges(3,1);
t269 = pkin(1) * mrSges(3,2);
t263 = t133 / 0.2e1;
t252 = mrSges(4,3) * t107;
t78 = mrSges(4,2) * t133 - t252;
t81 = mrSges(5,1) * t107 + mrSges(5,3) * t133;
t255 = -t78 + t81;
t251 = mrSges(4,3) * t108;
t79 = -mrSges(4,1) * t133 - t251;
t83 = mrSges(5,1) * t108 - mrSges(5,2) * t133;
t254 = -t79 + t83;
t82 = -mrSges(6,1) * t107 - mrSges(6,2) * t133;
t253 = -t81 + t82;
t112 = t197 * qJD(2);
t240 = t147 * t112 + t114 * t222;
t237 = qJ(4) * t107;
t235 = qJD(2) * mrSges(3,2);
t234 = t147 * t148;
t232 = t148 * t149;
t228 = pkin(3) * t234 + t148 * pkin(6);
t136 = pkin(6) * t230;
t86 = t147 * t114 + t136;
t225 = qJD(2) * t148;
t135 = pkin(6) * t233;
t50 = -t74 * mrSges(5,1) + mrSges(5,2) * t210;
t49 = -t75 * mrSges(6,1) + mrSges(6,2) * t210;
t208 = -qJ(4) * t147 - pkin(2);
t85 = t149 * t114 - t135;
t205 = pkin(6) * t210;
t204 = t288 * pkin(3) + pkin(6) * t224 + qJ(4) * t214;
t203 = m(4) * t117 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t107 + mrSges(4,2) * t108 + mrSges(3,3) * t227;
t99 = qJD(1) * t112;
t199 = t118 * t222 - t149 * t99 + t98 * t223;
t198 = t218 * t148;
t166 = t118 * t223 - t147 * t99 - t98 * t222;
t157 = qJD(4) * t133 + t166;
t162 = t207 * t227;
t5 = qJD(2) * t162 + t157;
t164 = qJD(2) * t198;
t9 = qJD(1) * t164 + t199;
t196 = t147 * t9 - t149 * t5;
t72 = qJ(4) * t150 - t86;
t195 = -qJD(4) * t150 + t240;
t194 = qJD(3) * t136 - t149 * t112 + t114 * t223;
t192 = mrSges(4,1) * t149 - mrSges(4,2) * t147;
t190 = mrSges(5,2) * t149 - mrSges(5,3) * t147;
t187 = mrSges(6,2) * t147 + mrSges(6,3) * t149;
t182 = Ifges(5,4) * t147 + Ifges(5,5) * t149;
t181 = Ifges(6,4) * t149 - Ifges(6,5) * t147;
t180 = Ifges(4,5) * t147 + Ifges(4,6) * t149;
t173 = -Ifges(6,3) * t147 + t243;
t11 = -t149 * t205 - t166;
t12 = t147 * t205 - t199;
t170 = t11 * t149 - t12 * t147;
t47 = -t74 * mrSges(6,1) - mrSges(6,3) * t210;
t160 = pkin(6) * t209 + qJ(4) * t74 - qJD(4) * t108;
t2 = -pkin(4) * t74 + qJD(5) * t133 + t165 * t220 + t199;
t3 = -pkin(4) * t75 - t193 * t220 - t157;
t154 = -t12 * mrSges(4,1) + t11 * mrSges(4,2) - t9 * mrSges(5,2) - t3 * mrSges(6,2) + t5 * mrSges(5,3) + t2 * mrSges(6,3);
t145 = t150 * pkin(3);
t131 = Ifges(5,1) * t210;
t130 = Ifges(6,1) * t210;
t129 = Ifges(4,3) * t210;
t119 = t271 * t147;
t116 = mrSges(3,3) * t226 - t235;
t113 = -pkin(3) * t149 + t208;
t97 = -t149 * t256 + t208;
t89 = -qJ(4) * t232 + t228;
t80 = mrSges(6,1) * t108 + mrSges(6,3) * t133;
t77 = -pkin(6) * t216 + t95;
t76 = pkin(6) * t217 + t231;
t73 = t145 - t85;
t71 = Ifges(5,4) * t74;
t70 = Ifges(6,4) * t75;
t69 = Ifges(4,5) * t74;
t68 = Ifges(5,5) * t75;
t67 = Ifges(6,5) * t74;
t66 = Ifges(4,6) * t75;
t64 = t148 * t167 + t228;
t62 = -mrSges(5,2) * t107 - mrSges(5,3) * t108;
t60 = pkin(3) * t108 + t237;
t59 = -mrSges(6,2) * t108 + mrSges(6,3) * t107;
t58 = qJD(1) * t198 - t231;
t57 = -t95 + t162;
t53 = -pkin(4) * t234 - t72;
t52 = -mrSges(4,2) * t210 - mrSges(4,3) * t75;
t51 = mrSges(4,1) * t210 + mrSges(4,3) * t74;
t48 = mrSges(5,1) * t75 - mrSges(5,3) * t210;
t46 = qJ(5) * t150 + t135 + t145 + (pkin(4) * t148 - t114) * t149;
t31 = t108 * t256 + t237;
t28 = t225 * t257 - t194;
t27 = t240 - t283;
t26 = (-qJ(4) * t224 - qJD(4) * t148) * t149 + t204;
t25 = t164 + t194;
t24 = mrSges(6,2) * t74 + mrSges(6,3) * t75;
t23 = mrSges(4,1) * t75 - mrSges(4,2) * t74;
t22 = -mrSges(5,2) * t75 + mrSges(5,3) * t74;
t19 = -qJ(4) * t225 - t195 + t283;
t17 = -t74 * Ifges(4,4) - t75 * Ifges(4,2) + Ifges(4,6) * t210;
t16 = Ifges(5,4) * t210 + t74 * Ifges(5,2) + t75 * Ifges(5,6);
t8 = qJD(2) * t159 + (qJD(5) * t147 + (qJ(5) * qJD(3) - qJD(4)) * t149) * t148 + t204;
t7 = (-pkin(4) * t232 - t135) * qJD(3) + t156 * qJD(2) + t195;
t6 = pkin(3) * t75 + t160;
t4 = -pkin(4) * t214 + qJD(2) * t155 + qJD(5) * t150 + t194;
t1 = qJD(5) * t107 + t256 * t75 + t160;
t13 = [t85 * t51 + t86 * t52 + t89 * t22 + t27 * t78 + t28 * t79 + t4 * t80 + t19 * t81 + t7 * t82 + t25 * t83 + t72 * t48 + t73 * t50 + t8 * t59 + t26 * t62 + t64 * t24 + t46 * t47 + t53 * t49 + m(5) * (t19 * t34 + t25 * t33 + t26 * t36 + t5 * t72 + t6 * t89 + t73 * t9) + m(4) * (t11 * t86 + t12 * t85 + t56 * t27 - t55 * t28) + m(6) * (t1 * t64 + t10 * t4 + t2 * t46 + t20 * t7 + t21 * t8 + t3 * t53) + (-t130 / 0.2e1 - t131 / 0.2e1 - t129 / 0.2e1 - t70 / 0.2e1 - t71 / 0.2e1 + t67 / 0.2e1 - t68 / 0.2e1 + t69 / 0.2e1 + t66 / 0.2e1 + (t203 * pkin(6) + (-0.2e1 * t269 + 0.3e1 / 0.2e1 * Ifges(3,4) * t150) * qJD(1) + t212 + t285) * qJD(2) + t200 * t75 - t201 * t74 + t154) * t150 + (t6 * t189 + t16 * t260 + t17 * t262 + pkin(6) * t23 + t184 * t273 + t179 * t274 + t1 * t188 + (-t11 * t147 - t12 * t149) * mrSges(4,3) + (-t147 * t3 + t149 * t2) * mrSges(6,1) + (t147 * t5 + t149 * t9) * mrSges(5,1) + (-pkin(6) * t116 + t211 - t276) * qJD(2) + (t180 * t263 + t117 * t192 + t173 * t265 - t36 * t190 + t21 * t187 + t42 * t260 + (-t147 * t55 - t149 * t56) * mrSges(4,3) + (-t10 * t147 - t149 * t20) * mrSges(6,1) + (-t147 * t33 + t149 * t34) * mrSges(5,1) + t289 * t266 + (t181 + t182) * t264 + t282 * t262 + t284 * t267) * qJD(3) + (-0.3e1 / 0.2e1 * t250 - 0.2e1 * t270 - t201 * t232 - t200 * t234 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + (m(4) * pkin(6) + t191) * pkin(6) - t202) * t150) * t220 + t291 * t275 + t290 * t75 / 0.2e1 + (t41 * qJD(3) + t297) * t261 + (t281 * qJD(3) + t296) * t259) * t148; -t1 * t187 + t6 * t190 + t170 * mrSges(4,3) + t151 * qJD(3) + ((-t48 + t52) * t149 + (t50 - t51) * t147 + m(4) * t170 + m(5) * t196 + (-m(4) * t168 + m(5) * t169 + t147 * t255 + t149 * t254) * qJD(3)) * pkin(7) + t292 * t80 + t293 * t82 + t289 * t275 + t284 * t273 + ((t211 + (t116 + t235) * pkin(6) + (t250 / 0.2e1 + t270) * qJD(1) + (-t182 / 0.2e1 + t180 / 0.2e1 - t181 / 0.2e1) * qJD(2) + t276) * t148 + (((-m(4) * pkin(2) - mrSges(3,1) - t192) * qJD(2) - t203) * pkin(6) - t138 / 0.2e1 + (t269 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t148) * qJD(1) + t212 - t285) * t150) * qJD(1) + t113 * t22 + t119 * t47 + t120 * t49 + (t1 * t97 + t292 * t10 + t119 * t2 + t120 * t3 + t293 * t20 + t295 * t21) * m(6) + t295 * t59 + t296 * t261 + t297 * t260 + (t113 * t6 + t294 * t36 - t33 * t58 - t34 * t57) * m(5) + t294 * t62 + t97 * t24 - t77 * t78 - t76 * t79 - t57 * t81 - t58 * t83 - pkin(2) * t23 + (t147 * t2 + t149 * t3) * mrSges(6,1) - m(4) * (-t55 * t76 + t56 * t77) + t196 * mrSges(5,1) + t17 * t259 + t16 * t262 + t173 * t274; (-t303 * t107 + t102 - t242 + t281) * t266 + (t10 * t107 + t108 * t20) * mrSges(6,1) + (t107 * t33 - t108 * t34) * mrSges(5,1) + (t302 * t108 + t103 - t245 + t41) * t268 + (qJ(4) * t3 + t279 * t10 - t256 * t2 + t287 * t20 - t21 * t31) * m(6) - t117 * (mrSges(4,1) * t108 - mrSges(4,2) * t107) - t36 * (-mrSges(5,2) * t108 + mrSges(5,3) * t107) - t21 * (mrSges(6,2) * t107 + mrSges(6,3) * t108) + t130 + t131 + t129 + t70 + t71 - t67 + t68 - t69 - t66 + (-Ifges(4,2) * t108 - t104 + t282) * t267 - t31 * t59 - t60 * t62 - pkin(3) * t50 + t279 * t80 + (-pkin(3) * t9 - qJ(4) * t5 + t280 * t34 - t33 * t56 - t36 * t60) * m(5) + (t278 * t107 + t277 * t108) * t263 + (Ifges(5,2) * t107 + t241 + t42) * t265 + t163 * t82 - t256 * t47 - t154 + (-t48 + t49) * qJ(4) + t253 * qJD(4) + (t251 - t254) * t56 + (t252 - t255) * t55; t253 * t133 + (t59 + t62) * t108 + t47 + t50 + (t108 * t21 + t133 * t20 + t2) * m(6) + (t108 * t36 - t133 * t34 + t9) * m(5); -t107 * t59 - t133 * t80 + 0.2e1 * (t3 / 0.2e1 + t10 * t264 + t21 * t268) * m(6) + t49;];
tauc = t13(:);
