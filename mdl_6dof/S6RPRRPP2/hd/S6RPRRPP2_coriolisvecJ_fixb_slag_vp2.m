% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:35
% EndTime: 2019-03-09 04:31:48
% DurationCPUTime: 7.20s
% Computational Cost: add. (3404->507), mult. (8191->607), div. (0->0), fcn. (4522->6), ass. (0->227)
t315 = Ifges(7,4) + Ifges(6,5);
t294 = Ifges(6,4) + Ifges(5,5);
t316 = Ifges(7,5) - t294;
t311 = Ifges(7,2) + Ifges(6,3);
t314 = Ifges(5,6) - Ifges(6,6);
t313 = Ifges(6,6) - Ifges(7,6);
t307 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t312 = qJD(3) / 0.2e1;
t154 = cos(qJ(4));
t310 = t315 * t154;
t152 = sin(qJ(4));
t309 = t315 * t152;
t143 = sin(pkin(9)) * pkin(1) + pkin(7);
t131 = t143 * qJD(1);
t153 = sin(qJ(3));
t155 = cos(qJ(3));
t105 = t155 * qJD(2) - t153 * t131;
t226 = qJD(1) * t155;
t146 = Ifges(4,4) * t226;
t205 = Ifges(4,5) * t312;
t142 = -qJD(4) + t226;
t277 = pkin(4) + pkin(5);
t224 = qJD(3) * t152;
t227 = qJD(1) * t153;
t125 = t154 * t227 + t224;
t206 = -cos(pkin(9)) * pkin(1) - pkin(2);
t117 = -pkin(3) * t155 - pkin(8) * t153 + t206;
t102 = t117 * qJD(1);
t225 = qJD(2) * t153;
t106 = t131 * t155 + t225;
t99 = qJD(3) * pkin(8) + t106;
t36 = t154 * t102 - t152 * t99;
t17 = qJ(6) * t125 + t36;
t291 = qJD(5) - t17;
t11 = t142 * t277 + t291;
t138 = t142 * qJ(5);
t215 = t154 * qJD(3);
t124 = t152 * t227 - t215;
t37 = t152 * t102 + t154 * t99;
t18 = qJ(6) * t124 + t37;
t12 = -t138 + t18;
t197 = qJD(3) * pkin(3) + t105;
t162 = qJ(5) * t125 + t197;
t16 = -t124 * t277 + qJD(6) + t162;
t166 = t152 * t37 + t154 * t36;
t290 = qJD(5) - t36;
t27 = pkin(4) * t142 + t290;
t28 = -t138 + t37;
t167 = t152 * t28 - t154 * t27;
t248 = Ifges(5,4) * t154;
t179 = -Ifges(5,2) * t152 + t248;
t186 = -mrSges(7,1) * t152 + mrSges(7,2) * t154;
t188 = mrSges(6,1) * t152 - mrSges(6,3) * t154;
t190 = mrSges(5,1) * t152 + mrSges(5,2) * t154;
t267 = t154 / 0.2e1;
t269 = t152 / 0.2e1;
t270 = -t152 / 0.2e1;
t271 = t142 / 0.2e1;
t272 = -t142 / 0.2e1;
t273 = t125 / 0.2e1;
t275 = t124 / 0.2e1;
t276 = -t124 / 0.2e1;
t249 = Ifges(5,4) * t152;
t284 = t154 * t307 - t249 + t309;
t121 = Ifges(5,4) * t124;
t297 = t315 * t124;
t285 = t307 * t125 + t142 * t316 - t121 + t297;
t287 = t152 * t311 + t310;
t300 = t315 * t125;
t292 = t311 * t124 - t313 * t142 + t300;
t35 = pkin(4) * t124 - t162;
t250 = Ifges(5,4) * t125;
t49 = -t124 * Ifges(5,2) - t142 * Ifges(5,6) + t250;
t282 = t167 * mrSges(6,2) + t166 * mrSges(5,3) + (t11 * t154 - t12 * t152) * mrSges(7,3) - t49 * t270 - (Ifges(7,5) * t154 + Ifges(7,6) * t152) * t271 - t179 * t276 + t197 * t190 - t35 * t188 - t16 * t186 - t287 * t275 - (-t152 * t314 + t154 * t294) * t272 - t292 * t269 - t284 * t273 - t285 * t267;
t306 = t227 / 0.2e1;
t308 = Ifges(4,1) * t306 + t146 / 0.2e1 + t205 - t105 * mrSges(4,3) - t282;
t305 = -qJD(3) / 0.2e1;
t223 = qJD(3) * t153;
t202 = qJD(1) * t223;
t214 = qJD(3) * qJD(4);
t221 = qJD(4) * t152;
t84 = t154 * t214 + (-t153 * t221 + t155 * t215) * qJD(1);
t220 = qJD(4) * t154;
t222 = qJD(3) * t155;
t85 = t152 * t214 + (t152 * t222 + t153 * t220) * qJD(1);
t304 = t313 * t202 + t311 * t85 + t315 * t84;
t234 = qJ(5) * t154;
t163 = -t152 * t277 + t234;
t219 = qJD(5) * t152;
t303 = qJD(4) * t163 + t219 + t225 - (qJD(1) * t163 - t131) * t155;
t216 = qJD(6) * t154;
t260 = pkin(8) - qJ(6);
t196 = pkin(3) * t153 - pkin(8) * t155;
t128 = t196 * qJD(1);
t54 = t154 * t105 + t152 * t128;
t40 = qJ(5) * t227 + t54;
t302 = -qJ(6) * t152 * t226 - t221 * t260 - t216 - t40;
t136 = t260 * t154;
t230 = t154 * t155;
t53 = -t152 * t105 + t128 * t154;
t301 = qJD(4) * t136 - qJD(6) * t152 - (-qJ(6) * t230 - t153 * t277) * qJD(1) + t53;
t299 = -t311 * t154 + t309;
t298 = t294 * t152 + t314 * t154;
t100 = qJD(3) * t105;
t296 = (-Ifges(5,4) + t315) * t85 + t307 * t84 - t316 * t202;
t295 = t307 * t152 + t248 - t310;
t274 = -t125 / 0.2e1;
t204 = Ifges(4,6) * t305;
t101 = qJD(3) * t106;
t289 = -qJ(5) * t84 - qJD(5) * t125 + t101;
t288 = qJ(5) * t223 - qJD(5) * t155;
t129 = t196 * qJD(3);
t116 = qJD(1) * t129;
t8 = t154 * t100 + t102 * t220 + t152 * t116 - t221 * t99;
t9 = -t152 * t100 - t102 * t221 + t116 * t154 - t99 * t220;
t193 = -t152 * t9 + t154 * t8;
t4 = qJ(5) * t202 - t142 * qJD(5) + t8;
t5 = -pkin(4) * t202 - t9;
t194 = t152 * t5 + t154 * t4;
t235 = qJ(5) * t152;
t283 = -t154 * t277 - t235;
t133 = t206 * qJD(1);
t209 = Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t198 = Ifges(7,6) / 0.2e1 + t209;
t211 = -Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1;
t200 = Ifges(7,5) / 0.2e1 + t211;
t210 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t251 = Ifges(4,4) * t153;
t281 = -t198 * t124 - t200 * t125 - (Ifges(7,3) / 0.2e1 + t210) * t142 + t12 * mrSges(7,2) + t133 * mrSges(4,1) + t28 * mrSges(6,3) + t36 * mrSges(5,1) + t204 - (Ifges(4,2) * t155 + t251) * qJD(1) / 0.2e1 + Ifges(7,5) * t274 + Ifges(6,6) * t275 - t11 * mrSges(7,1) - t27 * mrSges(6,1) - t37 * mrSges(5,2) + (Ifges(7,6) + Ifges(5,6)) * t276 + t294 * t273 + (Ifges(7,3) + Ifges(5,3) + Ifges(6,2)) * t272;
t280 = t84 / 0.2e1;
t279 = -t85 / 0.2e1;
t278 = t85 / 0.2e1;
t268 = -t154 / 0.2e1;
t262 = t84 * mrSges(7,3);
t58 = -mrSges(6,2) * t85 + mrSges(6,3) * t202;
t63 = -mrSges(5,2) * t202 - mrSges(5,3) * t85;
t259 = t58 + t63;
t60 = mrSges(5,1) * t202 - mrSges(5,3) * t84;
t73 = t84 * mrSges(6,2);
t61 = -mrSges(6,1) * t202 + t73;
t258 = -t60 + t61;
t68 = mrSges(6,1) * t124 - mrSges(6,3) * t125;
t69 = -mrSges(7,1) * t124 + mrSges(7,2) * t125;
t257 = t68 - t69;
t89 = -mrSges(7,2) * t142 + mrSges(7,3) * t124;
t94 = -mrSges(6,2) * t124 - mrSges(6,3) * t142;
t256 = t89 + t94;
t253 = mrSges(5,3) * t124;
t90 = mrSges(5,2) * t142 - t253;
t255 = t90 + t94;
t252 = mrSges(5,3) * t125;
t92 = -mrSges(5,1) * t142 - t252;
t93 = mrSges(6,1) * t142 + mrSges(6,2) * t125;
t254 = -t92 + t93;
t240 = t133 * mrSges(4,2);
t208 = mrSges(4,3) * t227;
t239 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t124 + mrSges(5,2) * t125 + t208;
t236 = qJ(5) * t124;
t233 = t143 * t152;
t232 = t152 * t153;
t231 = t153 * t154;
t127 = t143 * t230;
t229 = qJD(4) * t127 + t117 * t221;
t228 = t117 * t220 + t152 * t129;
t66 = t152 * t117 + t127;
t218 = qJD(5) * t154;
t213 = t89 + t255;
t91 = mrSges(7,1) * t142 - mrSges(7,3) * t125;
t212 = t91 + t254;
t207 = mrSges(4,3) * t226;
t203 = m(4) * t143 + mrSges(4,3);
t31 = -t85 * mrSges(7,1) + t84 * mrSges(7,2);
t201 = -pkin(4) - t233;
t126 = t155 * t233;
t65 = t117 * t154 - t126;
t1 = qJ(6) * t85 + qJD(6) * t124 + t4;
t2 = -qJ(6) * t84 - qJD(6) * t125 - t202 * t277 - t9;
t195 = -t1 * t154 - t152 * t2;
t56 = -qJ(5) * t155 + t66;
t192 = -t129 * t154 + t229;
t191 = mrSges(5,1) * t154 - mrSges(5,2) * t152;
t189 = mrSges(6,1) * t154 + mrSges(6,3) * t152;
t187 = mrSges(7,1) * t154 + mrSges(7,2) * t152;
t178 = Ifges(5,2) * t154 + t249;
t171 = Ifges(7,5) * t152 - Ifges(7,6) * t154;
t170 = pkin(4) * t154 + t235;
t169 = pkin(4) * t152 - t234;
t164 = t143 + t169;
t161 = -t143 + t163;
t19 = (-t153 * t215 - t155 * t221) * t143 + t228;
t160 = -t9 * mrSges(5,1) + t5 * mrSges(6,1) + t2 * mrSges(7,1) + t8 * mrSges(5,2) - t1 * mrSges(7,2) - t4 * mrSges(6,3);
t149 = t155 * pkin(4);
t141 = Ifges(6,2) * t202;
t140 = Ifges(5,3) * t202;
t135 = t260 * t152;
t134 = -qJD(3) * mrSges(4,2) + t207;
t130 = -pkin(3) - t170;
t118 = pkin(3) - t283;
t110 = qJD(4) * t169 - t219;
t86 = t164 * t153;
t77 = Ifges(6,4) * t84;
t76 = Ifges(5,5) * t84;
t75 = Ifges(5,6) * t85;
t74 = Ifges(6,6) * t85;
t67 = pkin(4) * t125 + t236;
t64 = t161 * t153;
t62 = mrSges(7,2) * t202 + mrSges(7,3) * t85;
t59 = -mrSges(7,1) * t202 - t262;
t57 = t149 - t65;
t55 = t225 + (qJD(1) * t169 + t131) * t155;
t43 = qJ(6) * t232 + t56;
t42 = -pkin(4) * t227 - t53;
t39 = -t125 * t277 - t236;
t38 = pkin(5) * t155 + t126 + t149 + (-qJ(6) * t153 - t117) * t154;
t33 = (qJD(4) * t170 - t218) * t153 + t164 * t222;
t32 = mrSges(5,1) * t85 + mrSges(5,2) * t84;
t30 = mrSges(6,1) * t85 - mrSges(6,3) * t84;
t23 = t84 * Ifges(5,4) - t85 * Ifges(5,2) + Ifges(5,6) * t202;
t20 = t223 * t233 - t192;
t15 = t201 * t223 + t192;
t14 = (qJD(4) * t283 + t218) * t153 + t161 * t222;
t13 = t19 + t288;
t10 = pkin(4) * t85 + t289;
t7 = (qJ(6) * qJD(4) - qJD(3) * t143) * t231 + (qJD(6) * t153 + (qJ(6) * qJD(3) - qJD(4) * t143) * t155) * t152 + t228 + t288;
t6 = (-qJ(6) * t222 - t129) * t154 + (qJ(6) * t221 - t216 + (-pkin(5) + t201) * qJD(3)) * t153 + t229;
t3 = -t277 * t85 - t289;
t21 = [t86 * t30 + t7 * t89 + t19 * t90 + t6 * t91 + t20 * t92 + t15 * t93 + t13 * t94 + t33 * t68 + t14 * t69 + t56 * t58 + t38 * t59 + t57 * t61 + t43 * t62 + t64 * t31 + t65 * t60 + t66 * t63 + m(5) * (t19 * t37 + t20 * t36 + t65 * t9 + t66 * t8) + m(6) * (t10 * t86 + t13 * t28 + t15 * t27 + t33 * t35 + t4 * t56 + t5 * t57) + m(7) * (t1 * t43 + t11 * t6 + t12 * t7 + t14 * t16 + t2 * t38 + t3 * t64) + (-t140 / 0.2e1 - t141 / 0.2e1 - t77 / 0.2e1 - t74 / 0.2e1 + t75 / 0.2e1 - t76 / 0.2e1 + t203 * t100 + (Ifges(7,6) + t209) * t85 + (Ifges(7,5) + t211) * t84 + (0.3e1 / 0.2e1 * t146 + 0.2e1 * t240 + (-m(4) * t105 - m(5) * t197 + t239) * t143 + t205 + t308) * qJD(3) + t160) * t155 + (t10 * t188 + t3 * t186 + t179 * t279 + (t1 * t152 - t154 * t2) * mrSges(7,3) + (-t152 * t8 - t154 * t9) * mrSges(5,3) + (-t152 * t4 + t154 * t5) * mrSges(6,2) + (mrSges(4,3) + t190) * t101 + (-t203 * t106 + (t206 * mrSges(4,1) - 0.3e1 / 0.2e1 * t251 - t200 * t231 - t198 * t232 + (-Ifges(7,3) - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - t210) * t155) * qJD(1) + t204 + t281) * qJD(3) + t287 * t278 + t304 * t269 + t284 * t280 + t23 * t270 + t296 * t267 + (t32 + (m(4) + m(5)) * t101 - qJD(3) * t134) * t143 + (t49 * t268 + t171 * t272 + t178 * t275 - t197 * t191 + t35 * t189 - t16 * t187 + (t11 * t152 + t12 * t154) * mrSges(7,3) + (t152 * t36 - t154 * t37) * mrSges(5,3) + (-t152 * t27 - t154 * t28) * mrSges(6,2) + t299 * t276 + t298 * t271 + t295 * t274 + t285 * t270 + t292 * t267) * qJD(4)) * t153; (-t30 + t31 - t32 + (t152 * t212 + t154 * t213 + t134 - t207) * qJD(3) + m(5) * (t215 * t37 - t224 * t36 - t101) + m(6) * (t215 * t28 + t224 * t27 - t10) + m(7) * (t11 * t224 + t12 * t215 + t3)) * t155 + ((t62 + t259) * t154 + (t59 + t258) * t152 + (-t208 + t239 + t257) * qJD(3) + (-t152 * t213 + t154 * t212) * qJD(4) + m(5) * (-qJD(3) * t197 - t220 * t36 - t221 * t37 + t193) + m(6) * (qJD(3) * t35 + t220 * t27 - t221 * t28 + t194) + m(7) * (-qJD(3) * t16 + t11 * t220 - t12 * t221 - t195)) * t153; (-t282 + (-m(5) * t166 - m(6) * t167 - t152 * t255 + t154 * t254) * pkin(8)) * qJD(4) + t130 * t30 - t105 * t134 + t135 * t59 + t136 * t62 + t118 * t31 - t54 * t90 - t53 * t92 - t42 * t93 - t40 * t94 - t100 * mrSges(4,2) - pkin(3) * t32 + t299 * t278 + t304 * t268 - t239 * t106 + t301 * t91 + t302 * t89 + (t152 * t258 + t154 * t259) * pkin(8) + ((t106 * mrSges(4,3) + Ifges(4,4) * t306 + t171 * t305 + t298 * t312 + t204 - t281) * t153 + ((Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t227 - t146 / 0.2e1 - t240 + t205 - t308) * t155) * qJD(1) + t193 * mrSges(5,3) + t194 * mrSges(6,2) + t195 * mrSges(7,3) + (-mrSges(4,1) - t191) * t101 + t3 * t187 - t10 * t189 - m(6) * (t27 * t42 + t28 * t40 + t35 * t55) - m(5) * (-t106 * t197 + t36 * t53 + t37 * t54) + t295 * t280 + t296 * t269 + (t1 * t136 + t11 * t301 + t118 * t3 + t12 * t302 + t135 * t2 + t16 * t303) * m(7) + t303 * t69 + t178 * t279 + t23 * t267 + m(6) * (pkin(8) * t194 + t10 * t130 + t110 * t35) + m(5) * (-pkin(3) * t101 + pkin(8) * t193) + (t110 - t55) * t68; t140 + t141 + t77 + t74 - t75 + t76 - t35 * (mrSges(6,1) * t125 + mrSges(6,3) * t124) - t16 * (-mrSges(7,1) * t125 - mrSges(7,2) * t124) - t17 * t89 - t18 * t91 - Ifges(7,5) * t84 - Ifges(7,6) * t85 - t67 * t68 - t39 * t69 - pkin(4) * t61 + (t252 - t254) * t37 + (-t253 - t255) * t36 + t256 * qJD(5) + (-t124 * t307 - t250 + t292 + t300) * t274 - t160 + (t58 + t62) * qJ(5) + t197 * (mrSges(5,1) * t125 - mrSges(5,2) * t124) - t277 * t59 + (t1 * qJ(5) - t11 * t18 + t12 * t291 - t16 * t39 - t2 * t277) * m(7) + (-t11 * t124 - t12 * t125) * mrSges(7,3) + (t124 * t27 + t125 * t28) * mrSges(6,2) + (-Ifges(7,5) * t124 + Ifges(7,6) * t125) * t272 + t49 * t273 + (-pkin(4) * t5 + qJ(5) * t4 - t27 * t37 + t28 * t290 - t35 * t67) * m(6) + (-Ifges(5,2) * t125 - t121 + t285) * t275 + Ifges(7,3) * t202 + (t125 * t311 - t297) * t276 + (-t124 * t294 - t125 * t314) * t271; -t262 + t73 + t256 * t142 + t257 * t125 + (-mrSges(6,1) - mrSges(7,1)) * t202 + (t12 * t142 - t125 * t16 + t2) * m(7) + (t125 * t35 + t142 * t28 + t5) * m(6); -t124 * t89 + t125 * t91 + 0.2e1 * (t3 / 0.2e1 + t11 * t273 + t12 * t276) * m(7) + t31;];
tauc  = t21(:);
