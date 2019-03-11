% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:01:21
% EndTime: 2019-03-09 03:01:34
% DurationCPUTime: 7.88s
% Computational Cost: add. (5036->469), mult. (12427->614), div. (0->0), fcn. (8271->8), ass. (0->223)
t321 = Ifges(6,4) + Ifges(7,4);
t322 = Ifges(6,1) + Ifges(7,1);
t315 = Ifges(7,5) + Ifges(6,5);
t320 = Ifges(6,2) + Ifges(7,2);
t319 = Ifges(6,6) + Ifges(7,6);
t155 = cos(qJ(5));
t324 = t321 * t155;
t153 = sin(qJ(5));
t323 = t321 * t153;
t154 = sin(qJ(3));
t156 = cos(qJ(3));
t246 = sin(pkin(10));
t247 = cos(pkin(10));
t135 = t247 * t154 + t246 * t156;
t125 = t135 * qJD(1);
t101 = qJD(3) * t155 - t125 * t153;
t318 = t321 * t101;
t297 = -t153 * t319 + t155 * t315;
t296 = -t153 * t320 + t324;
t295 = t155 * t322 - t323;
t102 = qJD(3) * t153 + t125 * t155;
t317 = t321 * t102;
t159 = -t154 * t246 + t156 * t247;
t124 = t159 * qJD(1);
t279 = -t124 / 0.2e1;
t316 = Ifges(5,2) * t279 - Ifges(5,6) * qJD(3) / 0.2e1;
t309 = qJD(5) - t124;
t304 = t101 * t320 + t309 * t319 + t317;
t303 = t102 * t322 + t315 * t309 + t318;
t313 = t153 * t315 + t155 * t319;
t312 = t155 * t320 + t323;
t311 = t153 * t322 + t324;
t227 = qJD(5) * t155;
t228 = qJD(5) * t153;
t146 = sin(pkin(9)) * pkin(1) + pkin(7);
t138 = t146 * qJD(1);
t192 = qJ(4) * qJD(1) + t138;
t226 = t154 * qJD(2);
t104 = t156 * t192 + t226;
t225 = t154 * qJD(4);
t157 = -qJD(1) * t225 - qJD(3) * t104;
t150 = t156 * qJD(2);
t230 = qJD(3) * t154;
t111 = qJD(3) * t150 - t138 * t230;
t229 = qJD(4) * t156;
t92 = (-qJ(4) * t230 + t229) * qJD(1) + t111;
t33 = t157 * t246 + t247 * t92;
t198 = t247 * t104;
t103 = -t154 * t192 + t150;
t97 = qJD(3) * pkin(3) + t103;
t49 = t246 * t97 + t198;
t46 = qJD(3) * pkin(8) + t49;
t210 = -cos(pkin(9)) * pkin(1) - pkin(2);
t137 = -pkin(3) * t156 + t210;
t233 = qJD(1) * t137;
t123 = qJD(4) + t233;
t60 = -t124 * pkin(4) - t125 * pkin(8) + t123;
t126 = t135 * qJD(3);
t116 = qJD(1) * t126;
t127 = t159 * qJD(3);
t117 = qJD(1) * t127;
t224 = qJD(1) * qJD(3);
t203 = t154 * t224;
t191 = pkin(3) * t203;
t70 = pkin(4) * t116 - pkin(8) * t117 + t191;
t4 = t153 * t70 + t155 * t33 + t60 * t227 - t228 * t46;
t20 = t153 * t60 + t155 * t46;
t5 = -qJD(5) * t20 - t153 * t33 + t155 * t70;
t189 = -t153 * t5 + t155 * t4;
t19 = -t153 * t46 + t155 * t60;
t310 = -t19 * t227 - t20 * t228 + t189;
t244 = Ifges(5,5) * qJD(3);
t308 = -t244 / 0.2e1 - t123 * mrSges(5,2);
t68 = qJD(5) * t101 + t117 * t155;
t69 = -qJD(5) * t102 - t117 * t153;
t306 = t116 * t319 + t320 * t69 + t321 * t68;
t305 = t315 * t116 + t321 * t69 + t322 * t68;
t208 = t246 * pkin(3);
t145 = t208 + pkin(8);
t234 = qJ(6) + t145;
t193 = qJD(5) * t234;
t94 = t246 * t104;
t51 = t103 * t247 - t94;
t232 = qJD(1) * t154;
t221 = pkin(3) * t232;
t82 = pkin(4) * t125 - pkin(8) * t124 + t221;
t22 = -t153 * t51 + t155 * t82;
t239 = t124 * t155;
t301 = -pkin(5) * t125 + qJ(6) * t239 - qJD(6) * t153 - t155 * t193 - t22;
t23 = t153 * t82 + t155 * t51;
t240 = t124 * t153;
t300 = qJ(6) * t240 + qJD(6) * t155 - t153 * t193 - t23;
t50 = t103 * t246 + t198;
t299 = -t50 + (t228 - t240) * pkin(5);
t252 = t125 * mrSges(5,3);
t248 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t101 + mrSges(6,2) * t102 + t252;
t183 = mrSges(7,1) * t153 + mrSges(7,2) * t155;
t185 = mrSges(6,1) * t153 + mrSges(6,2) * t155;
t48 = t247 * t97 - t94;
t45 = -qJD(3) * pkin(4) - t48;
t28 = -t101 * pkin(5) + qJD(6) + t45;
t298 = -t28 * t183 - t45 * t185;
t1 = pkin(5) * t116 - qJ(6) * t68 - qJD(6) * t102 + t5;
t11 = qJ(6) * t101 + t20;
t2 = qJ(6) * t69 + qJD(6) * t101 + t4;
t10 = -qJ(6) * t102 + t19;
t9 = pkin(5) * t309 + t10;
t292 = -t1 * t153 - t11 * t228 + t2 * t155 - t9 * t227;
t250 = t125 * Ifges(5,4);
t287 = -t250 / 0.2e1 + t316;
t283 = t102 / 0.2e1;
t291 = -t298 + t296 * t101 / 0.2e1 + t295 * t283 + t297 * t309 / 0.2e1;
t290 = t123 * mrSges(5,1) + t19 * mrSges(6,1) + t9 * mrSges(7,1) - t20 * mrSges(6,2) - t11 * mrSges(7,2) + t287;
t289 = t68 / 0.2e1;
t288 = t69 / 0.2e1;
t286 = -t101 / 0.2e1;
t284 = -t102 / 0.2e1;
t282 = t116 / 0.2e1;
t281 = -t309 / 0.2e1;
t278 = -t125 / 0.2e1;
t277 = -t153 / 0.2e1;
t274 = t155 / 0.2e1;
t32 = -t247 * t157 + t246 * t92;
t235 = qJ(4) + t146;
t131 = t235 * t154;
t133 = t235 * t156;
t89 = t247 * t131 + t133 * t246;
t269 = t32 * t89;
t41 = mrSges(7,1) * t116 - mrSges(7,3) * t68;
t42 = mrSges(6,1) * t116 - mrSges(6,3) * t68;
t268 = -t41 - t42;
t43 = -mrSges(7,2) * t116 + mrSges(7,3) * t69;
t44 = -mrSges(6,2) * t116 + mrSges(6,3) * t69;
t267 = t43 + t44;
t71 = -mrSges(7,2) * t309 + mrSges(7,3) * t101;
t72 = -mrSges(6,2) * t309 + mrSges(6,3) * t101;
t266 = t71 + t72;
t73 = mrSges(7,1) * t309 - mrSges(7,3) * t102;
t74 = mrSges(6,1) * t309 - mrSges(6,3) * t102;
t265 = t73 + t74;
t90 = -t131 * t246 + t133 * t247;
t84 = t155 * t90;
t85 = -pkin(4) * t159 - pkin(8) * t135 + t137;
t30 = t153 * t85 + t84;
t264 = Ifges(4,4) * t154;
t121 = Ifges(5,4) * t124;
t256 = t116 * mrSges(5,3);
t255 = t117 * mrSges(5,3);
t254 = t124 * mrSges(5,3);
t251 = t125 * Ifges(5,1);
t249 = t159 * t32;
t245 = Ifges(4,5) * qJD(3);
t243 = Ifges(4,6) * qJD(3);
t241 = qJ(6) * t135;
t238 = t127 * t155;
t237 = t135 * t153;
t236 = t153 * t127;
t231 = qJD(1) * t156;
t194 = qJD(3) * t235;
t107 = -t154 * t194 + t229;
t108 = -t156 * t194 - t225;
t59 = t107 * t247 + t108 * t246;
t83 = pkin(3) * t230 + pkin(4) * t126 - pkin(8) * t127;
t223 = t153 * t83 + t155 * t59 + t227 * t85;
t56 = -mrSges(7,1) * t101 + mrSges(7,2) * t102;
t222 = t56 + t248;
t218 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t217 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t216 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t215 = mrSges(4,3) * t232;
t214 = mrSges(4,3) * t231;
t149 = Ifges(4,4) * t231;
t209 = t247 * pkin(3);
t207 = t135 * t227;
t204 = m(4) * t146 + mrSges(4,3);
t24 = -mrSges(7,1) * t69 + mrSges(7,2) * t68;
t199 = -t153 * t59 + t155 * t83;
t29 = -t153 * t90 + t155 * t85;
t195 = mrSges(5,1) * t116 + mrSges(5,2) * t117;
t147 = -t209 - pkin(4);
t190 = -t1 * t155 - t153 * t2;
t188 = -t153 * t4 - t155 * t5;
t140 = t210 * qJD(1);
t187 = t11 * t155 - t153 * t9;
t186 = mrSges(6,1) * t155 - mrSges(6,2) * t153;
t184 = mrSges(7,1) * t155 - mrSges(7,2) * t153;
t170 = -t153 * t20 - t155 * t19;
t169 = t153 * t19 - t155 * t20;
t58 = t107 * t246 - t108 * t247;
t168 = -qJ(6) * t127 - qJD(6) * t135;
t119 = t138 * t156 + t226;
t158 = t5 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2);
t141 = -qJD(3) * mrSges(4,2) + t214;
t139 = qJD(3) * mrSges(4,1) - t215;
t136 = -pkin(5) * t155 + t147;
t132 = t234 * t155;
t130 = t234 * t153;
t129 = Ifges(4,1) * t232 + t149 + t245;
t128 = t243 + (Ifges(4,2) * t156 + t264) * qJD(1);
t118 = -t138 * t154 + t150;
t115 = Ifges(6,3) * t116;
t114 = Ifges(7,3) * t116;
t112 = t119 * qJD(3);
t109 = -qJD(3) * mrSges(5,2) + t254;
t91 = -mrSges(5,1) * t124 + mrSges(5,2) * t125;
t87 = t121 + t244 + t251;
t67 = Ifges(6,5) * t68;
t66 = Ifges(7,5) * t68;
t65 = Ifges(6,6) * t69;
t64 = Ifges(7,6) * t69;
t55 = pkin(5) * t237 + t89;
t36 = t102 * Ifges(6,5) + t101 * Ifges(6,6) + Ifges(6,3) * t309;
t35 = t102 * Ifges(7,5) + t101 * Ifges(7,6) + Ifges(7,3) * t309;
t27 = (t207 + t236) * pkin(5) + t58;
t26 = -qJ(6) * t237 + t30;
t25 = -mrSges(6,1) * t69 + mrSges(6,2) * t68;
t21 = -pkin(5) * t159 - t155 * t241 + t29;
t13 = -pkin(5) * t69 + t32;
t8 = -qJD(5) * t30 + t199;
t7 = -t228 * t90 + t223;
t6 = -qJ(6) * t207 + (-qJD(5) * t90 + t168) * t153 + t223;
t3 = pkin(5) * t126 + t168 * t155 + (-t84 + (-t85 + t241) * t153) * qJD(5) + t199;
t12 = [t21 * t41 + t29 * t42 + t26 * t43 + t30 * t44 + t55 * t24 + t27 * t56 + t6 * t71 + t7 * t72 + t3 * t73 + t8 * t74 + t89 * t25 + t59 * t109 + t137 * t195 + t248 * t58 + (-t116 * t90 + t117 * t89) * mrSges(5,3) + m(5) * (t33 * t90 - t48 * t58 + t49 * t59 + t269) + m(7) * (t1 * t21 + t11 * t6 + t13 * t55 + t2 * t26 + t27 * t28 + t3 * t9) + m(6) * (t19 * t8 + t20 * t7 + t29 * t5 + t30 * t4 + t45 * t58 + t269) + (-t49 * mrSges(5,3) + t35 / 0.2e1 + t36 / 0.2e1 + t216 * t309 + t218 * t102 + t217 * t101 + t290 + t287) * t126 + (-t48 * mrSges(5,3) + t251 / 0.2e1 + t304 * t277 + t303 * t274 + t121 / 0.2e1 + t87 / 0.2e1 + t170 * mrSges(6,3) + t291 + (-t11 * t153 - t155 * t9) * mrSges(7,3) - t308) * t127 + (t204 * t111 + (-t146 * t139 + t129 / 0.2e1 + 0.3e1 / 0.2e1 * t149 + t245 / 0.2e1 - t204 * t118 + 0.2e1 * t140 * mrSges(4,2)) * qJD(3)) * t156 - (-t33 * mrSges(5,3) + t66 / 0.2e1 + t64 / 0.2e1 + t114 / 0.2e1 + t67 / 0.2e1 + t65 / 0.2e1 + t115 / 0.2e1 - Ifges(5,4) * t117 + t217 * t69 + t218 * t68 + (Ifges(5,2) + t216) * t116 + t158) * t159 + (t13 * t183 + Ifges(5,1) * t117 - Ifges(5,4) * t116 + (mrSges(5,3) + t185) * t32 + t190 * mrSges(7,3) + t188 * mrSges(6,3) + (mrSges(6,3) * t169 - mrSges(7,3) * t187 + t184 * t28 + t186 * t45 + t312 * t286 + t311 * t284 + t313 * t281 - t304 * t155 / 0.2e1) * qJD(5) + t295 * t289 + t296 * t288 + t297 * t282 + t305 * t274 + (qJD(5) * t303 + t306) * t277) * t135 + (t204 * t112 + (-t146 * t141 - t128 / 0.2e1 + t140 * mrSges(4,1) - t243 / 0.2e1 - t204 * t119 + (t210 * mrSges(4,1) - 0.3e1 / 0.2e1 * t264 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t156) * qJD(1) + (t91 + qJD(1) * (-mrSges(5,1) * t159 + mrSges(5,2) * t135) + m(5) * (t123 + t233)) * pkin(3)) * qJD(3)) * t154; -(t24 + t25 + t255) * t159 + t222 * t126 + (-t154 * t139 + t156 * t141 + (-t154 ^ 2 - t156 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + (-t153 * t265 + t155 * t266 + t109) * t127 + m(5) * (-t126 * t48 + t127 * t49 - t249) + m(7) * (t11 * t238 + t126 * t28 - t13 * t159 - t236 * t9) + m(6) * (t126 * t45 - t19 * t236 + t20 * t238 - t249) + m(4) * (t111 * t154 - t112 * t156 + (-t118 * t154 + t119 * t156) * qJD(3)) + (-t256 + t267 * t155 + t268 * t153 + (-t153 * t266 - t155 * t265) * qJD(5) + m(5) * t33 + m(7) * t292 + m(6) * t310) * t135; (-t250 + t36 + t35) * t278 + (Ifges(5,1) * t278 + t297 * t281 + t295 * t284 + t296 * t286 + t298 + t308) * t124 + (t121 + t87) * t279 - m(6) * (t19 * t22 + t20 * t23) + ((t246 * t33 - t247 * t32) * pkin(3) - t123 * t221 + t48 * t50 - t49 * t51) * m(5) + (t215 + t139) * t119 + (t19 * t239 + t20 * t240 + t310) * mrSges(6,3) + t311 * t289 + t312 * t288 + t313 * t282 - (-Ifges(4,2) * t232 + t129 + t149) * t231 / 0.2e1 + t305 * t153 / 0.2e1 + t306 * t274 + (m(6) * t145 * t170 + t291) * qJD(5) + (-m(6) * t45 - t248) * t50 + t299 * t56 + t300 * t71 + t301 * t73 + (-t1 * t130 + t11 * t300 + t13 * t136 + t132 * t2 + t28 * t299 + t301 * t9) * m(7) + t224 * Ifges(4,5) * t156 / 0.2e1 + (t214 - t141) * t118 + (m(6) * t147 - mrSges(5,1) - t186) * t32 + (t11 * t240 + t239 * t9 + t292) * mrSges(7,3) + (-t140 * (mrSges(4,1) * t154 + mrSges(4,2) * t156) - (Ifges(4,1) * t156 - t264) * t232 / 0.2e1) * qJD(1) - (-t319 * t286 - t315 * t284 + (-Ifges(7,3) - Ifges(6,3)) * t281 + t290 + t316) * t125 + (m(6) * t189 - t153 * t42 + t155 * t44 - t227 * t74 - t228 * t72) * t145 + t49 * t252 - t209 * t255 - t208 * t256 - t13 * t184 - t33 * mrSges(5,2) + t48 * t254 - t23 * t72 - t22 * t74 - t51 * t109 - t111 * mrSges(4,2) - t112 * mrSges(4,1) - Ifges(5,6) * t116 + Ifges(5,5) * t117 - t130 * t41 + t132 * t43 + t136 * t24 + t147 * t25 + (t227 / 0.2e1 - t239 / 0.2e1) * t303 + (-t228 / 0.2e1 + t240 / 0.2e1) * t304 - Ifges(4,6) * t203 / 0.2e1 - t91 * t221 + t128 * t232 / 0.2e1; -t124 * t109 - t222 * t125 + (t266 * t309 - t268) * t155 + (-t265 * t309 + t267) * t153 + t195 + (-t125 * t28 + t187 * t309 - t190) * m(7) + (-t125 * t45 - t169 * t309 - t188) * m(6) + (-t124 * t49 + t125 * t48 + t191) * m(5); (-t102 * t56 + t41) * pkin(5) + (t101 * t9 + t102 * t11) * mrSges(7,3) + (t101 * t19 + t102 * t20) * mrSges(6,3) + t158 + t115 + t114 + t66 + t67 + t65 + t64 + (-(t10 - t9) * t11 + (-t102 * t28 + t1) * pkin(5)) * m(7) - t10 * t71 - t19 * t72 + t11 * t73 + t20 * t74 - t28 * (mrSges(7,1) * t102 + mrSges(7,2) * t101) - t45 * (mrSges(6,1) * t102 + mrSges(6,2) * t101) + (t101 * t322 - t317) * t284 + t304 * t283 + (t101 * t315 - t102 * t319) * t281 + (-t102 * t320 + t303 + t318) * t286; -t101 * t71 + t102 * t73 + 0.2e1 * (t13 / 0.2e1 + t11 * t286 + t9 * t283) * m(7) + t24;];
tauc  = t12(:);
