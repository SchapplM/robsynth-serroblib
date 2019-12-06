% Calculate vector of inverse dynamics joint torques for
% S5RPRPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:55:53
% EndTime: 2019-12-05 17:56:22
% DurationCPUTime: 10.27s
% Computational Cost: add. (4906->497), mult. (12066->686), div. (0->0), fcn. (8654->14), ass. (0->264)
t213 = cos(pkin(8));
t278 = qJD(1) * t213;
t188 = qJD(3) - t278;
t271 = qJDD(1) * qJ(2);
t273 = qJD(1) * qJD(2);
t183 = t271 + t273;
t354 = t183 + t273;
t353 = Ifges(4,3) + Ifges(5,3);
t210 = sin(pkin(9));
t212 = cos(pkin(9));
t216 = sin(qJ(3));
t219 = cos(qJ(3));
t166 = t210 * t219 + t212 * t216;
t153 = t166 * qJD(3);
t232 = qJD(1) * t166;
t335 = t213 * t232 - t153;
t236 = t210 * t216 - t212 * t219;
t334 = t188 * t236;
t211 = sin(pkin(8));
t207 = t211 ^ 2;
t352 = t354 * t207;
t318 = m(5) + m(6);
t330 = t318 + m(3) + m(4);
t209 = qJ(3) + pkin(9);
t200 = sin(t209);
t307 = pkin(3) * t216;
t170 = pkin(4) * t200 + t307;
t351 = -m(6) * t170 - qJ(2) * t330 + mrSges(2,2) - mrSges(3,3);
t322 = m(5) * pkin(3);
t269 = qJDD(1) * t213;
t186 = qJDD(3) - t269;
t312 = t186 / 0.2e1;
t349 = m(5) * t307;
t272 = qJD(1) * qJD(3);
t147 = (qJDD(1) * t219 - t216 * t272) * t211;
t148 = (-qJDD(1) * t216 - t219 * t272) * t211;
t82 = t147 * t212 + t148 * t210;
t348 = Ifges(5,5) * t82;
t215 = sin(qJ(5));
t218 = cos(qJ(5));
t279 = qJD(1) * t211;
t256 = t219 * t279;
t257 = t216 * t279;
t123 = t210 * t257 - t212 * t256;
t124 = t211 * t232;
t248 = t123 * t215 - t218 * t124;
t81 = -t147 * t210 + t148 * t212;
t21 = qJD(5) * t248 + t215 * t81 + t218 * t82;
t347 = Ifges(6,5) * t21;
t346 = Ifges(5,6) * t81;
t64 = t123 * t218 + t124 * t215;
t22 = qJD(5) * t64 - t215 * t82 + t218 * t81;
t345 = Ifges(6,6) * t22;
t94 = -t166 * t215 - t218 * t236;
t344 = qJD(5) * t94 + t215 * t335 - t218 * t334;
t95 = t166 * t218 - t215 * t236;
t343 = -qJD(5) * t95 + t215 * t334 + t218 * t335;
t342 = Ifges(4,5) * t147;
t341 = Ifges(4,6) * t148;
t178 = qJDD(5) + t186;
t340 = Ifges(6,3) * t178;
t339 = mrSges(4,1) + t322;
t308 = pkin(3) * t212;
t196 = pkin(4) + t308;
t309 = pkin(3) * t210;
t150 = t196 * t215 + t218 * t309;
t305 = pkin(7) * t124;
t174 = -pkin(2) * t213 - pkin(6) * t211 - pkin(1);
t151 = qJD(1) * t174 + qJD(2);
t260 = qJ(2) * t278;
t105 = t151 * t216 + t219 * t260;
t89 = -qJ(4) * t257 + t105;
t297 = t212 * t89;
t132 = t219 * t151;
t292 = t211 * t219;
t263 = qJ(4) * t292;
t291 = t213 * t216;
t264 = qJ(2) * t291;
t230 = -t263 - t264;
t88 = qJD(1) * t230 + t132;
t45 = -t210 * t88 - t297;
t29 = t45 + t305;
t306 = pkin(7) * t123;
t83 = t210 * t89;
t46 = t212 * t88 - t83;
t30 = t46 + t306;
t338 = -t150 * qJD(5) + t215 * t30 - t218 * t29;
t149 = t196 * t218 - t215 * t309;
t337 = t149 * qJD(5) - t215 * t29 - t218 * t30;
t217 = sin(qJ(1));
t220 = cos(qJ(1));
t336 = t211 * (g(2) * t217 - g(3) * t220);
t289 = t213 * t219;
t187 = qJ(2) * t289;
t122 = t216 * t174 + t187;
t208 = t213 ^ 2;
t333 = t207 + t208;
t332 = t186 * t353 + t341 + t342 + t346 + t348;
t331 = 0.2e1 * t312;
t329 = t208 * t354;
t299 = Ifges(4,4) * t219;
t300 = Ifges(4,4) * t216;
t327 = (qJ(2) * (mrSges(4,1) * t219 - mrSges(4,2) * t216) - t216 * (-Ifges(4,2) * t219 - t300) / 0.2e1 + t219 * (-Ifges(4,1) * t216 - t299) / 0.2e1) * t207;
t146 = qJDD(1) * t174 + qJDD(2);
t131 = t219 * t146;
t259 = qJ(2) * qJD(3) * t213;
t274 = qJD(4) * t211;
t225 = qJD(1) * (-t259 - t274);
t36 = pkin(3) * t186 - qJ(4) * t147 + t131 + (-qJD(3) * t151 - t183 * t213) * t216 + t219 * t225;
t275 = qJD(3) * t219;
t262 = t216 * t146 + t151 * t275 + t183 * t289;
t41 = qJ(4) * t148 + t216 * t225 + t262;
t12 = -t210 * t41 + t212 * t36;
t6 = pkin(4) * t186 - pkin(7) * t82 + t12;
t13 = t210 * t36 + t212 * t41;
t7 = pkin(7) * t81 + t13;
t75 = pkin(3) * t188 + t88;
t42 = t212 * t75 - t83;
t25 = pkin(4) * t188 + t306 + t42;
t43 = t210 * t75 + t297;
t26 = t43 - t305;
t8 = -t215 * t26 + t218 * t25;
t2 = qJD(5) * t8 + t215 * t6 + t218 * t7;
t9 = t215 * t25 + t218 * t26;
t3 = -qJD(5) * t9 - t215 * t7 + t218 * t6;
t326 = -t3 * mrSges(6,1) + t2 * mrSges(6,2);
t201 = cos(t209);
t205 = t219 * pkin(3);
t171 = pkin(4) * t201 + t205;
t214 = -qJ(4) - pkin(6);
t241 = -mrSges(3,1) * t213 + mrSges(3,2) * t211;
t325 = -m(6) * (-(pkin(2) + t171) * t213 - pkin(1)) - m(4) * t174 + m(3) * pkin(1) - t241 + mrSges(2,1) - m(5) * (-(t205 + pkin(2)) * t213 - pkin(1)) + (-m(6) * (-pkin(7) + t214) + mrSges(6,3) + mrSges(4,3) - m(5) * t214 + mrSges(5,3)) * t211;
t233 = (-Ifges(4,2) * t216 + t299) * t211;
t234 = (t219 * Ifges(4,1) - t300) * t211;
t324 = t219 * (t188 * Ifges(4,6) + qJD(1) * t233) + t216 * (t188 * Ifges(4,5) + qJD(1) * t234);
t244 = t216 * t259;
t52 = -qJD(1) * t244 + t262;
t53 = -qJD(3) * t105 - t183 * t291 + t131;
t323 = -t53 * mrSges(4,1) - t12 * mrSges(5,1) + t52 * mrSges(4,2) + t13 * mrSges(5,2);
t221 = qJD(1) ^ 2;
t321 = -t248 / 0.2e1;
t320 = t64 / 0.2e1;
t319 = -t64 / 0.2e1;
t317 = t8 * mrSges(6,3);
t316 = t9 * mrSges(6,3);
t315 = -t123 / 0.2e1;
t180 = qJD(5) + t188;
t313 = -t180 / 0.2e1;
t310 = Ifges(6,4) * t64;
t304 = g(1) * t211;
t277 = qJD(2) * t213;
t281 = t174 * t275 + t219 * t277;
t72 = qJD(3) * t230 - t216 * t274 + t281;
t255 = t216 * t277;
t73 = -t255 - t219 * t274 + (-t187 + (qJ(4) * t211 - t174) * t216) * qJD(3);
t35 = t210 * t73 + t212 * t72;
t293 = t211 * t216;
t106 = -qJ(4) * t293 + t122;
t164 = t219 * t174;
t93 = -t263 + t164 + (-qJ(2) * t216 - pkin(3)) * t213;
t51 = t212 * t106 + t210 * t93;
t302 = mrSges(5,3) * t123;
t301 = mrSges(5,3) * t124;
t298 = Ifges(5,4) * t123;
t296 = qJ(2) * t221;
t202 = t211 * qJ(2);
t172 = t211 * t183;
t290 = t213 * t217;
t288 = t213 * t220;
t287 = t216 * t217;
t286 = t216 * t220;
t285 = t217 * t219;
t284 = t219 * t220;
t203 = qJ(5) + t209;
t194 = sin(t203);
t195 = cos(t203);
t117 = t194 * t290 + t195 * t220;
t118 = -t194 * t220 + t195 * t290;
t283 = t117 * mrSges(6,1) + t118 * mrSges(6,2);
t119 = t194 * t288 - t195 * t217;
t120 = -t194 * t217 - t195 * t288;
t282 = -t119 * mrSges(6,1) + t120 * mrSges(6,2);
t280 = t352 * qJ(2);
t199 = t211 * qJD(2);
t254 = t211 * t275;
t161 = pkin(3) * t254 + t199;
t167 = pkin(3) * t293 + t202;
t276 = qJD(3) * t211;
t270 = qJDD(1) * t211;
t268 = t340 + t345 + t347;
t266 = mrSges(4,3) * t293;
t152 = pkin(3) * t257 + qJ(2) * t279 + qJD(4);
t253 = -t81 * mrSges(5,1) + t82 * mrSges(5,2);
t252 = -t22 * mrSges(6,1) + t21 * mrSges(6,2);
t34 = -t210 * t72 + t212 * t73;
t50 = -t106 * t210 + t212 * t93;
t247 = pkin(3) * t256;
t246 = mrSges(4,3) * t257;
t245 = mrSges(4,3) * t256;
t100 = -pkin(3) * t148 + qJDD(4) + t172;
t242 = -mrSges(3,1) * t269 + mrSges(3,2) * t270;
t240 = mrSges(4,1) * t216 + mrSges(4,2) * t219;
t239 = -mrSges(6,1) * t194 - mrSges(6,2) * t195;
t142 = t236 * t211;
t40 = -pkin(4) * t213 + pkin(7) * t142 + t50;
t141 = t166 * t211;
t44 = -pkin(7) * t141 + t51;
t14 = -t215 * t44 + t218 * t40;
t15 = t215 * t40 + t218 * t44;
t76 = -t141 * t218 + t142 * t215;
t77 = -t141 * t215 - t142 * t218;
t143 = -mrSges(4,2) * t188 - t246;
t144 = mrSges(4,1) * t188 - t245;
t237 = t143 * t219 - t144 * t216;
t235 = t268 - t326;
t157 = t213 * t286 - t285;
t155 = t213 * t287 + t284;
t228 = t188 * t211 * (-Ifges(4,5) * t216 - Ifges(4,6) * t219);
t198 = -qJDD(1) * pkin(1) + qJDD(2);
t191 = t207 * t296;
t159 = t240 * t211;
t158 = -t213 * t284 - t287;
t156 = t213 * t285 - t286;
t145 = t240 * t279;
t138 = -t200 * t217 - t201 * t288;
t137 = t200 * t288 - t201 * t217;
t136 = -t200 * t220 + t201 * t290;
t135 = t200 * t290 + t201 * t220;
t128 = t211 * t153;
t126 = t236 * t276;
t121 = t164 - t264;
t116 = Ifges(5,4) * t124;
t108 = -mrSges(4,2) * t186 + mrSges(4,3) * t148;
t107 = mrSges(4,1) * t186 - mrSges(4,3) * t147;
t104 = -t216 * t260 + t132;
t103 = -qJD(3) * t122 - t255;
t102 = -t244 + t281;
t101 = -pkin(4) * t123 + t247;
t98 = mrSges(5,1) * t188 + t302;
t97 = -mrSges(5,2) * t188 - t301;
t96 = pkin(4) * t141 + t167;
t90 = -pkin(4) * t126 + t161;
t87 = pkin(4) * t124 + t152;
t74 = mrSges(5,1) * t124 - mrSges(5,2) * t123;
t60 = Ifges(6,4) * t248;
t59 = mrSges(5,1) * t186 - mrSges(5,3) * t82;
t58 = -mrSges(5,2) * t186 + mrSges(5,3) * t81;
t57 = -t123 * Ifges(5,1) + t188 * Ifges(5,5) - t116;
t56 = -t124 * Ifges(5,2) + t188 * Ifges(5,6) - t298;
t55 = mrSges(6,1) * t180 + mrSges(6,3) * t64;
t54 = -mrSges(6,2) * t180 + mrSges(6,3) * t248;
t47 = -pkin(4) * t81 + t100;
t39 = -qJD(5) * t77 + t126 * t218 + t128 * t215;
t38 = qJD(5) * t76 + t126 * t215 - t128 * t218;
t33 = -mrSges(6,1) * t248 - mrSges(6,2) * t64;
t28 = -Ifges(6,1) * t64 + Ifges(6,5) * t180 + t60;
t27 = Ifges(6,2) * t248 + Ifges(6,6) * t180 - t310;
t24 = pkin(7) * t126 + t35;
t23 = pkin(7) * t128 + t34;
t17 = -mrSges(6,2) * t178 + mrSges(6,3) * t22;
t16 = mrSges(6,1) * t178 - mrSges(6,3) * t21;
t5 = -qJD(5) * t15 - t215 * t24 + t218 * t23;
t4 = qJD(5) * t14 + t215 * t23 + t218 * t24;
t1 = [(-t158 * mrSges(4,1) - t138 * mrSges(5,1) - t120 * mrSges(6,1) - t157 * mrSges(4,2) - t137 * mrSges(5,2) - t119 * mrSges(6,2) + t325 * t220 + (t349 - t351) * t217) * g(2) + (t156 * mrSges(4,1) + t136 * mrSges(5,1) + t118 * mrSges(6,1) - t155 * mrSges(4,2) - t135 * mrSges(5,2) - t117 * mrSges(6,2) + t217 * t325 + t351 * t220 - t286 * t322) * g(3) + (Ifges(3,2) * t269 + Ifges(3,4) * t270 - t341 / 0.2e1 - t342 / 0.2e1 - t347 / 0.2e1 - t345 / 0.2e1 - t346 / 0.2e1 - t348 / 0.2e1 - t340 / 0.2e1 - t353 * t312 + t323 + t326) * t213 + (-Ifges(5,1) * t128 + Ifges(5,4) * t126) * t315 + (t126 * t43 + t128 * t42) * mrSges(5,3) - t124 * (-Ifges(5,4) * t128 + Ifges(5,2) * t126) / 0.2e1 + t152 * (-mrSges(5,1) * t126 - mrSges(5,2) * t128) + t188 * (-Ifges(5,5) * t128 + Ifges(5,6) * t126) / 0.2e1 + (t104 * t266 + t228 / 0.2e1) * qJD(3) + (-t105 * t254 - t292 * t53) * mrSges(4,3) + t148 * t233 / 0.2e1 + t39 * t316 + (Ifges(6,1) * t38 + Ifges(6,4) * t39) * t319 + t327 * t272 + m(3) * (-pkin(1) * t198 + qJ(2) * t329 + t280) - t52 * t266 + (t271 * t333 + t329 + t352) * mrSges(3,3) + m(5) * (t100 * t167 + t12 * t50 + t13 * t51 + t152 * t161 + t34 * t42 + t35 * t43) + m(6) * (t14 * t3 + t15 * t2 + t4 * t9 + t47 * t96 + t5 * t8 + t87 * t90) - (t100 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,1) * t82 + Ifges(5,4) * t81 + Ifges(5,5) * t331) * t142 - (-t100 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t82 + Ifges(5,2) * t81 + Ifges(5,6) * t331) * t141 + t248 * (Ifges(6,4) * t38 + Ifges(6,2) * t39) / 0.2e1 + (-mrSges(4,1) * t148 + mrSges(4,2) * t147) * t202 + (-mrSges(6,1) * t47 + mrSges(6,3) * t2 + Ifges(6,4) * t21 + Ifges(6,2) * t22 + Ifges(6,6) * t178) * t76 - (t268 + t332) * t213 / 0.2e1 - t324 * t276 / 0.2e1 + (mrSges(6,2) * t47 - mrSges(6,3) * t3 + Ifges(6,1) * t21 + Ifges(6,4) * t22 + Ifges(6,5) * t178) * t77 + t145 * t199 + t14 * t16 + t15 * t17 + t198 * t241 + t38 * t28 / 0.2e1 + t39 * t27 / 0.2e1 - pkin(1) * t242 + t4 * t54 + t5 * t55 + t51 * t58 + t50 * t59 + ((Ifges(4,5) * t219 - Ifges(4,6) * t216) * t312 + Ifges(3,4) * t269 + Ifges(3,1) * t270) * t211 + t87 * (-mrSges(6,1) * t39 + mrSges(6,2) * t38) + t90 * t33 + t35 * t97 + t34 * t98 + t147 * t234 / 0.2e1 + t121 * t107 + t122 * t108 + t126 * t56 / 0.2e1 - t128 * t57 / 0.2e1 + t102 * t143 + t103 * t144 + t161 * t74 + t180 * (Ifges(6,5) * t38 + Ifges(6,6) * t39) / 0.2e1 + t96 * t252 + t167 * t253 + m(4) * (t102 * t105 + t103 * t104 + t121 * t53 + t122 * t52 + t280) + (Ifges(4,1) * t147 + Ifges(4,4) * t148 + Ifges(4,5) * t186) * t292 / 0.2e1 - (Ifges(4,4) * t147 + Ifges(4,2) * t148 + Ifges(4,6) * t186) * t293 / 0.2e1 + t159 * t172 + Ifges(2,3) * qJDD(1) - t38 * t317; t335 * t98 - t334 * t97 + t343 * t55 + t344 * t54 - t333 * t221 * mrSges(3,3) + t242 + t237 * qJD(3) + (-t237 * t213 + (-t145 - t33 - t74) * t211) * qJD(1) + t94 * t16 + t95 * t17 - t236 * t59 + t166 * t58 + t216 * t108 + t219 * t107 - (g(2) * t220 + g(3) * t217) * t330 + (t2 * t95 - t87 * t279 + t3 * t94 + t343 * t8 + t344 * t9) * m(6) + (-t12 * t236 + t13 * t166 - t152 * t279 - t334 * t43 + t335 * t42) * m(5) + (t52 * t216 + t53 * t219 - t191 + t188 * (-t104 * t216 + t105 * t219)) * m(4) + (-t208 * t296 - t191 + t198) * m(3); (Ifges(5,2) * t123 - t116 + t57) * t124 / 0.2e1 + (Ifges(6,6) * t313 + Ifges(6,4) * t320 + Ifges(6,2) * t321 - t27 / 0.2e1 + t87 * mrSges(6,1) - t316) * t64 - t152 * (-mrSges(5,1) * t123 - mrSges(5,2) * t124) - t188 * (-Ifges(5,5) * t124 + Ifges(5,6) * t123) / 0.2e1 + t123 * (-Ifges(5,1) * t124 + t298) / 0.2e1 + t58 * t309 + t56 * t315 + (t12 * t212 + t13 * t210) * t322 - t327 * t221 - t323 + (Ifges(6,5) * t313 + t317 + Ifges(6,1) * t320 + Ifges(6,4) * t321 - t28 / 0.2e1 - t87 * mrSges(6,2)) * t248 + (t245 + t144) * t105 + t235 - t42 * t301 + t59 * t308 + t332 + t337 * t54 + t338 * t55 + (-t101 * t87 + t149 * t3 + t150 * t2 + t170 * t304 + t337 * t9 + t338 * t8) * m(6) + (-mrSges(4,2) * t158 - m(6) * (-t170 * t288 + t171 * t217) - t282 + t137 * mrSges(5,1) - t138 * mrSges(5,2) + t339 * t157) * g(3) + (-mrSges(4,2) * t156 - m(6) * (t170 * t290 + t171 * t220) - t283 - t135 * mrSges(5,1) - t136 * mrSges(5,2) - t339 * t155) * g(2) + t324 * t279 / 0.2e1 - qJD(1) * t228 / 0.2e1 + (mrSges(5,1) * t200 + mrSges(5,2) * t201 - t239 + t349) * t304 + (-t143 - t246) * t104 - t46 * t97 - t45 * t98 - t101 * t33 + t149 * t16 + t150 * t17 + g(1) * t159 - t74 * t247 - m(5) * (t152 * t247 + t42 * t45 + t43 * t46) - t43 * t302; t318 * t213 * g(1) - t123 * t98 + t124 * t97 - t248 * t54 - t64 * t55 + t252 + t253 + (-t248 * t9 - t8 * t64 + t336 + t47) * m(6) + (-t123 * t42 + t124 * t43 + t100 + t336) * m(5); -t87 * (-mrSges(6,1) * t64 + mrSges(6,2) * t248) + (Ifges(6,1) * t248 + t310) * t320 + t27 * t319 + (Ifges(6,5) * t248 + Ifges(6,6) * t64) * t313 - t8 * t54 + t9 * t55 - g(2) * t283 - g(3) * t282 - t239 * t304 + (t248 * t8 - t64 * t9) * mrSges(6,3) + t235 + (Ifges(6,2) * t64 + t28 + t60) * t321;];
tau = t1;
