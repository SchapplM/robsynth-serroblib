% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:13
% EndTime: 2019-03-08 19:51:33
% DurationCPUTime: 11.76s
% Computational Cost: add. (3206->544), mult. (6589->721), div. (0->0), fcn. (4263->10), ass. (0->256)
t362 = -mrSges(5,1) + mrSges(6,2);
t322 = -m(5) - m(4);
t350 = Ifges(6,4) - Ifges(5,5);
t349 = -Ifges(5,6) + Ifges(6,5);
t354 = m(7) + m(6);
t355 = m(7) * pkin(9);
t360 = pkin(4) * t354 + t355 - t362;
t144 = sin(qJ(6));
t147 = cos(qJ(6));
t265 = qJD(4) * t144;
t145 = sin(qJ(4));
t269 = qJD(2) * t145;
t101 = t147 * t269 - t265;
t148 = cos(qJ(4));
t253 = qJD(2) * qJD(4);
t108 = qJDD(2) * t145 + t148 * t253;
t31 = qJD(6) * t101 + qJDD(4) * t147 + t108 * t144;
t107 = -t148 * qJDD(2) + t145 * t253;
t99 = qJDD(6) - t107;
t22 = mrSges(7,1) * t99 - mrSges(7,3) * t31;
t263 = qJD(4) * t147;
t102 = t144 * t269 + t263;
t32 = -qJD(6) * t102 - qJDD(4) * t144 + t108 * t147;
t23 = -mrSges(7,2) * t99 + mrSges(7,3) * t32;
t254 = t148 * qJD(2);
t130 = qJD(6) + t254;
t61 = -mrSges(7,2) * t130 + mrSges(7,3) * t101;
t62 = mrSges(7,1) * t130 - mrSges(7,3) * t102;
t339 = -t144 * t62 + t147 * t61;
t331 = -qJD(6) * t339 - t144 * t23 - t147 * t22;
t239 = mrSges(5,3) * t254;
t241 = mrSges(6,1) * t254;
t273 = qJD(4) * t362 + t239 + t241;
t200 = t147 * mrSges(7,1) - t144 * mrSges(7,2);
t312 = mrSges(3,1) - mrSges(4,2);
t361 = t200 + mrSges(6,1) + mrSges(5,3) + t312;
t321 = pkin(4) + pkin(9);
t143 = cos(pkin(6));
t270 = qJD(1) * t143;
t231 = t145 * t270;
t151 = -pkin(2) - pkin(8);
t149 = cos(qJ(2));
t142 = sin(pkin(6));
t271 = qJD(1) * t142;
t232 = t149 * t271;
t186 = qJD(3) - t232;
t89 = qJD(2) * t151 + t186;
t54 = -t148 * t89 + t231;
t162 = pkin(5) * t254 + t54;
t357 = t162 + qJD(5);
t24 = -qJD(4) * t321 + t357;
t209 = -qJ(5) * t148 + qJ(3);
t317 = pkin(9) * t145;
t179 = t209 + t317;
t146 = sin(qJ(2));
t233 = t146 * t271;
t203 = pkin(4) * t269 + t233;
t48 = qJD(2) * t179 + t203;
t10 = t144 * t24 + t147 * t48;
t250 = qJDD(2) * qJ(3);
t251 = qJDD(1) * t146;
t168 = t142 * t251 + t250;
t185 = qJD(3) + t232;
t260 = qJD(5) * t148;
t153 = qJ(5) * t107 + (t185 - t260) * qJD(2) + t168;
t16 = t108 * t321 + t153;
t222 = qJD(4) * t270;
t252 = qJDD(1) * t143;
t264 = qJD(4) * t145;
t268 = qJD(2) * t146;
t229 = t142 * t268;
t122 = qJD(1) * t229;
t280 = t142 * t149;
t76 = qJDD(1) * t280 - t122;
t176 = qJDD(3) - t76;
t56 = qJDD(2) * t151 + t176;
t18 = -t145 * t252 - t89 * t264 + (-t222 + t56) * t148;
t169 = qJDD(5) - t18;
t5 = -pkin(5) * t107 - qJDD(4) * t321 + t169;
t9 = -t144 * t48 + t147 * t24;
t1 = qJD(6) * t9 + t144 * t5 + t147 * t16;
t2 = -qJD(6) * t10 - t144 * t16 + t147 * t5;
t204 = t1 * t144 + t147 * t2;
t257 = qJD(6) * t147;
t259 = qJD(6) * t144;
t359 = -t10 * t257 + t9 * t259 - t204;
t358 = t354 - t322;
t199 = mrSges(7,1) * t144 + mrSges(7,2) * t147;
t166 = m(7) * qJ(5) + t199;
t201 = mrSges(5,1) * t145 + mrSges(5,2) * t148;
t351 = mrSges(3,2) - mrSges(4,3);
t356 = t166 * t148 - t201 + t351;
t325 = t31 / 0.2e1;
t324 = t32 / 0.2e1;
t323 = t99 / 0.2e1;
t275 = t147 * t148;
t310 = pkin(5) - t151;
t111 = t310 * t148;
t139 = t145 * pkin(4);
t95 = t139 + t179;
t36 = t111 * t144 + t147 * t95;
t262 = qJD(4) * t148;
t224 = pkin(4) * t262 + qJ(5) * t264 + qJD(3);
t58 = (qJD(4) * pkin(9) - qJD(5)) * t148 + t224;
t215 = qJD(4) * t310;
t97 = t145 * t215;
t353 = -qJD(6) * t36 - t144 * t58 - t147 * t97 - (-t144 * t149 - t146 * t275) * t271;
t278 = t144 * t148;
t35 = t111 * t147 - t144 * t95;
t352 = qJD(6) * t35 - t144 * t97 + t147 * t58 - (-t146 * t278 + t147 * t149) * t271;
t81 = -t107 * mrSges(6,1) + qJDD(4) * mrSges(6,2);
t348 = qJDD(4) * mrSges(5,1) + mrSges(5,3) * t107 - t81;
t11 = -mrSges(7,1) * t32 + mrSges(7,2) * t31;
t80 = mrSges(6,1) * t108 - qJDD(4) * mrSges(6,3);
t347 = -t80 + t11;
t308 = Ifges(5,4) * t145;
t197 = t148 * Ifges(5,1) - t308;
t346 = Ifges(5,5) * qJD(4) + t102 * Ifges(7,5) + t101 * Ifges(7,6) + t130 * Ifges(7,3) + qJD(2) * t197;
t242 = mrSges(6,1) * t269;
t118 = -qJD(4) * mrSges(6,3) + t242;
t41 = -mrSges(7,1) * t101 + mrSges(7,2) * t102;
t289 = t118 - t41;
t15 = -qJDD(4) * pkin(4) + t169;
t266 = qJD(4) * qJ(5);
t230 = t148 * t270;
t55 = t145 * t89 + t230;
t40 = -t55 - t266;
t345 = -qJD(4) * t40 - t15;
t240 = mrSges(5,3) * t269;
t116 = -qJD(4) * mrSges(5,2) - t240;
t344 = t116 - t118;
t303 = Ifges(6,6) * t148;
t343 = t145 * (-Ifges(5,2) * t148 - t308) + t148 * (Ifges(6,2) * t145 + t303);
t342 = t145 * t350 + t148 * t349;
t245 = t145 * t56 + t148 * t252 + t89 * t262;
t17 = -t145 * t222 + t245;
t337 = t145 * t17 + t148 * t18;
t14 = -qJDD(4) * qJ(5) + qJD(4) * (-qJD(5) + t231) - t245;
t336 = -t14 * t145 - t148 * t15;
t187 = t145 * Ifges(6,3) - t303;
t302 = t102 * Ifges(7,4);
t26 = t101 * Ifges(7,2) + t130 * Ifges(7,6) + t302;
t96 = Ifges(7,4) * t101;
t27 = t102 * Ifges(7,1) + t130 * Ifges(7,5) + t96;
t334 = Ifges(6,5) * qJD(4) + qJD(2) * t187 + t144 * t27 + t147 * t26;
t333 = -m(6) * qJ(5) - mrSges(6,3);
t329 = -m(7) * (-pkin(5) - pkin(8)) + t361;
t328 = -mrSges(5,2) + t166 - t333;
t141 = sin(pkin(10));
t284 = t141 * t142;
t288 = cos(pkin(10));
t213 = t288 * t146;
t282 = t141 * t149;
t85 = t143 * t282 + t213;
t42 = t145 * t284 - t85 * t148;
t214 = t142 * t288;
t212 = t288 * t149;
t283 = t141 * t146;
t83 = -t143 * t212 + t283;
t44 = t145 * t214 + t83 * t148;
t87 = t143 * t145 + t148 * t280;
t170 = -g(1) * t42 + g(2) * t44 - g(3) * t87;
t297 = t145 * mrSges(6,2);
t198 = -t148 * mrSges(6,3) - t297;
t327 = -m(6) * t209 - m(7) * t317 - t145 * mrSges(7,3) - t198 + t356 + (-m(7) + t322) * qJ(3);
t152 = qJD(2) ^ 2;
t326 = Ifges(7,1) * t325 + Ifges(7,4) * t324 + Ifges(7,5) * t323;
t319 = t102 / 0.2e1;
t307 = Ifges(5,4) * t148;
t306 = Ifges(7,4) * t144;
t305 = Ifges(7,4) * t147;
t304 = Ifges(6,6) * t145;
t34 = t230 + (-pkin(5) * qJD(2) + t89) * t145;
t30 = t34 + t266;
t287 = qJD(4) * t30;
t106 = qJD(2) * qJ(3) + t233;
t285 = t106 * t149;
t281 = t142 * t146;
t279 = t144 * t145;
t277 = t145 * t147;
t272 = pkin(2) * t280 + qJ(3) * t281;
t103 = pkin(4) * t254 + qJ(5) * t269;
t267 = qJD(2) * t149;
t261 = qJD(4) * t151;
t258 = qJD(6) * t145;
t255 = qJDD(2) * mrSges(4,2);
t79 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t108;
t247 = t79 + t347;
t246 = Ifges(7,5) * t31 + Ifges(7,6) * t32 + Ifges(7,3) * t99;
t243 = t116 - t289;
t237 = t145 * t280;
t236 = pkin(8) * t280 + t272;
t74 = t83 * pkin(2);
t235 = -pkin(8) * t83 - t74;
t75 = t85 * pkin(2);
t234 = -pkin(8) * t85 - t75;
t228 = t142 * t267;
t223 = qJD(1) * t267;
t217 = -t259 / 0.2e1;
t211 = -t253 / 0.2e1;
t205 = t142 * t223;
t202 = t10 * t147 - t144 * t9;
t196 = Ifges(7,1) * t147 - t306;
t195 = Ifges(7,1) * t144 + t305;
t194 = -t145 * Ifges(5,2) + t307;
t192 = -Ifges(7,2) * t144 + t305;
t191 = Ifges(7,2) * t147 + t306;
t189 = Ifges(7,5) * t147 - Ifges(7,6) * t144;
t188 = Ifges(7,5) * t144 + Ifges(7,6) * t147;
t57 = qJD(2) * t185 + t168;
t181 = qJ(3) * t57 + qJD(3) * t106;
t104 = t198 * qJD(2);
t178 = -t104 - t339;
t50 = t144 * t87 + t147 * t281;
t49 = -t144 * t281 + t147 * t87;
t177 = t106 * t267 + t146 * t57;
t63 = qJD(2) * t209 + t203;
t175 = t63 * (-mrSges(6,2) * t148 + mrSges(6,3) * t145);
t174 = t106 * (mrSges(5,1) * t148 - mrSges(5,2) * t145);
t172 = t148 * (-Ifges(5,1) * t145 - t307);
t165 = -t144 * t258 + t147 * t262;
t164 = t144 * t262 + t145 * t257;
t159 = -Ifges(7,5) * t145 + t148 * t195;
t158 = -Ifges(7,6) * t145 + t148 * t191;
t157 = -Ifges(7,3) * t145 + t148 * t188;
t155 = qJD(6) * t202 + t204;
t133 = Ifges(6,6) * t269;
t110 = t310 * t145;
t109 = t139 + t209;
t105 = t201 * qJD(2);
t100 = -qJD(2) * pkin(2) + t186;
t98 = t148 * t215;
t93 = Ifges(6,4) * qJD(4) - Ifges(6,2) * t254 + t133;
t90 = Ifges(5,6) * qJD(4) + qJD(2) * t194;
t88 = t143 * t148 - t237;
t86 = -t143 * t283 + t212;
t84 = t143 * t213 + t282;
t77 = (t223 + t251) * t142;
t73 = pkin(9) * t254 + t103;
t72 = t224 - t260;
t66 = t86 * t139;
t65 = t84 * t139;
t64 = -qJDD(2) * pkin(2) + t176;
t53 = mrSges(5,1) * t108 - mrSges(5,2) * t107;
t52 = -mrSges(6,2) * t108 + mrSges(6,3) * t107;
t47 = -qJD(4) * t237 + t143 * t262 - t148 * t229;
t46 = qJD(4) * t87 - t145 * t229;
t37 = -qJD(4) * pkin(4) + qJD(5) + t54;
t21 = pkin(4) * t108 + t153;
t20 = t144 * t34 + t147 * t73;
t19 = -t144 * t73 + t147 * t34;
t13 = qJD(6) * t49 + t144 * t47 + t147 * t228;
t12 = -qJD(6) * t50 - t144 * t228 + t147 * t47;
t6 = -pkin(5) * t108 - t14;
t3 = t31 * Ifges(7,4) + t32 * Ifges(7,2) + t99 * Ifges(7,6);
t4 = [t12 * t62 + t13 * t61 + t49 * t22 + t50 * t23 - t348 * t87 + t273 * t47 + t247 * t88 - t243 * t46 + (-m(2) - m(3) - t358) * g(3) + m(7) * (t1 * t50 + t10 * t13 + t12 * t9 + t2 * t49 - t30 * t46 + t6 * t88) + m(5) * (t17 * t88 - t18 * t87 - t46 * t55 + t47 * t54) + m(6) * (-t14 * t88 + t15 * t87 + t37 * t47 + t40 * t46) + ((-t351 * t152 + t312 * qJDD(2) + (t104 + t105) * qJD(2)) * t149 + (-qJDD(2) * t351 - t152 * t312 + t52 + t53) * t146 + m(4) * (t100 * t268 - t149 * t64 + t177) + m(3) * (t146 * t77 + t149 * t76) + m(5) * t177 + m(6) * (t146 * t21 + t267 * t63)) * t142 + (m(2) + 0.2e1 * (m(4) / 0.2e1 + m(3) / 0.2e1) * t143 ^ 2) * qJDD(1); (t145 * t349 - t148 * t350) * qJDD(4) / 0.2e1 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + t342 * qJD(4) ^ 2 / 0.2e1 + t343 * t211 + t186 * t105 + (m(4) * t74 - m(7) * (t65 - t74) - m(6) * (t235 + t65) - m(5) * t235 + t327 * t84 + t329 * t83) * g(2) + (m(4) * t75 - m(7) * (t66 - t75) - m(6) * (t234 + t66) - m(5) * t234 + t327 * t86 + t329 * t85) * g(1) + (t1 * t36 + t10 * t352 - t110 * t6 + t2 * t35 - t30 * t98 + t353 * t9) * m(7) + t145 * (Ifges(6,5) * qJDD(4) + Ifges(6,6) * t107 + Ifges(6,3) * t108) / 0.2e1 + (-m(6) * (t149 * t63 + (-t145 * t40 - t148 * t37) * t146) - m(5) * (t285 + (t145 * t55 - t148 * t54) * t146) - m(4) * (t100 * t146 + t285)) * t271 - t145 * (-Ifges(5,4) * t107 - Ifges(5,2) * t108 + Ifges(5,6) * qJDD(4)) / 0.2e1 + t109 * t52 - t110 * t11 - t98 * t41 + qJ(3) * t53 + t35 * t22 + t36 * t23 + (-t346 / 0.2e1 - t37 * mrSges(6,1) - t54 * mrSges(5,3) + t93 / 0.2e1 - t9 * mrSges(7,1) + t10 * mrSges(7,2)) * t264 + qJD(4) * t174 + qJD(4) * t175 + (t72 - t232) * t104 + t145 * t26 * t217 - t6 * t200 * t145 + (t64 - t122) * mrSges(4,2) + (-t77 + t205) * mrSges(3,2) + m(6) * (t109 * t21 + t63 * t72) + (t145 * (Ifges(6,3) * t148 + t304) + t172) * t253 / 0.2e1 + (qJD(6) * t27 + t3) * t277 / 0.2e1 + (-Ifges(5,1) * t107 - Ifges(5,4) * t108 + Ifges(5,5) * qJDD(4) + t246) * t148 / 0.2e1 + (t40 * mrSges(6,1) - t55 * mrSges(5,3) + t334 / 0.2e1 - t90 / 0.2e1) * t262 + m(5) * t181 + t10 * mrSges(7,3) * t165 - t9 * mrSges(7,3) * t164 + m(4) * (-pkin(2) * t64 + t181) + t30 * (-mrSges(7,1) * t165 + mrSges(7,2) * t164) + (-m(4) * t272 - m(5) * t236 - t354 * (t281 * t139 + t236) + ((-m(7) * pkin(5) - t361) * t149 + (t297 - t333 * t148 - (mrSges(7,3) + t355) * t145 + t356) * t146) * t142) * g(3) + t273 * (t145 * t261 + t148 * t233) + (-m(7) * t30 - t344 - t41) * t145 * t233 + t344 * t148 * t261 + (m(6) * ((t145 * t37 - t148 * t40) * qJD(4) + t336) + m(5) * ((t145 * t54 + t148 * t55) * qJD(4) + t337) + (t79 - t80) * t145 + t348 * t148) * t151 - t336 * mrSges(6,1) - t337 * mrSges(5,3) + t2 * (mrSges(7,1) * t148 - mrSges(7,3) * t279) + t1 * (-mrSges(7,2) * t148 + mrSges(7,3) * t277) + t130 * (qJD(4) * t157 + t189 * t258) / 0.2e1 + t101 * (qJD(4) * t158 + t192 * t258) / 0.2e1 - pkin(2) * t255 + t352 * t61 + t353 * t62 + t108 * t187 / 0.2e1 - t108 * t194 / 0.2e1 - t107 * t197 / 0.2e1 + t21 * t198 + t57 * t201 - t148 * (Ifges(6,4) * qJDD(4) + Ifges(6,6) * t108) / 0.2e1 + t107 * t304 / 0.2e1 - t107 * t148 * Ifges(6,2) + (qJD(2) * qJD(3) - t205 + t250 + t57) * mrSges(4,3) + (t76 + t122) * mrSges(3,1) + (qJD(4) * t159 + t196 * t258) * t319 + (Ifges(7,3) * t148 + t145 * t188) * t323 + (Ifges(7,6) * t148 + t145 * t191) * t324 + (Ifges(7,5) * t148 + t145 * t195) * t325 + t279 * t326; t255 - t152 * mrSges(4,3) + m(4) * t64 + ((t144 * t61 + t147 * t62 + t273) * qJD(4) + m(7) * (t10 * t265 + t263 * t9 + t6) + m(5) * (qJD(4) * t54 + t17) + m(6) * (qJD(4) * t37 - t14) + t247) * t145 + (t243 * qJD(4) + m(7) * (t287 + t359) + m(5) * (qJD(4) * t55 + t18) + m(6) * t345 + t348 + t331) * t148 + (-m(6) * t63 - m(7) * t202 + t106 * t322 - t105 + t178) * qJD(2) + t358 * (-g(1) * t85 - g(2) * t83 + g(3) * t280); t349 * t108 + t350 * t107 + t130 * t30 * t200 + t346 * t269 / 0.2e1 + t347 * qJ(5) - t289 * qJD(5) + t342 * t211 + (-m(6) * t37 + t239 - t273) * t55 + (-m(6) * t40 + t240 + t344) * t54 - t334 * t254 / 0.2e1 + t162 * t41 - t144 * t3 / 0.2e1 - t103 * t104 - pkin(4) * t81 - t20 * t61 - t19 * t62 + t37 * t242 - t17 * mrSges(5,2) + t18 * mrSges(5,1) - t14 * mrSges(6,3) + t15 * mrSges(6,2) + (-pkin(4) * t15 - qJ(5) * t14 - qJD(5) * t40 - t103 * t63) * m(6) + (Ifges(6,1) + Ifges(5,3)) * qJDD(4) + (-t10 * (mrSges(7,2) * t145 + mrSges(7,3) * t275) - t9 * (-mrSges(7,1) * t145 - mrSges(7,3) * t278) - t174 - t175) * qJD(2) - (Ifges(6,3) * t254 + t133 + t93) * t269 / 0.2e1 - (t101 * t191 + t102 * t195 + t130 * t188) * qJD(6) / 0.2e1 - (t101 * t158 + t102 * t159 + t130 * t157) * qJD(2) / 0.2e1 + (-t170 + t359) * mrSges(7,3) + (-t328 * t88 + t360 * t87) * g(3) + (-t328 * (t145 * t85 + t148 * t284) + t360 * t42) * g(1) + (qJ(5) * t6 - t10 * t20 - t19 * t9 + t30 * t357) * m(7) + (t328 * (-t83 * t145 + t148 * t214) - t360 * t44) * g(2) + (-m(7) * t155 + t331) * t321 + (t343 / 0.2e1 - t172 / 0.2e1) * t152 - t26 * t257 / 0.2e1 + t90 * t254 / 0.2e1 - t40 * t241 + t6 * t199 + t27 * t217 + t189 * t323 + t192 * t324 + t196 * t325 + t147 * t326; t289 * qJD(4) - t178 * t254 + t81 + (t202 * t254 + t155 + t170 - t287) * m(7) + (t254 * t63 + t170 - t345) * m(6) - t331; -t1 * mrSges(7,2) + t2 * mrSges(7,1) - t30 * (mrSges(7,1) * t102 + mrSges(7,2) * t101) - t102 * (Ifges(7,1) * t101 - t302) / 0.2e1 + t26 * t319 - t130 * (Ifges(7,5) * t101 - Ifges(7,6) * t102) / 0.2e1 - t9 * t61 + t10 * t62 - g(1) * ((-t144 * t86 + t147 * t42) * mrSges(7,1) + (-t144 * t42 - t147 * t86) * mrSges(7,2)) - g(2) * ((-t144 * t84 - t147 * t44) * mrSges(7,1) + (t144 * t44 - t147 * t84) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t49 - mrSges(7,2) * t50) + (t10 * t102 + t101 * t9) * mrSges(7,3) + t246 - (-Ifges(7,2) * t102 + t27 + t96) * t101 / 0.2e1;];
tau  = t4;
