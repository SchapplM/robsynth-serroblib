% Calculate vector of inverse dynamics joint torques for
% S5RRPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:36
% EndTime: 2019-12-31 19:32:04
% DurationCPUTime: 14.95s
% Computational Cost: add. (5829->588), mult. (13485->796), div. (0->0), fcn. (9786->14), ass. (0->269)
t196 = sin(pkin(9));
t198 = cos(pkin(9));
t227 = -mrSges(5,1) * t198 + mrSges(5,2) * t196;
t356 = -m(5) * pkin(3) - mrSges(4,1) + t227;
t308 = m(5) + m(6);
t355 = m(4) + t308;
t354 = mrSges(4,2) - mrSges(6,3);
t202 = sin(qJ(2));
t205 = cos(qJ(2));
t178 = -mrSges(3,1) * t205 + mrSges(3,2) * t202;
t195 = qJ(2) + pkin(8);
t191 = sin(t195);
t193 = cos(t195);
t353 = t191 * t354 + t356 * t193 + t178;
t352 = m(5) * qJ(4) + mrSges(5,3);
t203 = sin(qJ(1));
t206 = cos(qJ(1));
t351 = g(1) * t206 + g(2) * t203;
t350 = -m(3) * pkin(1) - mrSges(2,1) + t353;
t197 = sin(pkin(8));
t271 = cos(pkin(8));
t166 = t197 * t205 + t202 * t271;
t149 = t166 * qJD(1);
t132 = qJD(2) * t196 + t149 * t198;
t201 = sin(qJ(5));
t204 = cos(qJ(5));
t232 = t198 * qJD(2) - t149 * t196;
t348 = -t132 * t201 + t204 * t232;
t73 = t132 * t204 + t201 * t232;
t250 = qJD(1) * qJD(2);
t239 = t202 * t250;
t249 = qJDD(1) * t205;
t171 = -t239 + t249;
t172 = qJDD(1) * t202 + t205 * t250;
t126 = t197 * t171 + t172 * t271;
t96 = qJDD(2) * t198 - t126 * t196;
t97 = qJDD(2) * t196 + t126 * t198;
t23 = qJD(5) * t348 + t201 * t96 + t204 * t97;
t319 = t23 / 0.2e1;
t24 = -qJD(5) * t73 - t201 * t97 + t204 * t96;
t318 = t24 / 0.2e1;
t310 = t96 / 0.2e1;
t309 = t97 / 0.2e1;
t235 = t271 * t205;
t256 = qJD(1) * t202;
t147 = -qJD(1) * t235 + t197 * t256;
t199 = -qJ(3) - pkin(6);
t177 = t199 * t202;
t169 = qJD(1) * t177;
t161 = qJD(2) * pkin(2) + t169;
t179 = t199 * t205;
t170 = qJD(1) * t179;
t236 = t271 * t170;
t118 = t197 * t161 - t236;
t106 = qJD(2) * qJ(4) + t118;
t293 = pkin(2) * t205;
t186 = pkin(1) + t293;
t174 = -qJD(1) * t186 + qJD(3);
t82 = pkin(3) * t147 - qJ(4) * t149 + t174;
t44 = -t106 * t196 + t198 * t82;
t29 = pkin(4) * t147 - pkin(7) * t132 + t44;
t45 = t198 * t106 + t196 * t82;
t34 = pkin(7) * t232 + t45;
t9 = -t201 * t34 + t204 * t29;
t347 = t9 * mrSges(6,1);
t125 = -t271 * t171 + t172 * t197;
t122 = qJDD(5) + t125;
t307 = t122 / 0.2e1;
t306 = t125 / 0.2e1;
t346 = t171 / 0.2e1;
t10 = t201 * t29 + t204 * t34;
t345 = t10 * mrSges(6,2);
t344 = t44 * mrSges(5,1);
t343 = t45 * mrSges(5,2);
t287 = qJD(2) / 0.2e1;
t295 = pkin(2) * t197;
t182 = qJ(4) + t295;
t286 = pkin(7) + t182;
t158 = t286 * t196;
t159 = t286 * t198;
t108 = -t158 * t201 + t159 * t204;
t167 = t196 * t204 + t198 * t201;
t267 = t147 * t198;
t153 = t197 * t170;
t124 = t169 * t271 + t153;
t246 = pkin(2) * t256;
t93 = pkin(3) * t149 + qJ(4) * t147 + t246;
t52 = -t124 * t196 + t198 * t93;
t33 = pkin(4) * t149 + pkin(7) * t267 + t52;
t268 = t147 * t196;
t53 = t198 * t124 + t196 * t93;
t43 = pkin(7) * t268 + t53;
t342 = -qJD(4) * t167 - qJD(5) * t108 + t201 * t43 - t204 * t33;
t107 = -t158 * t204 - t159 * t201;
t218 = t196 * t201 - t198 * t204;
t341 = -qJD(4) * t218 + qJD(5) * t107 - t201 * t33 - t204 * t43;
t144 = qJD(5) + t147;
t338 = t232 * Ifges(5,6);
t339 = t132 * Ifges(5,5);
t340 = t73 * Ifges(6,5) + Ifges(6,6) * t348 + t147 * Ifges(5,3) + t144 * Ifges(6,3) + t338 + t339;
t226 = mrSges(5,1) * t196 + mrSges(5,2) * t198;
t117 = t161 * t271 + t153;
t99 = -qJD(2) * pkin(3) + qJD(4) - t117;
t337 = t99 * t226;
t272 = qJDD(2) / 0.2e1;
t280 = t149 * mrSges(4,3);
t336 = qJD(2) * mrSges(4,1) + mrSges(5,1) * t232 - mrSges(5,2) * t132 - t280;
t152 = t167 * qJD(5);
t88 = t167 * t147;
t335 = -t88 - t152;
t151 = t218 * qJD(5);
t89 = t218 * t147;
t334 = -t89 - t151;
t333 = t191 * t351;
t187 = pkin(6) * t249;
t162 = -pkin(6) * t239 + t187;
t163 = t172 * pkin(6);
t330 = t162 * t205 + t163 * t202;
t83 = -mrSges(5,2) * t147 + mrSges(5,3) * t232;
t84 = mrSges(5,1) * t147 - mrSges(5,3) * t132;
t329 = -t196 * t84 + t198 * t83;
t269 = qJDD(1) * pkin(1);
t142 = -pkin(2) * t171 + qJDD(3) - t269;
t42 = pkin(3) * t125 - qJ(4) * t126 - qJD(4) * t149 + t142;
t253 = qJD(3) * t202;
t113 = qJDD(2) * pkin(2) - qJ(3) * t172 - qJD(1) * t253 - t163;
t254 = qJD(2) * t202;
t244 = pkin(6) * t254;
t252 = qJD(3) * t205;
t120 = qJ(3) * t171 + t187 + (-t244 + t252) * qJD(1);
t58 = t197 * t113 + t271 * t120;
t51 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t58;
t16 = -t196 * t51 + t198 * t42;
t17 = t196 * t42 + t198 * t51;
t328 = -t16 * t196 + t17 * t198;
t327 = 0.2e1 * t272;
t326 = -m(3) * pkin(6) + m(5) * t199 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - t226;
t11 = pkin(7) * t96 + t17;
t8 = pkin(4) * t125 - pkin(7) * t97 + t16;
t1 = qJD(5) * t9 + t11 * t204 + t201 * t8;
t2 = -qJD(5) * t10 - t11 * t201 + t204 * t8;
t323 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t321 = Ifges(6,4) * t319 + Ifges(6,2) * t318 + Ifges(6,6) * t307;
t320 = Ifges(6,1) * t319 + Ifges(6,4) * t318 + Ifges(6,5) * t307;
t296 = Ifges(6,4) * t73;
t27 = Ifges(6,2) * t348 + Ifges(6,6) * t144 + t296;
t317 = t27 / 0.2e1;
t70 = Ifges(6,4) * t348;
t28 = Ifges(6,1) * t73 + Ifges(6,5) * t144 + t70;
t316 = t28 / 0.2e1;
t315 = Ifges(5,1) * t309 + Ifges(5,4) * t310 + Ifges(5,5) * t306;
t314 = -t348 / 0.2e1;
t313 = t348 / 0.2e1;
t312 = -t73 / 0.2e1;
t311 = t73 / 0.2e1;
t305 = -t144 / 0.2e1;
t304 = t144 / 0.2e1;
t303 = -t147 / 0.2e1;
t302 = t147 / 0.2e1;
t299 = t149 / 0.2e1;
t297 = t198 / 0.2e1;
t292 = pkin(4) * t196;
t291 = pkin(6) * t205;
t288 = t198 * pkin(4);
t148 = t166 * qJD(2);
t215 = -t197 * t202 + t235;
t150 = t215 * qJD(2);
t245 = pkin(2) * t254;
t75 = pkin(3) * t148 - qJ(4) * t150 - qJD(4) * t166 + t245;
t237 = qJD(2) * t199;
t145 = t202 * t237 + t252;
t146 = t205 * t237 - t253;
t95 = t145 * t271 + t197 * t146;
t38 = t196 * t75 + t198 * t95;
t285 = Ifges(3,4) * t202;
t284 = Ifges(3,4) * t205;
t283 = Ifges(5,4) * t196;
t282 = Ifges(5,4) * t198;
t281 = t117 * mrSges(4,3);
t279 = t149 * Ifges(4,4);
t276 = t191 * mrSges(5,3);
t270 = qJ(4) * t191;
t266 = t150 * t196;
t265 = t150 * t198;
t262 = t166 * t196;
t261 = t166 * t198;
t200 = -pkin(7) - qJ(4);
t260 = t191 * t200;
t259 = t193 * t206;
t194 = pkin(9) + qJ(5);
t190 = sin(t194);
t258 = t203 * t190;
t192 = cos(t194);
t257 = t203 * t192;
t116 = -pkin(3) * t215 - qJ(4) * t166 - t186;
t129 = t197 * t177 - t179 * t271;
t65 = t196 * t116 + t198 * t129;
t255 = qJD(1) * t205;
t248 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t256) * t291;
t247 = Ifges(6,5) * t23 + Ifges(6,6) * t24 + Ifges(6,3) * t122;
t242 = -t196 * (t132 * Ifges(5,4) + Ifges(5,2) * t232 + Ifges(5,6) * t147) / 0.2e1;
t241 = (t132 * Ifges(5,1) + Ifges(5,4) * t232 + Ifges(5,5) * t147) * t297;
t240 = t271 * pkin(2);
t50 = -t96 * mrSges(5,1) + t97 * mrSges(5,2);
t7 = -t24 * mrSges(6,1) + t23 * mrSges(6,2);
t37 = -t196 * t95 + t198 * t75;
t233 = t125 * mrSges(4,1) + t126 * mrSges(4,2);
t64 = t198 * t116 - t129 * t196;
t94 = t145 * t197 - t271 * t146;
t123 = t169 * t197 - t236;
t128 = -t271 * t177 - t179 * t197;
t185 = -t240 - pkin(3);
t229 = mrSges(3,1) * t202 + mrSges(3,2) * t205;
t225 = Ifges(5,1) * t198 - t283;
t224 = t205 * Ifges(3,2) + t285;
t223 = -Ifges(5,2) * t196 + t282;
t222 = Ifges(3,5) * t205 - Ifges(3,6) * t202;
t221 = Ifges(5,5) * t198 - Ifges(5,6) * t196;
t220 = t196 * t44 - t198 * t45;
t41 = -pkin(4) * t215 - pkin(7) * t261 + t64;
t46 = -pkin(7) * t262 + t65;
t14 = -t201 * t46 + t204 * t41;
t15 = t201 * t41 + t204 * t46;
t57 = t113 * t271 - t197 * t120;
t184 = pkin(3) + t288;
t219 = t184 * t193 - t260;
t217 = pkin(1) * t229;
t216 = t202 * (Ifges(3,1) * t205 - t285);
t211 = m(6) * t184 + t192 * mrSges(6,1) - t190 * mrSges(6,2);
t54 = -qJDD(2) * pkin(3) + qJDD(4) - t57;
t188 = Ifges(3,4) * t255;
t180 = t206 * t186;
t176 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t255;
t173 = t185 - t288;
t157 = Ifges(3,1) * t256 + Ifges(3,5) * qJD(2) + t188;
t156 = Ifges(3,6) * qJD(2) + qJD(1) * t224;
t143 = Ifges(4,4) * t147;
t138 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t147;
t137 = t192 * t259 + t258;
t136 = -t190 * t259 + t257;
t135 = t190 * t206 - t193 * t257;
t134 = t192 * t206 + t193 * t258;
t114 = mrSges(4,1) * t147 + mrSges(4,2) * t149;
t110 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t126;
t109 = -qJDD(2) * mrSges(4,2) - mrSges(4,3) * t125;
t103 = t149 * Ifges(4,1) + Ifges(4,5) * qJD(2) - t143;
t102 = -t147 * Ifges(4,2) + Ifges(4,6) * qJD(2) + t279;
t101 = t218 * t166;
t100 = t167 * t166;
t90 = pkin(4) * t262 + t128;
t79 = -pkin(4) * t268 + t123;
t67 = pkin(4) * t266 + t94;
t66 = -pkin(4) * t232 + t99;
t60 = mrSges(5,1) * t125 - mrSges(5,3) * t97;
t59 = -mrSges(5,2) * t125 + mrSges(5,3) * t96;
t56 = mrSges(6,1) * t144 - mrSges(6,3) * t73;
t55 = -mrSges(6,2) * t144 + mrSges(6,3) * t348;
t49 = -t150 * t167 + t151 * t166;
t48 = -t150 * t218 - t152 * t166;
t35 = t97 * Ifges(5,4) + t96 * Ifges(5,2) + t125 * Ifges(5,6);
t32 = -mrSges(6,1) * t348 + mrSges(6,2) * t73;
t31 = -t96 * pkin(4) + t54;
t30 = -pkin(7) * t266 + t38;
t25 = pkin(4) * t148 - pkin(7) * t265 + t37;
t19 = -mrSges(6,2) * t122 + mrSges(6,3) * t24;
t18 = mrSges(6,1) * t122 - mrSges(6,3) * t23;
t4 = -qJD(5) * t15 - t201 * t30 + t204 * t25;
t3 = qJD(5) * t14 + t201 * t25 + t204 * t30;
t5 = [t205 * t157 * t287 + (-m(4) * t57 + m(5) * t54 - t110 + t50) * t128 + (t171 * t291 + t330) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t330) + (-Ifges(4,4) * t299 + Ifges(5,3) * t302 - Ifges(4,2) * t303 + Ifges(6,3) * t304 - Ifges(4,6) * t287 + t338 / 0.2e1 + t339 / 0.2e1 + Ifges(6,5) * t311 + Ifges(6,6) * t313 - t118 * mrSges(4,3) + t344 - t343 + t347 - t345 - t102 / 0.2e1 + t174 * mrSges(4,1) + t340 / 0.2e1) * t148 + (-t16 * t261 - t17 * t262 - t265 * t44 - t266 * t45) * mrSges(5,3) + (t142 * mrSges(4,2) - t57 * mrSges(4,3) + Ifges(4,1) * t126 - Ifges(4,4) * t125 + Ifges(4,5) * t327 + t221 * t306 + t223 * t310 + t225 * t309 + t54 * t226) * t166 - (t142 * mrSges(4,1) + t16 * mrSges(5,1) - t17 * mrSges(5,2) - t58 * mrSges(4,3) - Ifges(4,4) * t126 + Ifges(5,5) * t309 + Ifges(6,5) * t319 + Ifges(4,2) * t125 - t327 * Ifges(4,6) + Ifges(5,6) * t310 + Ifges(6,6) * t318 + Ifges(5,3) * t306 + Ifges(6,3) * t307 + t323) * t215 - (Ifges(5,5) * t97 + Ifges(5,6) * t96 + Ifges(5,3) * t125 + t247) * t215 / 0.2e1 + t172 * t284 / 0.2e1 + (-pkin(6) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t172) + Ifges(3,1) * t172 + Ifges(3,4) * t346 + t327 * Ifges(3,5)) * t202 + (-m(5) * t180 - t137 * mrSges(6,1) - t136 * mrSges(6,2) + (-m(4) - m(6)) * (-t203 * t199 + t180) + (-m(6) * t292 + t326) * t203 + (-m(6) * t219 - t191 * t352 + t350) * t206) * g(2) + (Ifges(6,4) * t48 + Ifges(6,2) * t49) * t313 + (-m(4) * t117 + m(5) * t99 - t336) * t94 + (t337 + Ifges(4,1) * t299 + t221 * t302 + Ifges(4,4) * t303 + t241 + t242 + Ifges(4,5) * t287 + t232 * t223 / 0.2e1 + t132 * t225 / 0.2e1 + t103 / 0.2e1 - t281 + t174 * mrSges(4,2)) * t150 + t224 * t346 + (t222 * t287 - t248) * qJD(2) + (Ifges(6,5) * t48 + Ifges(6,6) * t49) * t304 + (Ifges(6,1) * t48 + Ifges(6,4) * t49) * t311 + Ifges(3,6) * t205 * t272 + m(4) * (t118 * t95 + t129 * t58 - t142 * t186 + t174 * t245) + m(5) * (t16 * t64 + t17 * t65 + t37 * t44 + t38 * t45) + (t216 + t205 * (-Ifges(3,2) * t202 + t284)) * t250 / 0.2e1 + (-t1 * t100 + t10 * t49 + t101 * t2 - t48 * t9) * mrSges(6,3) + (-Ifges(6,4) * t101 - Ifges(6,2) * t100) * t318 + (-Ifges(6,5) * t101 - Ifges(6,6) * t100) * t307 + (-Ifges(6,1) * t101 - Ifges(6,4) * t100) * t319 + t31 * (mrSges(6,1) * t100 - mrSges(6,2) * t101) + t114 * t245 + t261 * t315 + t48 * t316 + (-t135 * mrSges(6,1) - t134 * mrSges(6,2) + (m(4) * t199 - m(6) * (-t199 + t292) + t326) * t206 + (m(4) * t186 - m(5) * (-t186 - t270) + t276 - m(6) * (-t186 - t219) - t350) * t203) * g(1) + t49 * t317 - t101 * t320 - t100 * t321 - t186 * t233 + t14 * t18 + t15 * t19 + t3 * t55 + t4 * t56 + t64 * t60 + t65 * t59 + t66 * (-mrSges(6,1) * t49 + mrSges(6,2) * t48) + t67 * t32 - qJDD(2) * mrSges(3,2) * t291 - t176 * t244 + t38 * t83 + t37 * t84 + t90 * t7 - t217 * t250 - t156 * t254 / 0.2e1 - t35 * t262 / 0.2e1 + t129 * t109 - t178 * t269 + t95 * t138 - pkin(1) * (-mrSges(3,1) * t171 + mrSges(3,2) * t172) + m(6) * (t1 * t15 + t10 * t3 + t14 * t2 + t31 * t90 + t4 * t9 + t66 * t67) + t205 * (Ifges(3,4) * t172 + Ifges(3,2) * t171 + Ifges(3,6) * qJDD(2)) / 0.2e1 + Ifges(2,3) * qJDD(1); t109 * t295 + t35 * t297 + t102 * t299 + (Ifges(5,3) * t149 - t147 * t221) * t303 + (Ifges(6,5) * t89 + Ifges(6,6) * t88 + Ifges(6,3) * t149) * t305 + (Ifges(5,5) * t196 + Ifges(5,6) * t198) * t306 + (Ifges(5,1) * t196 + t282) * t309 + (-m(4) * t293 - m(6) * (-t260 + t293) - t211 * t193 - m(5) * (t270 + t293) - t276 + t353) * g(3) + (pkin(6) * t176 + t156 / 0.2e1) * t256 + (t248 + (t217 - t216 / 0.2e1) * qJD(1)) * qJD(1) + t147 * t241 + (-t267 * t44 - t268 * t45 + t328) * mrSges(5,3) + (-t220 * qJD(4) - t123 * t99 + t182 * t328 + t185 * t54 - t44 * t52 - t45 * t53) * m(5) + t329 * qJD(4) + (-Ifges(6,5) * t151 - Ifges(6,6) * t152) * t304 + (-Ifges(6,1) * t151 - Ifges(6,4) * t152) * t311 + (-Ifges(6,4) * t151 - Ifges(6,2) * t152) * t313 - t149 * t344 + t118 * t280 + t110 * t240 + (Ifges(6,5) * t167 - Ifges(6,6) * t218) * t307 + (Ifges(6,4) * t167 - Ifges(6,2) * t218) * t318 + (Ifges(6,1) * t167 - Ifges(6,4) * t218) * t319 + t31 * (mrSges(6,1) * t218 + mrSges(6,2) * t167) + (-t1 * t218 + t10 * t335 - t167 * t2 - t334 * t9) * mrSges(6,3) - t218 * t321 - t149 * t347 + t147 * t242 + (-mrSges(6,1) * t335 + mrSges(6,2) * t334) * t66 + t336 * t123 + t147 * t337 + t149 * t343 + t149 * t345 - (-Ifges(4,1) * t147 - t279 + t340) * t149 / 0.2e1 + t341 * t55 + t342 * t56 + (t1 * t108 + t10 * t341 + t107 * t2 + t173 * t31 + t342 * t9 - t66 * t79) * m(6) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + ((t197 * t58 + t271 * t57) * pkin(2) + t117 * t123 - t118 * t124 - t174 * t246) * m(4) + (-Ifges(4,2) * t149 + t103 - t143) * t302 - (-Ifges(3,2) * t256 + t157 + t188) * t255 / 0.2e1 - t232 * (Ifges(5,6) * t149 - t147 * t223) / 0.2e1 + (-t196 * t60 + t198 * t59) * t182 - t132 * (Ifges(5,5) * t149 - t147 * t225) / 0.2e1 + t54 * t227 + (Ifges(5,2) * t198 + t283) * t310 + (Ifges(6,1) * t89 + Ifges(6,4) * t88 + Ifges(6,5) * t149) * t312 + (Ifges(6,4) * t89 + Ifges(6,2) * t88 + Ifges(6,6) * t149) * t314 + t196 * t315 - t151 * t316 + t351 * (t229 + t355 * pkin(2) * t202 + (m(6) * t200 - t352 + t354) * t193 + (t211 - t356) * t191) - t152 * t317 + t167 * t320 + t57 * mrSges(4,1) - t58 * mrSges(4,2) - t114 * t246 - t79 * t32 - t53 * t83 - t52 * t84 - t88 * t27 / 0.2e1 - t89 * t28 / 0.2e1 - t222 * t250 / 0.2e1 + t107 * t18 + t108 * t19 - Ifges(4,6) * t125 + Ifges(4,5) * t126 - t124 * t138 - qJD(2) * (-Ifges(4,5) * t147 - Ifges(4,6) * t149) / 0.2e1 - t147 * t281 - t162 * mrSges(3,2) - t163 * mrSges(3,1) + Ifges(3,6) * t171 + Ifges(3,5) * t172 + t173 * t7 - t174 * (mrSges(4,1) * t149 - mrSges(4,2) * t147) + t185 * t50; -t218 * t18 + t167 * t19 + t196 * t59 + t198 * t60 + t335 * t56 + t334 * t55 + (-t32 + t336) * t149 + (t138 + t329) * t147 + t233 + (-g(1) * t203 + g(2) * t206) * t355 + (t1 * t167 + t10 * t334 - t149 * t66 - t2 * t218 + t335 * t9) * m(6) + (-t147 * t220 - t149 * t99 + t16 * t198 + t17 * t196) * m(5) + (t117 * t149 + t118 * t147 + t142) * m(4); t308 * t193 * g(3) + t132 * t84 - t232 * t83 - t348 * t55 + t73 * t56 + t50 + t7 + (-t10 * t348 + t73 * t9 + t31 - t333) * m(6) + (t132 * t44 - t232 * t45 - t333 + t54) * m(5); -t66 * (mrSges(6,1) * t73 + mrSges(6,2) * t348) + (Ifges(6,1) * t348 - t296) * t312 + t27 * t311 + (Ifges(6,5) * t348 - Ifges(6,6) * t73) * t305 - t9 * t55 + t10 * t56 - g(1) * (mrSges(6,1) * t136 - mrSges(6,2) * t137) - g(2) * (-mrSges(6,1) * t134 + mrSges(6,2) * t135) - g(3) * (-mrSges(6,1) * t190 - mrSges(6,2) * t192) * t191 + (t10 * t73 + t348 * t9) * mrSges(6,3) + t247 + (-Ifges(6,2) * t73 + t28 + t70) * t314 + t323;];
tau = t5;
