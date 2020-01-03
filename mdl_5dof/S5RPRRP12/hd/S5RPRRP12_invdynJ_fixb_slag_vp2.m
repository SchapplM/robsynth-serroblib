% Calculate vector of inverse dynamics joint torques for
% S5RPRRP12
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP12_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:16
% EndTime: 2019-12-31 18:56:39
% DurationCPUTime: 14.10s
% Computational Cost: add. (2630->436), mult. (5261->564), div. (0->0), fcn. (2871->6), ass. (0->205)
t338 = Ifges(5,4) + Ifges(6,4);
t347 = -Ifges(4,2) / 0.2e1;
t339 = Ifges(5,1) + Ifges(6,1);
t321 = Ifges(5,5) + Ifges(6,5);
t337 = Ifges(5,2) + Ifges(6,2);
t320 = Ifges(5,6) + Ifges(6,6);
t124 = sin(qJ(4));
t127 = cos(qJ(4));
t230 = qJD(3) * t127;
t128 = cos(qJ(3));
t233 = qJD(1) * t128;
t94 = -t124 * t233 + t230;
t346 = -t94 / 0.2e1;
t232 = qJD(3) * t124;
t95 = t127 * t233 + t232;
t345 = -t95 / 0.2e1;
t125 = sin(qJ(3));
t222 = t125 * qJD(1);
t115 = qJD(4) + t222;
t344 = -t115 / 0.2e1;
t252 = Ifges(4,4) * t128;
t343 = t125 * t347 + t252 / 0.2e1;
t342 = t338 * t127;
t341 = t338 * t124;
t220 = qJD(1) * qJD(3);
t101 = qJDD(1) * t128 - t125 * t220;
t40 = qJD(4) * t94 + qJDD(3) * t124 + t101 * t127;
t284 = t40 / 0.2e1;
t41 = -qJD(4) * t95 + qJDD(3) * t127 - t101 * t124;
t283 = t41 / 0.2e1;
t102 = -qJDD(1) * t125 - t128 * t220;
t91 = qJDD(4) - t102;
t282 = t91 / 0.2e1;
t116 = pkin(4) * t127 + pkin(3);
t340 = m(6) * t116;
t336 = Ifges(5,3) + Ifges(6,3);
t253 = Ifges(4,4) * t125;
t165 = t128 * Ifges(4,1) - t253;
t329 = t338 * t94;
t297 = t321 * t115 + t339 * t95 + t329;
t335 = Ifges(4,5) * qJD(3) + qJD(1) * t165 + t127 * t297;
t313 = -t320 * t124 + t321 * t127;
t311 = -t337 * t124 + t342;
t309 = t339 * t127 - t341;
t334 = t336 * t344 + t320 * t346 + t321 * t345 + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t343;
t123 = -qJ(5) - pkin(7);
t333 = -m(6) * t123 + mrSges(5,3) + mrSges(6,3);
t226 = qJD(4) * t127;
t227 = qJD(4) * t124;
t221 = qJD(1) * qJD(2);
t114 = qJDD(1) * qJ(2) + t221;
t42 = -pkin(3) * t102 - pkin(7) * t101 + t114;
t130 = -pkin(1) - pkin(6);
t111 = qJDD(1) * t130 + qJDD(2);
t113 = qJD(1) * t130 + qJD(2);
t229 = qJD(3) * t128;
t54 = t125 * t111 + t113 * t229;
t52 = qJDD(3) * pkin(7) + t54;
t174 = pkin(3) * t125 - pkin(7) * t128;
t104 = qJ(2) + t174;
t73 = t104 * qJD(1);
t103 = t125 * t113;
t80 = qJD(3) * pkin(7) + t103;
t3 = t124 * t42 + t127 * t52 + t73 * t226 - t227 * t80;
t2 = qJ(5) * t41 + qJD(5) * t94 + t3;
t32 = t124 * t73 + t127 * t80;
t4 = -qJD(4) * t32 - t124 * t52 + t127 * t42;
t332 = t4 * mrSges(5,1) - t3 * mrSges(5,2) - t2 * mrSges(6,2);
t331 = t128 / 0.2e1;
t330 = t321 * t282 + t338 * t283 + t339 * t284;
t327 = t338 * t95;
t126 = sin(qJ(1));
t129 = cos(qJ(1));
t173 = g(1) * t126 - g(2) * t129;
t326 = m(6) * pkin(4);
t325 = -m(6) - m(5);
t324 = t320 * t91 + t337 * t41 + t338 * t40;
t322 = -mrSges(5,1) - mrSges(6,1);
t300 = mrSges(5,2) + mrSges(6,2);
t12 = -mrSges(5,1) * t41 + mrSges(5,2) * t40;
t319 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t101 - t12;
t298 = t320 * t115 + t337 * t94 + t327;
t294 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t94 - mrSges(5,2) * t95 - mrSges(4,3) * t233;
t317 = -t313 * t125 + t336 * t128;
t316 = -t311 * t125 + t320 * t128;
t315 = -t309 * t125 + t321 * t128;
t314 = t321 * t124 + t320 * t127;
t312 = t337 * t127 + t341;
t310 = t339 * t124 + t342;
t172 = -t124 * t4 + t127 * t3;
t31 = -t124 * t80 + t127 * t73;
t308 = -t31 * t226 - t32 * t227 + t172;
t307 = t320 * t41 + t321 * t40 + t336 * t91;
t231 = qJD(3) * t125;
t53 = t111 * t128 - t113 * t231;
t306 = -t125 * t54 - t128 * t53;
t305 = -m(4) + t325;
t171 = mrSges(4,1) * t128 - mrSges(4,2) * t125;
t304 = (-Ifges(4,1) * t125 - t252) * t331 + qJ(2) * t171;
t303 = -mrSges(4,3) - mrSges(2,1) + mrSges(3,2);
t170 = mrSges(4,1) * t125 + mrSges(4,2) * t128;
t302 = -m(5) * t174 - t125 * t340 + t128 * t333 + mrSges(2,2) - mrSges(3,3) - t170;
t184 = qJD(4) * t123;
t207 = t124 * t222;
t224 = qJD(5) * t127;
t237 = t127 * t128;
t175 = pkin(3) * t128 + pkin(7) * t125;
t99 = t175 * qJD(1);
t48 = t113 * t237 + t124 * t99;
t299 = -qJ(5) * t207 + t124 * t184 + t224 - t48;
t240 = t125 * t127;
t242 = t124 * t128;
t47 = -t113 * t242 + t127 * t99;
t296 = -qJD(5) * t124 + t127 * t184 - (pkin(4) * t128 + qJ(5) * t240) * qJD(1) - t47;
t295 = mrSges(6,1) + t326;
t293 = -t103 + (t207 + t227) * pkin(4);
t167 = -mrSges(6,1) * t127 + mrSges(6,2) * t124;
t169 = -mrSges(5,1) * t127 + mrSges(5,2) * t124;
t292 = m(5) * pkin(3) - t167 - t169 + t340;
t291 = -m(5) * pkin(7) - t333;
t289 = -pkin(4) * t124 + t130;
t288 = -mrSges(5,1) - t295;
t1 = pkin(4) * t91 - qJ(5) * t40 - qJD(5) * t95 + t4;
t19 = -qJ(5) * t95 + t31;
t14 = pkin(4) * t115 + t19;
t20 = qJ(5) * t94 + t32;
t287 = -t1 * t124 + t127 * t2 - t14 * t226 - t20 * t227;
t131 = qJD(1) ^ 2;
t278 = t95 / 0.2e1;
t271 = mrSges(5,3) * t94;
t270 = mrSges(5,3) * t95;
t269 = mrSges(6,3) * t94;
t268 = mrSges(6,3) * t95;
t56 = -mrSges(6,2) * t115 + t269;
t57 = -mrSges(5,2) * t115 + t271;
t256 = t56 + t57;
t58 = mrSges(6,1) * t115 - t268;
t59 = mrSges(5,1) * t115 - t270;
t255 = -t58 - t59;
t254 = mrSges(6,2) * t127;
t239 = t125 * t130;
t61 = t124 * t104 + t127 * t239;
t245 = t113 * t128;
t244 = t124 * t125;
t243 = t124 * t126;
t241 = t124 * t129;
t238 = t126 * t127;
t236 = t127 * t129;
t234 = t129 * pkin(1) + t126 * qJ(2);
t228 = qJD(3) * t130;
t225 = qJD(4) * t128;
t223 = qJDD(1) * mrSges(3,2);
t204 = t128 * t228;
t92 = qJD(3) * t175 + qJD(2);
t216 = t104 * t226 + t124 * t92 + t127 * t204;
t209 = t124 * t239;
t206 = t124 * t231;
t205 = t125 * t228;
t203 = t127 * t225;
t202 = t125 * t230;
t11 = -t41 * mrSges(6,1) + t40 * mrSges(6,2);
t185 = -t124 * t130 + pkin(4);
t183 = -t220 / 0.2e1;
t182 = (t114 + t221) * qJ(2);
t168 = mrSges(5,1) * t124 + mrSges(5,2) * t127;
t166 = mrSges(6,1) * t124 + t254;
t155 = -Ifges(4,5) * t125 - Ifges(4,6) * t128;
t150 = t32 * t124 + t31 * t127;
t81 = -qJD(3) * pkin(3) - t245;
t51 = -qJDD(3) * pkin(3) - t53;
t76 = t125 * t241 + t238;
t74 = -t125 * t243 + t236;
t148 = t125 * (-Ifges(4,2) * t128 - t253);
t144 = -t131 * qJ(2) - t173;
t143 = t124 * t225 + t202;
t142 = -t203 + t206;
t141 = -t124 * t256 + t127 * t255;
t139 = -qJD(4) * t61 + t127 * t92;
t117 = -qJDD(1) * pkin(1) + qJDD(2);
t110 = t123 * t127;
t108 = t123 * t124;
t106 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t222;
t98 = t170 * qJD(1);
t93 = t289 * t128;
t90 = t127 * t104;
t85 = t168 * t128;
t77 = t125 * t236 - t243;
t75 = t125 * t238 + t241;
t66 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t102;
t60 = t90 - t209;
t55 = -pkin(4) * t142 + t205;
t49 = -mrSges(6,1) * t94 + mrSges(6,2) * t95;
t46 = -qJ(5) * t242 + t61;
t45 = -pkin(4) * t94 + qJD(5) + t81;
t43 = -qJ(5) * t237 + t125 * t185 + t90;
t22 = -t124 * t204 + t139;
t21 = -qJD(4) * t209 + t216;
t18 = -mrSges(5,2) * t91 + mrSges(5,3) * t41;
t17 = -mrSges(6,2) * t91 + mrSges(6,3) * t41;
t16 = mrSges(5,1) * t91 - mrSges(5,3) * t40;
t15 = mrSges(6,1) * t91 - mrSges(6,3) * t40;
t13 = -pkin(4) * t41 + qJDD(5) + t51;
t10 = -qJ(5) * t203 + (-qJD(5) * t128 + (qJ(5) * qJD(3) - qJD(4) * t130) * t125) * t124 + t216;
t9 = qJ(5) * t202 + (qJ(5) * t227 + qJD(3) * t185 - t224) * t128 + t139;
t5 = [-t335 * t231 / 0.2e1 + (Ifges(2,3) + Ifges(3,1)) * qJDD(1) + m(5) * (t130 * t231 * t81 + t21 * t32 + t22 * t31 + t3 * t61 + t4 * t60) + (t170 + 0.2e1 * mrSges(3,3)) * t114 + m(3) * (-pkin(1) * t117 + t182) + m(6) * (t1 * t43 + t10 * t20 - t13 * t93 + t14 * t9 + t2 * t46 + t45 * t55) + (Ifges(4,5) * qJDD(3) + t13 * t166 + t313 * t282 + t311 * t283 + t309 * t284 + (-m(5) * t51 + t319) * t130) * t128 + t117 * mrSges(3,2) + qJ(2) * (-mrSges(4,1) * t102 + mrSges(4,2) * t101) + qJD(2) * t98 - t93 * t11 + t51 * t85 + t55 * t49 + t10 * t56 + t21 * t57 + t9 * t58 + t22 * t59 + t60 * t16 + t61 * t18 + t43 * t15 + t46 * t17 - t294 * t205 + t304 * t220 + t306 * mrSges(4,3) + m(4) * (-t130 * t306 + t182) + (t31 * mrSges(5,1) + t14 * mrSges(6,1) - t32 * mrSges(5,2) - t20 * mrSges(6,2) - t334) * t229 + (t142 * t32 + t143 * t31 - t237 * t4 - t242 * t3) * mrSges(5,3) - (t124 * t297 + t127 * t298) * t225 / 0.2e1 + t102 * t343 + t298 * t206 / 0.2e1 + (-t1 * t237 + t14 * t143 + t142 * t20 - t2 * t242) * mrSges(6,3) + (qJD(3) * t315 - t225 * t310) * t278 + (qJD(3) * t316 - t225 * t312) * t94 / 0.2e1 + (qJD(3) * t317 - t225 * t314) * t115 / 0.2e1 + t106 * t204 - pkin(1) * t223 + t237 * t330 + (Ifges(4,1) * t101 + Ifges(4,4) * t102) * t331 - t324 * t242 / 0.2e1 + (-t241 * t326 - m(3) * t234 + t322 * t75 - t300 * t74 + t305 * (t129 * pkin(6) + t234) + t303 * t129 + t302 * t126) * g(2) + (t322 * t77 + t300 * t76 + (m(3) * pkin(1) - t289 * m(6) + (-m(4) - m(5)) * t130 - t303) * t126 + ((-m(3) + t305) * qJ(2) + t302) * t129) * g(1) + (-Ifges(4,6) * qJDD(3) - Ifges(4,4) * t101 / 0.2e1 + t102 * t347 + t307 / 0.2e1 + t1 * mrSges(6,1) + t336 * t282 + t320 * t283 + t321 * t284 + t332) * t125 + t148 * t183 + t66 * t239 + t45 * (-mrSges(6,1) * t142 - mrSges(6,2) * t143) + t81 * (-mrSges(5,1) * t142 - mrSges(5,2) * t143) + qJD(3) ^ 2 * t155 / 0.2e1 + t101 * t165 / 0.2e1; t223 - t131 * mrSges(3,3) + t144 * m(4) + (t117 + t144) * m(3) + (-t98 - m(5) * t150 - m(6) * (t20 * t124 + t14 * t127) + t141) * qJD(1) + (-t11 + (t124 * t255 + t127 * t256 + t106) * qJD(3) + m(5) * (t230 * t32 - t232 * t31 - t51) + m(6) * (-t14 * t232 + t20 * t230 - t13) + m(4) * t53 + t319) * t128 + (t66 + (t17 + t18) * t127 + (-t15 - t16) * t124 + (t49 - t294) * qJD(3) + t141 * qJD(4) + m(5) * (qJD(3) * t81 + t308) + m(6) * (qJD(3) * t45 + t287) + m(4) * t54) * t125 + t325 * t173; t335 * t222 / 0.2e1 - t106 * t245 + (t125 * t292 + t128 * t291 + t170) * g(3) + (t127 * t18 - t124 * t16 - t59 * t226 - t57 * t227 + m(5) * (-qJD(4) * t150 + t172)) * pkin(7) - t116 * t11 + Ifges(4,5) * t101 + Ifges(4,6) * t102 + t108 * t15 - t110 * t17 - t48 * t57 - t47 * t59 + t53 * mrSges(4,1) - t54 * mrSges(4,2) + (t148 / 0.2e1 - t304) * t131 + t115 * (t166 * t45 + t168 * t81) + t173 * (t125 * t291 - t128 * t292 - t171) + (-t207 / 0.2e1 - t227 / 0.2e1) * t298 + t334 * t233 + (-pkin(3) * t51 - t103 * t81 - t31 * t47 - t32 * t48) * m(5) - pkin(3) * t12 + t297 * t226 / 0.2e1 + t299 * t56 + (t1 * t108 - t110 * t2 - t116 * t13 + t14 * t296 + t20 * t299 + t293 * t45) * m(6) + t308 * mrSges(5,3) + t310 * t284 + t312 * t283 + (t313 * t115 + t309 * t95 + t311 * t94) * qJD(4) / 0.2e1 + t314 * t282 - (t317 * t115 + t315 * t95 + t316 * t94) * qJD(1) / 0.2e1 + t296 * t58 + (-t20 * (-mrSges(6,2) * t128 + mrSges(6,3) * t244) - t32 * (-mrSges(5,2) * t128 + mrSges(5,3) * t244) - t14 * (mrSges(6,1) * t128 + mrSges(6,3) * t240) - t31 * (mrSges(5,1) * t128 + mrSges(5,3) * t240)) * qJD(1) + t293 * t49 + t294 * t103 + t124 * t330 + t324 * t127 / 0.2e1 + Ifges(4,3) * qJDD(3) + t155 * t183 + t287 * mrSges(6,3) + t13 * t167 + t51 * t169; -t81 * (mrSges(5,1) * t95 + mrSges(5,2) * t94) - t45 * (mrSges(6,1) * t95 + mrSges(6,2) * t94) - t19 * t56 + t14 * t269 + (t59 + t270) * t32 + (-t57 + t271) * t31 + (t339 * t94 - t327) * t345 + t298 * t278 + (-t320 * t95 + t321 * t94) * t344 + (-m(6) * (-t14 + t19) + t58 + t268) * t20 + t295 * t1 + (t85 - (-t124 * t295 - t254) * t128) * g(3) + (t288 * t76 - t300 * t77) * g(2) + (t288 * t74 + t300 * t75) * g(1) + (-t337 * t95 + t297 + t329) * t346 + (t15 + (-m(6) * t45 - t49) * t95) * pkin(4) + t307 + t332; -t94 * t56 + t95 * t58 + (-g(3) * t125 + t128 * t173 + t14 * t95 - t20 * t94 + t13) * m(6) + t11;];
tau = t5;
