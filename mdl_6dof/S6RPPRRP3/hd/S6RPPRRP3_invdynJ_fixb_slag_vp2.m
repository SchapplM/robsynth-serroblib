% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:43
% EndTime: 2019-03-09 02:03:01
% DurationCPUTime: 11.52s
% Computational Cost: add. (3854->499), mult. (7093->641), div. (0->0), fcn. (3876->10), ass. (0->227)
t357 = Ifges(6,1) + Ifges(7,1);
t337 = Ifges(7,4) + Ifges(6,5);
t355 = Ifges(6,6) - Ifges(7,6);
t138 = sin(qJ(4));
t141 = cos(qJ(4));
t235 = qJD(1) * qJD(4);
t103 = qJDD(1) * t141 - t138 * t235;
t137 = sin(qJ(5));
t140 = cos(qJ(5));
t237 = t140 * qJD(4);
t246 = qJD(1) * t141;
t99 = t137 * t246 - t237;
t257 = qJD(5) * t99;
t51 = qJDD(4) * t137 + t103 * t140 - t257;
t305 = t51 / 0.2e1;
t244 = qJD(4) * t137;
t100 = t140 * t246 + t244;
t52 = qJD(5) * t100 - t140 * qJDD(4) + t103 * t137;
t303 = t52 / 0.2e1;
t104 = -qJDD(1) * t138 - t141 * t235;
t96 = qJDD(5) - t104;
t302 = t96 / 0.2e1;
t300 = t99 / 0.2e1;
t359 = -t100 / 0.2e1;
t247 = qJD(1) * t138;
t117 = qJD(5) + t247;
t358 = -t117 / 0.2e1;
t356 = Ifges(7,2) + Ifges(6,3);
t326 = -t355 * t137 + t337 * t140;
t266 = Ifges(7,5) * t137;
t268 = Ifges(6,4) * t137;
t324 = t357 * t140 + t266 - t268;
t304 = -t52 / 0.2e1;
t354 = t103 / 0.2e1;
t353 = t104 / 0.2e1;
t352 = t337 * t302 + (-Ifges(6,4) + Ifges(7,5)) * t303 + t357 * t305;
t136 = cos(pkin(9));
t123 = -pkin(1) * t136 - pkin(2);
t116 = -pkin(7) + t123;
t351 = -qJD(2) * qJD(4) + qJDD(1) * t116 + qJDD(3);
t188 = t140 * mrSges(7,1) + t137 * mrSges(7,3);
t190 = mrSges(6,1) * t140 - mrSges(6,2) * t137;
t350 = -t188 - t190;
t240 = qJD(5) * t140;
t241 = qJD(5) * t137;
t242 = qJD(4) * t141;
t94 = qJD(1) * t116 + qJD(3);
t29 = t141 * qJDD(2) + t138 * t351 + t94 * t242;
t27 = qJDD(4) * pkin(8) + t29;
t135 = sin(pkin(9));
t120 = pkin(1) * t135 + qJ(3);
t236 = qJD(1) * qJD(3);
t95 = qJDD(1) * t120 + t236;
t36 = -pkin(4) * t104 - pkin(8) * t103 + t95;
t245 = qJD(2) * t141;
t68 = t138 * t94 + t245;
t61 = qJD(4) * pkin(8) + t68;
t132 = t141 * pkin(8);
t287 = pkin(4) * t138;
t88 = -t132 + t120 + t287;
t69 = t88 * qJD(1);
t3 = t137 * t36 + t140 * t27 + t69 * t240 - t241 * t61;
t21 = t137 * t69 + t140 * t61;
t4 = -qJD(5) * t21 - t137 * t27 + t140 * t36;
t193 = -t137 * t4 + t140 * t3;
t20 = -t137 * t61 + t140 * t69;
t349 = -t20 * t240 - t21 * t241 + t193;
t16 = -pkin(5) * t117 + qJD(6) - t20;
t17 = qJ(6) * t117 + t21;
t1 = qJ(6) * t96 + qJD(6) * t117 + t3;
t2 = -pkin(5) * t96 + qJDD(6) - t4;
t194 = t1 * t140 + t137 * t2;
t348 = t16 * t240 - t17 * t241 + t194;
t23 = mrSges(6,1) * t96 - mrSges(6,3) * t51;
t24 = -t96 * mrSges(7,1) + t51 * mrSges(7,2);
t25 = -mrSges(6,2) * t96 - mrSges(6,3) * t52;
t26 = -mrSges(7,2) * t52 + mrSges(7,3) * t96;
t347 = (t25 + t26) * t140 + (-t23 + t24) * t137;
t271 = Ifges(5,4) * t138;
t186 = t141 * Ifges(5,1) - t271;
t289 = Ifges(7,5) * t99;
t92 = Ifges(6,4) * t99;
t334 = t100 * t357 + t337 * t117 + t289 - t92;
t91 = Ifges(7,5) * t100;
t38 = t117 * Ifges(7,6) + t99 * Ifges(7,3) + t91;
t346 = Ifges(5,5) * qJD(4) + qJD(1) * t186 + t137 * t38 + t140 * t334;
t270 = Ifges(5,4) * t141;
t181 = -t138 * Ifges(5,2) + t270;
t345 = t355 * t300 + t356 * t358 + t337 * t359 + Ifges(5,6) * qJD(4) / 0.2e1 + qJD(1) * t181 / 0.2e1;
t226 = mrSges(5,3) * t247;
t108 = -qJD(4) * mrSges(5,2) - t226;
t11 = mrSges(7,1) * t52 - mrSges(7,3) * t51;
t272 = mrSges(6,3) * t100;
t64 = mrSges(6,1) * t117 - t272;
t273 = mrSges(7,2) * t100;
t65 = -mrSges(7,1) * t117 + t273;
t275 = t64 - t65;
t290 = mrSges(6,3) * t99;
t63 = -mrSges(6,2) * t117 - t290;
t291 = mrSges(7,2) * t99;
t66 = mrSges(7,3) * t117 - t291;
t276 = t63 + t66;
t243 = qJD(4) * t138;
t30 = -t138 * qJDD(2) + t141 * t351 - t94 * t243;
t28 = -qJDD(4) * pkin(4) - t30;
t12 = mrSges(6,1) * t52 + mrSges(6,2) * t51;
t336 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t103 - t12;
t5 = pkin(5) * t52 - qJ(6) * t51 - qJD(6) * t100 + t28;
t344 = (t137 * t275 - t140 * t276 - t108) * qJD(4) - m(5) * (qJD(4) * t68 + t30) + m(6) * (t20 * t244 - t21 * t237 + t28) - m(7) * (t16 * t244 + t17 * t237 - t5) + t11 - t336;
t343 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t299 = -m(5) - m(4);
t342 = m(6) + m(7);
t338 = mrSges(6,3) + mrSges(7,2);
t170 = pkin(5) * t137 - qJ(6) * t140;
t333 = qJD(5) * t170 - qJD(6) * t137 - t245 - (-qJD(1) * t170 + t94) * t138;
t225 = mrSges(5,3) * t246;
t332 = -qJD(4) * mrSges(5,1) + mrSges(6,1) * t99 + mrSges(6,2) * t100 + t225;
t134 = qJ(1) + pkin(9);
t125 = sin(t134);
t255 = t125 * t138;
t254 = t125 * t141;
t171 = pkin(5) * t140 + qJ(6) * t137;
t106 = -pkin(4) - t171;
t330 = m(6) * pkin(4) - m(7) * t106 - t350;
t329 = -t138 * t326 + t141 * t356;
t328 = -t138 * t324 + t141 * t337;
t327 = t137 * t337 + t140 * t355;
t265 = Ifges(7,5) * t140;
t267 = Ifges(6,4) * t140;
t325 = t137 * t357 - t265 + t267;
t321 = t337 * t51 - t355 * t52 + t356 * t96;
t320 = t138 * t29 + t141 * t30;
t319 = -t51 * Ifges(7,5) / 0.2e1 - t96 * Ifges(7,6) / 0.2e1 + Ifges(6,4) * t305 + Ifges(6,6) * t302 + (Ifges(7,3) + Ifges(6,2)) * t304;
t318 = mrSges(4,2) - mrSges(3,1) - mrSges(5,3);
t317 = g(1) * t255;
t191 = mrSges(5,1) * t138 + mrSges(5,2) * t141;
t316 = -mrSges(4,3) + mrSges(3,2) - t191;
t250 = t138 * t140;
t274 = t116 * t250 + t137 * t88;
t195 = pkin(4) * t141 + pkin(8) * t138;
t97 = qJD(4) * t195 + qJD(3);
t312 = -qJD(5) * t274 + t140 * t97;
t311 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t310 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t301 = -t99 / 0.2e1;
t298 = -pkin(2) - pkin(7);
t296 = t100 / 0.2e1;
t139 = sin(qJ(1));
t288 = pkin(1) * t139;
t284 = g(1) * t125;
t126 = cos(t134);
t283 = g(2) * t126;
t142 = cos(qJ(1));
t133 = t142 * pkin(1);
t278 = -qJD(1) / 0.2e1;
t277 = qJD(5) / 0.2e1;
t102 = t195 * qJD(1);
t67 = -t138 * qJD(2) + t141 * t94;
t35 = t137 * t102 + t140 * t67;
t269 = Ifges(6,4) * t100;
t262 = t140 * t88;
t253 = t126 * t141;
t252 = t137 * t138;
t251 = t137 * t141;
t248 = pkin(4) * t254 + pkin(8) * t255;
t239 = qJD(5) * t141;
t238 = qJDD(1) * mrSges(4,2);
t217 = t116 * t242;
t231 = t137 * t97 + t140 * t217 + t88 * t240;
t227 = -t299 + t342;
t219 = t126 * pkin(2) + t125 * qJ(3) + t133;
t216 = t116 * t241;
t215 = t137 * t243;
t205 = t240 / 0.2e1;
t204 = -t239 / 0.2e1;
t203 = t126 * qJ(3) - t288;
t202 = t116 * t137 - pkin(5);
t201 = -t235 / 0.2e1;
t199 = t126 * pkin(7) + t219;
t192 = mrSges(5,1) * t141 - mrSges(5,2) * t138;
t189 = mrSges(6,1) * t137 + mrSges(6,2) * t140;
t187 = t137 * mrSges(7,1) - t140 * mrSges(7,3);
t180 = -Ifges(6,2) * t137 + t267;
t179 = Ifges(6,2) * t140 + t268;
t176 = -Ifges(5,5) * t138 - Ifges(5,6) * t141;
t173 = Ifges(7,3) * t137 + t265;
t172 = -Ifges(7,3) * t140 + t266;
t169 = t137 * t17 - t140 * t16;
t168 = t137 * t21 + t140 * t20;
t34 = t102 * t140 - t137 * t67;
t107 = t120 * qJD(1);
t166 = qJD(3) * t107 + t120 * t95;
t60 = -qJD(4) * pkin(4) - t67;
t164 = -t116 + t170;
t162 = t107 * t192;
t161 = t138 * (-Ifges(5,2) * t141 - t271);
t160 = t141 * (-Ifges(5,1) * t138 - t270);
t158 = t137 * t239 + t138 * t237;
t157 = -t140 * t239 + t215;
t156 = -t137 * t276 - t140 * t275;
t149 = Ifges(6,6) * t141 - t138 * t180;
t148 = Ifges(7,6) * t141 - t138 * t173;
t19 = pkin(5) * t99 - qJ(6) * t100 + t60;
t56 = mrSges(7,1) * t99 - mrSges(7,3) * t100;
t80 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t104;
t144 = t80 + (t56 + t332) * qJD(4) + t156 * qJD(5) + m(6) * (qJD(4) * t60 + t349) + m(7) * (qJD(4) * t19 + t348) + m(5) * (-qJD(4) * t67 + t29) + t347;
t105 = qJDD(1) * t123 + qJDD(3);
t101 = t191 * qJD(1);
t90 = t189 * t141;
t73 = -t125 * t137 + t126 * t250;
t72 = t125 * t140 + t126 * t252;
t71 = t125 * t250 + t126 * t137;
t70 = t125 * t252 - t126 * t140;
t62 = t164 * t141;
t55 = pkin(5) * t100 + qJ(6) * t99;
t53 = -t116 * t252 + t262;
t50 = t138 * t202 - t262;
t49 = qJ(6) * t138 + t274;
t41 = -t99 * Ifges(6,2) + t117 * Ifges(6,6) + t269;
t32 = -pkin(5) * t246 - t34;
t31 = qJ(6) * t246 + t35;
t22 = (qJD(5) * t171 - qJD(6) * t140) * t141 - t164 * t243;
t15 = -t137 * t217 + t312;
t14 = -t138 * t216 + t231;
t13 = t202 * t242 - t312;
t10 = qJ(6) * t242 + (qJD(6) - t216) * t138 + t231;
t6 = [t334 * t137 * t204 + (-m(3) * t133 - m(4) * t219 - m(5) * t199 - mrSges(2,1) * t142 + mrSges(2,2) * t139 - t342 * (pkin(4) * t255 - pkin(8) * t254 + t199) - t311 * t71 - t310 * t70 + t338 * t254 + t318 * t126 + t316 * t125) * g(2) + (m(3) * t288 + mrSges(2,1) * t139 + mrSges(2,2) * t142 - t342 * (-pkin(8) * t253 + t125 * t298 + t126 * t287 + t203) - t311 * t73 - t310 * t72 + t338 * t253 + t299 * t203 + t316 * t126 + (m(4) * pkin(2) - m(5) * t298 - t318) * t125) * g(1) + (qJD(4) * t328 - t239 * t325) * t296 + (qJD(4) * t329 - t239 * t327) * t117 / 0.2e1 + t120 * (-mrSges(5,1) * t104 + mrSges(5,2) * t103) + (Ifges(5,1) * t354 + Ifges(5,4) * t353 + Ifges(5,5) * qJDD(4) + t173 * t303 + t180 * t304 + t5 * t187 + t38 * t205 + t326 * t302 + t324 * t305) * t141 + t105 * mrSges(4,2) + qJD(3) * t101 + t28 * t90 - t319 * t251 + (t337 * t305 + t321 / 0.2e1 + t356 * t302 - Ifges(5,6) * qJDD(4) - Ifges(5,4) * t103 / 0.2e1 - Ifges(5,2) * t104 / 0.2e1 + Ifges(7,6) * t303 + Ifges(6,6) * t304 + t343) * t138 + (t95 + t236) * mrSges(4,3) + (t2 * mrSges(7,2) - t4 * mrSges(6,3) + t352) * t140 * t141 + m(5) * t166 + m(4) * (t105 * t123 + t166) + t181 * t353 + t186 * t354 + t95 * t191 + t62 * t11 + t14 * t63 + t15 * t64 + t13 * t65 + t10 * t66 - t346 * t243 / 0.2e1 + t53 * t23 + t22 * t56 + (t20 * mrSges(6,1) - t16 * mrSges(7,1) - t21 * mrSges(6,2) - t68 * mrSges(5,3) + t17 * mrSges(7,3) - t345) * t242 + t49 * t26 + t50 * t24 + t161 * t201 + qJD(4) ^ 2 * t176 / 0.2e1 + t160 * t235 / 0.2e1 + t274 * t25 + m(7) * (t1 * t49 + t10 * t17 + t13 * t16 + t19 * t22 + t2 * t50 + t5 * t62) + (t243 * t67 - t320) * mrSges(5,3) + (t157 * t21 + t158 * t20 - t251 * t3) * mrSges(6,3) + (t120 * mrSges(4,3) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + (0.2e1 * mrSges(3,1) * t136 - 0.2e1 * mrSges(3,2) * t135 + m(3) * (t135 ^ 2 + t136 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t19 * (-mrSges(7,1) * t157 + mrSges(7,3) * t158) + t60 * (-mrSges(6,1) * t157 - mrSges(6,2) * t158) + (qJD(4) * t149 - t179 * t239) * t301 + (qJD(4) * t148 - t172 * t239) * t300 + (-t1 * t251 + t157 * t17 - t158 * t16) * mrSges(7,2) + (t215 / 0.2e1 + t140 * t204) * t41 + m(6) * (t14 * t21 + t15 * t20 + t3 * t274 + t4 * t53) + (t336 * t141 + t332 * t243 + m(5) * ((-t138 * t67 + t141 * t68) * qJD(4) + t320) + t138 * t80 + m(6) * (-t141 * t28 + t243 * t60)) * t116 + t108 * t217 + qJD(4) * t162 + t123 * t238; (m(4) + m(3)) * qJDD(2) + (-m(3) - t227) * g(3) + t344 * t138 + t144 * t141; m(4) * t105 + t238 + (-m(6) * t168 - m(7) * t169 - mrSges(4,3) * qJD(1) + t107 * t299 - t101 + t156) * qJD(1) - t344 * t141 + t144 * t138 + (t283 - t284) * t227; -t192 * t284 + t333 * t56 + t334 * t205 + (-m(6) * t60 + t225 - t332) * t68 + t327 * t302 + t325 * t305 + t319 * t140 + Ifges(5,6) * t104 + t106 * t11 + Ifges(5,5) * t103 + (-t180 / 0.2e1 + t173 / 0.2e1) * t257 + t137 * t352 + (-t132 * t342 + t138 * t330 - t141 * t338 + t191) * g(3) + (-t317 + t348) * mrSges(7,2) + (-t317 + t349) * mrSges(6,3) + (-m(7) * t248 + (-m(7) * t171 + t350) * t254) * g(1) - t35 * t63 - t34 * t64 - t32 * t65 - t31 * t66 + t346 * t247 / 0.2e1 + (-t275 * t240 - t276 * t241 + t342 * t138 * t283 + (-qJD(5) * t169 + t194) * m(7) + (-qJD(5) * t168 + t193) * m(6) + t347) * pkin(8) + t345 * t246 - t29 * mrSges(5,2) + t30 * mrSges(5,1) - pkin(4) * t12 + t176 * t201 - t5 * t188 - t28 * t190 + (-t108 - t226) * t67 + (t277 * t324 + t278 * t328) * t100 + ((t149 / 0.2e1 - t148 / 0.2e1) * t99 - t162 - t16 * (-mrSges(7,1) * t141 - mrSges(7,2) * t250) - t20 * (mrSges(6,1) * t141 + mrSges(6,3) * t250) - t21 * (-mrSges(6,2) * t141 + mrSges(6,3) * t252) - t17 * (mrSges(7,2) * t252 + mrSges(7,3) * t141) + (t161 / 0.2e1 - t160 / 0.2e1) * qJD(1)) * qJD(1) + t38 * t241 / 0.2e1 + (t138 * t338 + t330 * t141 + t192) * t283 + (t106 * t5 - t16 * t32 - t17 * t31 + t333 * t19) * m(7) + (-pkin(4) * t28 - g(1) * t248 - t20 * t34 - t21 * t35) * m(6) + t172 * t303 + t179 * t304 + (t187 * t19 + t189 * t60 - t137 * t41 / 0.2e1 + t277 * t326 + t278 * t329) * t117 + Ifges(5,3) * qJDD(4); (-m(7) * t16 + t272 + t275) * t21 + (-m(7) * t17 - t276 - t290) * t20 + (-Ifges(6,2) * t100 + t334 - t92) * t300 + (-t100 * t355 - t337 * t99) * t358 + (-t310 * t71 + t311 * t70) * g(1) + (t310 * t73 - t311 * t72) * g(2) - t19 * (mrSges(7,1) * t100 + mrSges(7,3) * t99) - t60 * (mrSges(6,1) * t100 - mrSges(6,2) * t99) + qJD(6) * t66 - t55 * t56 - pkin(5) * t24 + qJ(6) * t26 + (-t357 * t99 - t269 + t38 + t91) * t359 + (t90 - (-m(7) * t170 - t187) * t141) * g(3) + t343 + (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t17 - t19 * t55) * m(7) + (Ifges(7,3) * t100 - t289) * t301 + t16 * t291 + t41 * t296 + t321 + t17 * t273; t100 * t56 - t117 * t66 + (-g(1) * t70 + g(2) * t72 - g(3) * t251 + t19 * t100 - t17 * t117 + t2) * m(7) + t24;];
tau  = t6;
