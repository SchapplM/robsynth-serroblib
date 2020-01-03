% Calculate vector of inverse dynamics joint torques for
% S5RPRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:10
% EndTime: 2019-12-31 18:42:30
% DurationCPUTime: 11.90s
% Computational Cost: add. (2777->445), mult. (5879->579), div. (0->0), fcn. (3402->10), ass. (0->201)
t345 = Ifges(5,4) + Ifges(6,4);
t346 = Ifges(5,1) + Ifges(6,1);
t320 = Ifges(5,5) + Ifges(6,5);
t344 = Ifges(5,2) + Ifges(6,2);
t343 = Ifges(5,6) + Ifges(6,6);
t352 = mrSges(5,3) + mrSges(6,3);
t138 = cos(qJ(4));
t351 = t345 * t138;
t135 = sin(qJ(4));
t350 = t345 * t135;
t231 = qJD(3) * t138;
t136 = sin(qJ(3));
t236 = qJD(1) * t136;
t101 = -t135 * t236 + t231;
t139 = cos(qJ(3));
t225 = qJD(1) * qJD(3);
t107 = qJDD(1) * t136 + t139 * t225;
t51 = qJD(4) * t101 + qJDD(3) * t135 + t107 * t138;
t289 = t51 / 0.2e1;
t233 = qJD(3) * t135;
t102 = t138 * t236 + t233;
t52 = -qJD(4) * t102 + qJDD(3) * t138 - t107 * t135;
t288 = t52 / 0.2e1;
t106 = qJDD(1) * t139 - t136 * t225;
t97 = qJDD(4) - t106;
t287 = t97 / 0.2e1;
t349 = -t101 / 0.2e1;
t348 = -t102 / 0.2e1;
t235 = qJD(1) * t139;
t329 = t235 - qJD(4);
t347 = t329 / 0.2e1;
t342 = Ifges(5,3) + Ifges(6,3);
t125 = pkin(4) * t138 + pkin(3);
t176 = -mrSges(6,1) * t138 + mrSges(6,2) * t135;
t178 = -mrSges(5,1) * t138 + mrSges(5,2) * t135;
t341 = m(5) * pkin(3) + m(6) * t125 - t176 - t178;
t305 = -t135 * t343 + t138 * t320;
t303 = -t135 * t344 + t351;
t301 = t138 * t346 - t350;
t134 = -qJ(5) - pkin(7);
t340 = -m(5) * pkin(7) + m(6) * t134 - t352;
t339 = t106 / 0.2e1;
t338 = t107 / 0.2e1;
t337 = t225 / 0.2e1;
t336 = t287 * t320 + t288 * t345 + t289 * t346;
t335 = t345 * t101;
t126 = Ifges(4,4) * t235;
t316 = t102 * t346 - t320 * t329 + t335;
t334 = Ifges(4,1) * t236 + Ifges(4,5) * qJD(3) + t138 * t316 + t126;
t132 = sin(pkin(8));
t121 = pkin(1) * t132 + pkin(6);
t333 = qJD(2) * qJD(3) + qJDD(1) * t121;
t227 = qJD(4) * t138;
t229 = qJD(4) * t135;
t110 = t121 * qJD(1);
t232 = qJD(3) * t136;
t33 = qJDD(2) * t136 - t110 * t232 + t139 * t333;
t28 = qJDD(3) * pkin(7) + t33;
t133 = cos(pkin(8));
t277 = pkin(1) * t133;
t122 = -pkin(2) - t277;
t109 = t122 * qJDD(1);
t50 = -pkin(3) * t106 - pkin(7) * t107 + t109;
t234 = qJD(2) * t136;
t75 = t110 * t139 + t234;
t67 = qJD(3) * pkin(7) + t75;
t183 = pkin(3) * t139 + pkin(7) * t136;
t159 = -pkin(2) - t183;
t92 = t159 - t277;
t68 = t92 * qJD(1);
t3 = t135 * t50 + t138 * t28 + t227 * t68 - t229 * t67;
t25 = t135 * t68 + t138 * t67;
t4 = -qJD(4) * t25 - t135 * t28 + t138 * t50;
t181 = -t135 * t4 + t138 * t3;
t24 = -t135 * t67 + t138 * t68;
t332 = -t227 * t24 - t229 * t25 + t181;
t331 = t345 * t102;
t131 = qJ(1) + pkin(8);
t127 = sin(t131);
t128 = cos(t131);
t330 = g(1) * t128 + g(2) * t127;
t256 = Ifges(4,4) * t136;
t170 = Ifges(4,2) * t139 + t256;
t328 = Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t170 / 0.2e1 + t342 * t347 + t320 * t348 + t343 * t349;
t2 = qJ(5) * t52 + qJD(5) * t101 + t3;
t327 = -t4 * mrSges(5,1) + t3 * mrSges(5,2) + t2 * mrSges(6,2);
t290 = m(6) * pkin(4);
t326 = t343 * t97 + t344 * t52 + t345 * t51;
t323 = -mrSges(6,1) - mrSges(5,1);
t322 = mrSges(3,2) - mrSges(4,3);
t321 = mrSges(5,2) + mrSges(6,2);
t317 = t101 * t344 - t329 * t343 + t331;
t13 = -mrSges(5,1) * t52 + mrSges(5,2) * t51;
t315 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t107 + t13;
t191 = qJD(4) * t134;
t209 = t135 * t235;
t226 = qJD(5) * t138;
t182 = pkin(3) * t136 - pkin(7) * t139;
t104 = t182 * qJD(1);
t74 = qJD(2) * t139 - t110 * t136;
t43 = t104 * t135 + t138 * t74;
t314 = qJ(5) * t209 + t135 * t191 + t226 - t43;
t237 = t138 * t139;
t158 = pkin(4) * t136 - qJ(5) * t237;
t42 = t104 * t138 - t135 * t74;
t313 = -qJD(1) * t158 - qJD(5) * t135 + t138 * t191 - t42;
t103 = t121 * t237;
t54 = t135 * t92 + t103;
t219 = mrSges(4,3) * t236;
t312 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t101 + mrSges(5,2) * t102 + t219;
t311 = t290 + mrSges(6,1);
t274 = pkin(4) * t135;
t310 = -t234 - (qJD(1) * t274 + t110) * t139 + pkin(4) * t229;
t309 = t136 * t342 + t139 * t305;
t308 = t136 * t343 + t139 * t303;
t307 = t136 * t320 + t139 * t301;
t304 = t135 * t320 + t138 * t343;
t302 = t138 * t344 + t350;
t300 = t135 * t346 + t351;
t298 = t320 * t51 + t342 * t97 + t343 * t52;
t230 = qJD(3) * t139;
t34 = qJDD(2) * t139 - t110 * t230 - t136 * t333;
t297 = -t136 * t34 + t139 * t33;
t296 = -m(6) - m(5) - m(4);
t295 = mrSges(5,1) + t311;
t180 = mrSges(4,1) * t139 - mrSges(4,2) * t136;
t292 = t136 * t352 + mrSges(3,1) + t180;
t1 = pkin(4) * t97 - qJ(5) * t51 - qJD(5) * t102 + t4;
t15 = -qJ(5) * t102 + t24;
t14 = -pkin(4) * t329 + t15;
t16 = qJ(5) * t101 + t25;
t291 = -t1 * t135 + t138 * t2 - t14 * t227 - t16 * t229;
t283 = t102 / 0.2e1;
t137 = sin(qJ(1));
t276 = pkin(1) * t137;
t140 = cos(qJ(1));
t130 = t140 * pkin(1);
t258 = mrSges(6,3) * t101;
t61 = mrSges(6,2) * t329 + t258;
t260 = mrSges(5,3) * t101;
t62 = mrSges(5,2) * t329 + t260;
t265 = t61 + t62;
t257 = mrSges(6,3) * t102;
t63 = -mrSges(6,1) * t329 - t257;
t259 = mrSges(5,3) * t102;
t64 = -mrSges(5,1) * t329 - t259;
t264 = -t63 - t64;
t105 = t182 * qJD(3);
t263 = t105 * t135 + t227 * t92;
t212 = t121 * t232;
t262 = t105 * t138 + t135 * t212;
t261 = mrSges(6,2) * t138;
t255 = Ifges(4,4) * t139;
t242 = t127 * t135;
t241 = t128 * t135;
t240 = t135 * t136;
t239 = t135 * t139;
t238 = t136 * t138;
t228 = qJD(4) * t136;
t218 = mrSges(4,3) * t235;
t211 = t121 * t230;
t210 = t135 * t230;
t12 = -mrSges(6,1) * t52 + mrSges(6,2) * t51;
t179 = mrSges(4,1) * t136 + mrSges(4,2) * t139;
t177 = mrSges(5,1) * t135 + mrSges(5,2) * t138;
t175 = mrSges(6,1) * t135 + t261;
t165 = Ifges(4,5) * t139 - Ifges(4,6) * t136;
t160 = t125 * t139 - t134 * t136;
t71 = t127 * t138 - t128 * t239;
t69 = t127 * t239 + t128 * t138;
t66 = -qJD(3) * pkin(3) - t74;
t155 = t122 * qJD(1) * t179;
t154 = t136 * (Ifges(4,1) * t139 - t256);
t151 = -t135 * t228 + t138 * t230;
t150 = t136 * t227 + t210;
t29 = -qJDD(3) * pkin(3) - t34;
t116 = t134 * t138;
t114 = t134 * t135;
t113 = -qJD(3) * mrSges(4,2) + t218;
t91 = t177 * t136;
t82 = (t121 + t274) * t136;
t78 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t106;
t77 = t138 * t92;
t72 = t128 * t237 + t242;
t70 = -t127 * t237 + t241;
t57 = pkin(4) * t150 + t211;
t55 = -mrSges(6,1) * t101 + mrSges(6,2) * t102;
t53 = -t121 * t239 + t77;
t49 = -pkin(4) * t101 + qJD(5) + t66;
t41 = -qJ(5) * t240 + t54;
t30 = -qJ(5) * t238 + t77 + (-t121 * t135 - pkin(4)) * t139;
t22 = -mrSges(5,2) * t97 + mrSges(5,3) * t52;
t21 = -mrSges(6,2) * t97 + mrSges(6,3) * t52;
t20 = mrSges(5,1) * t97 - mrSges(5,3) * t51;
t19 = mrSges(6,1) * t97 - mrSges(6,3) * t51;
t18 = -qJD(4) * t54 + t262;
t17 = (-t136 * t231 - t139 * t229) * t121 + t263;
t11 = -pkin(4) * t52 + qJDD(5) + t29;
t10 = (-qJ(5) * qJD(4) - qJD(3) * t121) * t238 + (-qJD(5) * t136 + (-qJ(5) * qJD(3) - qJD(4) * t121) * t139) * t135 + t263;
t5 = -t136 * t226 + t158 * qJD(3) + (-t103 + (qJ(5) * t136 - t92) * t135) * qJD(4) + t262;
t6 = [qJD(3) * t155 - (t135 * t316 + t138 * t317) * t228 / 0.2e1 + (m(4) * t122 - t180) * t109 + t154 * t337 + t255 * t338 + (t256 + t170) * t339 + m(5) * (t17 * t25 + t18 * t24 + t3 * t54 + t4 * t53) + (m(5) * (t136 * t29 + t230 * t66) + t139 * t78 + t315 * t136 + m(4) * ((-t136 * t75 - t139 * t74) * qJD(3) + t297)) * t121 - t113 * t212 + (-t150 * t25 - t151 * t24 - t238 * t4 - t240 * t3) * mrSges(5,3) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t133 - 0.2e1 * mrSges(3,2) * t132 + m(3) * (t132 ^ 2 + t133 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (Ifges(4,4) * t338 + Ifges(4,2) * t339 - t1 * mrSges(6,1) + Ifges(4,6) * qJDD(3) + (-Ifges(4,2) * t136 + t255) * t337 - t298 / 0.2e1 - t342 * t287 - t343 * t288 - t320 * t289 + t327) * t139 + (-t230 * t74 + t297) * mrSges(4,3) + t122 * (-mrSges(4,1) * t106 + mrSges(4,2) * t107) + m(6) * (t1 * t30 + t10 * t16 + t11 * t82 + t14 * t5 + t2 * t41 + t49 * t57) + t29 * t91 + t82 * t12 + t57 * t55 + t10 * t61 + t17 * t62 + t5 * t63 + t18 * t64 + t53 * t20 + t54 * t22 + t41 * t21 + t30 * t19 + t238 * t336 + (-t1 * t238 - t14 * t151 - t150 * t16 - t2 * t240) * mrSges(6,3) + t49 * (mrSges(6,1) * t150 + mrSges(6,2) * t151) + t66 * (mrSges(5,1) * t150 + mrSges(5,2) * t151) - (qJD(3) * t309 - t228 * t304) * t329 / 0.2e1 + qJD(3) ^ 2 * t165 / 0.2e1 + (Ifges(4,1) * t107 + Ifges(4,5) * qJDD(3) + t11 * t175 + t305 * t287 + t303 * t288 + t301 * t289) * t136 + (qJD(3) * t307 - t228 * t300) * t283 + (qJD(3) * t308 - t228 * t302) * t101 / 0.2e1 + (t24 * mrSges(5,1) + t14 * mrSges(6,1) - t25 * mrSges(5,2) - t16 * mrSges(6,2) - t75 * mrSges(4,3) - t328) * t232 + t312 * t211 + t334 * t230 / 0.2e1 - t317 * t210 / 0.2e1 - t326 * t240 / 0.2e1 + (-t242 * t290 - m(3) * t130 - mrSges(2,1) * t140 + mrSges(2,2) * t137 + t323 * t72 - t321 * t71 + t296 * (pkin(2) * t128 + pkin(6) * t127 + t130) + t322 * t127 + (-m(5) * t183 - m(6) * t160 - t292) * t128) * g(2) + (-t241 * t290 + m(3) * t276 + mrSges(2,1) * t137 + mrSges(2,2) * t140 + t323 * t70 - t321 * t69 + t296 * (pkin(6) * t128 - t276) + t322 * t128 + (-m(6) * (-pkin(2) - t160) - m(5) * t159 + m(4) * pkin(2) + t292) * t127) * g(1); m(3) * qJDD(2) + (-m(3) + t296) * g(3) + (-t12 + (t135 * t264 + t138 * t265 + t113) * qJD(3) + m(5) * (t231 * t25 - t233 * t24 - t29) + m(6) * (-t14 * t233 + t16 * t231 - t11) + m(4) * (qJD(3) * t75 + t34) - t315) * t139 + (t78 + (t21 + t22) * t138 + (-t19 - t20) * t135 + (t55 + t312) * qJD(3) + (-t135 * t265 + t138 * t264) * qJD(4) + m(5) * (qJD(3) * t66 + t332) + m(6) * (qJD(3) * t49 + t291) + m(4) * (-qJD(3) * t74 + t33)) * t136; (t138 * t22 - t64 * t227 - t135 * t20 - t62 * t229 + m(5) * ((-t135 * t25 - t138 * t24) * qJD(4) + t181)) * pkin(7) - t165 * t225 / 0.2e1 + (-t113 + t218) * t74 - t125 * t12 + t114 * t19 - t116 * t21 + Ifges(4,6) * t106 + Ifges(4,5) * t107 - t43 * t62 - t42 * t64 - t33 * mrSges(4,2) + t34 * mrSges(4,1) - pkin(3) * t13 + t135 * t336 + (-pkin(3) * t29 - t24 * t42 - t25 * t43) * m(5) + t300 * t289 + t302 * t288 + t304 * t287 + (t101 * t303 + t102 * t301 - t305 * t329) * qJD(4) / 0.2e1 - (t154 * qJD(1) + t308 * t101 + t102 * t307 - t309 * t329) * qJD(1) / 0.2e1 + (-t155 - t14 * (mrSges(6,1) * t136 - mrSges(6,3) * t237) - t24 * (mrSges(5,1) * t136 - mrSges(5,3) * t237) - t16 * (-mrSges(6,2) * t136 - mrSges(6,3) * t239) - t25 * (-mrSges(5,2) * t136 - mrSges(5,3) * t239)) * qJD(1) + (t209 / 0.2e1 - t229 / 0.2e1) * t317 + t291 * mrSges(6,3) + t11 * t176 + t328 * t236 + t29 * t178 + t329 * (-t49 * t175 - t66 * t177) + t332 * mrSges(5,3) + t310 * t55 + (-m(5) * t66 + t219 - t312) * t75 - (-Ifges(4,2) * t236 + t126 + t334) * t235 / 0.2e1 + t313 * t63 + t314 * t61 + (t1 * t114 - t11 * t125 - t116 * t2 + t14 * t313 + t16 * t314 + t310 * t49) * m(6) + t316 * t227 / 0.2e1 + t326 * t138 / 0.2e1 + t330 * (t136 * t341 + t139 * t340 + t179) + (t136 * t340 - t139 * t341 - t180) * g(3) + Ifges(4,3) * qJDD(3); (-t62 + t260) * t24 + t298 + t317 * t283 - t66 * (mrSges(5,1) * t102 + mrSges(5,2) * t101) - t49 * (mrSges(6,1) * t102 + mrSges(6,2) * t101) - t15 * t61 + (-m(6) * (-t14 + t15) + t63 + t257) * t16 + (t295 * t69 - t321 * t70) * g(2) + (-t295 * t71 + t321 * t72) * g(1) + (t101 * t320 - t102 * t343) * t347 + (t101 * t346 - t331) * t348 + (-t102 * t344 + t316 + t335) * t349 + (-(-t135 * t311 - t261) * t136 + t91) * g(3) + t14 * t258 + t311 * t1 + (t64 + t259) * t25 + (t19 + (-m(6) * t49 - t55) * t102) * pkin(4) - t327; -t101 * t61 + t102 * t63 + (g(3) * t139 - t16 * t101 + t14 * t102 - t136 * t330 + t11) * m(6) + t12;];
tau = t6;
