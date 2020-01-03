% Calculate vector of inverse dynamics joint torques for
% S5RRPPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR11_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR11_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:27
% EndTime: 2019-12-31 19:46:54
% DurationCPUTime: 15.34s
% Computational Cost: add. (3624->557), mult. (7727->751), div. (0->0), fcn. (4729->10), ass. (0->252)
t373 = qJD(2) / 0.2e1;
t188 = cos(pkin(8));
t194 = cos(qJ(5));
t187 = sin(pkin(8));
t191 = sin(qJ(5));
t280 = t187 * t191;
t348 = t188 * t194 - t280;
t109 = t348 * qJD(5);
t192 = sin(qJ(2));
t268 = qJD(1) * t192;
t251 = t188 * t268;
t85 = t194 * t251 - t268 * t280;
t355 = t109 + t85;
t217 = t194 * t187 + t191 * t188;
t110 = t217 * qJD(5);
t207 = t217 * t192;
t86 = qJD(1) * t207;
t372 = t110 + t86;
t195 = cos(qJ(2));
t261 = qJD(1) * qJD(2);
t138 = qJDD(1) * t192 + t195 * t261;
t137 = -t195 * qJDD(1) + t192 * t261;
t264 = qJD(3) * t192;
t215 = -qJD(4) * t195 - t264;
t285 = qJDD(1) * pkin(1);
t220 = -qJ(3) * t138 - t285;
t302 = pkin(2) + qJ(4);
t36 = qJD(1) * t215 + t137 * t302 + t220;
t123 = t138 * pkin(6);
t244 = qJDD(3) + t123;
t55 = pkin(3) * t138 - qJD(2) * qJD(4) - qJDD(2) * t302 + t244;
t16 = -t187 * t36 + t188 * t55;
t88 = qJDD(2) * t188 + t137 * t187;
t10 = pkin(4) * t138 - pkin(7) * t88 + t16;
t17 = t187 * t55 + t188 * t36;
t87 = -qJDD(2) * t187 + t137 * t188;
t11 = pkin(7) * t87 + t17;
t263 = t188 * qJD(2);
t267 = qJD(1) * t195;
t124 = t187 * t267 - t263;
t175 = t192 * qJ(3);
t243 = -pkin(1) - t175;
t89 = (-t195 * t302 + t243) * qJD(1);
t166 = pkin(6) * t268;
t133 = -pkin(3) * t268 - t166;
t97 = -qJD(2) * t302 + qJD(3) - t133;
t41 = -t187 * t89 + t188 * t97;
t29 = pkin(4) * t268 + pkin(7) * t124 + t41;
t125 = -qJD(2) * t187 - t188 * t267;
t42 = t187 * t97 + t188 * t89;
t30 = pkin(7) * t125 + t42;
t8 = -t191 * t30 + t194 * t29;
t1 = qJD(5) * t8 + t10 * t191 + t11 * t194;
t371 = t1 * mrSges(6,2);
t9 = t191 * t29 + t194 * t30;
t2 = -qJD(5) * t9 + t10 * t194 - t11 * t191;
t370 = t2 * mrSges(6,1);
t362 = m(5) + m(6);
t262 = m(4) + t362;
t360 = Ifges(3,1) + Ifges(5,3);
t359 = -Ifges(4,4) + Ifges(3,5);
t358 = Ifges(4,5) - Ifges(3,6);
t182 = pkin(8) + qJ(5);
t172 = sin(t182);
t173 = cos(t182);
t231 = mrSges(5,1) * t187 + mrSges(5,2) * t188;
t311 = pkin(4) * t187;
t369 = -m(6) * t311 - t172 * mrSges(6,1) - t173 * mrSges(6,2) - t231;
t189 = -pkin(7) - qJ(4);
t368 = m(5) * t302 + mrSges(5,3) - m(6) * (-pkin(2) + t189) + mrSges(6,3);
t295 = Ifges(4,6) * t195;
t223 = -t192 * Ifges(4,2) - t295;
t367 = t9 * mrSges(6,2) + Ifges(4,4) * t373 + qJD(1) * t223 / 0.2e1 - t8 * mrSges(6,1);
t193 = sin(qJ(1));
t366 = g(2) * t193;
t230 = t195 * mrSges(4,2) - t192 * mrSges(4,3);
t234 = mrSges(3,1) * t195 - mrSges(3,2) * t192;
t365 = t230 - t234;
t162 = pkin(4) * t188 + pkin(3);
t232 = mrSges(5,1) * t188 - mrSges(5,2) * t187;
t364 = -m(5) * pkin(3) - m(6) * t162 - mrSges(4,1) + mrSges(2,2) - mrSges(3,3) - t232;
t363 = -t1 * t217 - t2 * t348 - t355 * t9 + t372 * t8;
t240 = t124 * t191 + t194 * t125;
t23 = qJD(5) * t240 + t191 * t87 + t194 * t88;
t330 = t23 / 0.2e1;
t66 = t124 * t194 - t125 * t191;
t24 = qJD(5) * t66 - t191 * t88 + t194 * t87;
t329 = t24 / 0.2e1;
t126 = qJDD(5) + t138;
t317 = t126 / 0.2e1;
t361 = -t138 / 0.2e1;
t279 = t187 * t192;
t214 = pkin(4) * t195 - pkin(7) * t279;
t167 = pkin(2) * t268;
t286 = qJ(3) * t195;
t219 = qJ(4) * t192 - t286;
t102 = qJD(1) * t219 + t167;
t168 = pkin(6) * t267;
t134 = pkin(3) * t267 + t168;
t56 = -t102 * t187 + t188 * t134;
t38 = qJD(1) * t214 + t56;
t57 = t188 * t102 + t187 * t134;
t46 = pkin(7) * t251 + t57;
t301 = -pkin(7) - t302;
t139 = t301 * t187;
t140 = t301 * t188;
t75 = -t139 * t191 + t140 * t194;
t357 = -qJD(4) * t217 + qJD(5) * t75 - t191 * t38 - t194 * t46;
t76 = t139 * t194 + t140 * t191;
t356 = -qJD(4) * t348 - qJD(5) * t76 + t191 * t46 - t194 * t38;
t196 = cos(qJ(1));
t344 = g(1) * t196 + t366;
t353 = t195 * t344;
t253 = mrSges(4,1) * t267;
t147 = -qJD(2) * mrSges(4,3) - t253;
t352 = qJD(2) * mrSges(3,2) - mrSges(3,3) * t267 + t147;
t252 = mrSges(4,1) * t268;
t351 = mrSges(3,3) * t268 + t252 + (-mrSges(3,1) + mrSges(4,2)) * qJD(2);
t296 = Ifges(4,6) * t192;
t350 = t192 * (-Ifges(4,2) * t195 + t296) + t195 * (Ifges(4,3) * t192 - t295);
t349 = t192 * t358 + t195 * t359;
t122 = t137 * pkin(6);
t347 = -t122 * t195 + t123 * t192;
t101 = -qJDD(2) * pkin(2) + t244;
t90 = -qJDD(2) * qJ(3) - qJD(2) * qJD(3) + t122;
t346 = t101 * t192 - t195 * t90;
t58 = -mrSges(5,2) * t138 + mrSges(5,3) * t87;
t59 = mrSges(5,1) * t138 - mrSges(5,3) * t88;
t345 = t187 * t58 + t188 * t59;
t157 = qJD(5) + t268;
t165 = Ifges(3,4) * t267;
t342 = Ifges(3,5) * qJD(2) - t124 * Ifges(5,5) - Ifges(6,5) * t66 + t125 * Ifges(5,6) + Ifges(6,6) * t240 + t157 * Ifges(6,3) + t268 * t360 + t165;
t28 = -mrSges(6,1) * t240 - mrSges(6,2) * t66;
t74 = -mrSges(5,1) * t125 - mrSges(5,2) * t124;
t341 = t147 - t74 - t28;
t340 = -mrSges(2,1) + t365;
t290 = t195 * mrSges(4,3);
t339 = -t290 + t369 * t195 + (m(4) * pkin(2) - mrSges(4,2) + t368) * t192;
t222 = -t195 * Ifges(4,3) - t296;
t336 = t187 * (-t124 * Ifges(5,1) + t125 * Ifges(5,4) + Ifges(5,5) * t268) + t188 * (-t124 * Ifges(5,4) + t125 * Ifges(5,2) + Ifges(5,6) * t268) + Ifges(4,5) * qJD(2) + qJD(1) * t222;
t186 = qJD(2) * qJ(3);
t107 = qJD(4) + t186 + t134;
t179 = t195 * pkin(2);
t213 = t243 - t179;
t119 = t213 * qJD(1);
t141 = -qJD(2) * pkin(2) + qJD(3) + t166;
t145 = -t168 - t186;
t298 = Ifges(5,4) * t187;
t226 = Ifges(5,2) * t188 + t298;
t297 = Ifges(5,4) * t188;
t229 = Ifges(5,1) * t187 + t297;
t292 = t188 * mrSges(5,3);
t334 = t125 * (Ifges(5,6) * t195 + t192 * t226) / 0.2e1 - t124 * (Ifges(5,5) * t195 + t192 * t229) / 0.2e1 - t107 * t192 * t232 + t42 * (-mrSges(5,2) * t195 + t192 * t292) + t41 * (t195 * mrSges(5,1) - mrSges(5,3) * t279) + t119 * (-mrSges(4,2) * t192 - t290) + m(4) * (t141 * t195 + t145 * t192) * pkin(6);
t332 = Ifges(6,4) * t330 + Ifges(6,2) * t329 + Ifges(6,6) * t317;
t331 = Ifges(6,1) * t330 + Ifges(6,4) * t329 + Ifges(6,5) * t317;
t64 = Ifges(6,4) * t240;
t27 = -Ifges(6,1) * t66 + t157 * Ifges(6,5) + t64;
t327 = t27 / 0.2e1;
t326 = -t88 * Ifges(5,4) / 0.2e1 - t87 * Ifges(5,2) / 0.2e1 + Ifges(5,6) * t361;
t325 = -t240 / 0.2e1;
t324 = t240 / 0.2e1;
t323 = t66 / 0.2e1;
t322 = -t66 / 0.2e1;
t321 = t87 / 0.2e1;
t320 = t88 / 0.2e1;
t318 = pkin(3) + pkin(6);
t316 = t138 / 0.2e1;
t315 = -t157 / 0.2e1;
t314 = t157 / 0.2e1;
t312 = Ifges(6,4) * t66;
t310 = pkin(6) * t192;
t307 = g(3) * t195;
t177 = t195 * pkin(6);
t300 = Ifges(3,4) * t192;
t299 = Ifges(3,4) * t195;
t294 = t17 * t187;
t289 = t195 * mrSges(6,3);
t265 = qJD(2) * t195;
t136 = t318 * t265;
t266 = qJD(2) * t192;
t170 = pkin(2) * t266;
t77 = qJD(2) * t219 + t170 + t215;
t44 = t187 * t136 + t188 * t77;
t278 = t187 * t195;
t276 = t188 * t195;
t275 = t189 * t195;
t274 = t192 * t196;
t273 = t193 * t172;
t272 = t193 * t173;
t271 = t195 * t196;
t270 = t179 + t175;
t235 = qJ(4) * t195 + t270;
t121 = -pkin(1) - t235;
t150 = t318 * t192;
t71 = t188 * t121 + t187 * t150;
t151 = t195 * pkin(3) + t177;
t269 = t196 * pkin(1) + t193 * pkin(6);
t256 = Ifges(6,5) * t23 + Ifges(6,6) * t24 + Ifges(6,3) * t126;
t45 = -t87 * mrSges(5,1) + t88 * mrSges(5,2);
t7 = -t24 * mrSges(6,1) + t23 * mrSges(6,2);
t159 = qJ(3) + t311;
t242 = -t261 / 0.2e1;
t104 = t138 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t43 = t188 * t136 - t187 * t77;
t233 = mrSges(3,1) * t192 + mrSges(3,2) * t195;
t228 = t195 * Ifges(3,2) + t300;
t224 = Ifges(5,5) * t187 + Ifges(5,6) * t188;
t221 = t16 * t188 + t294;
t130 = t188 * t150;
t49 = pkin(4) * t192 + t130 + (pkin(7) * t195 - t121) * t187;
t54 = -pkin(7) * t276 + t71;
t18 = -t191 * t54 + t194 * t49;
t19 = t191 * t49 + t194 * t54;
t212 = pkin(1) * t233;
t210 = t192 * (Ifges(3,1) * t195 - t300);
t98 = t348 * t195;
t199 = t192 * (Ifges(5,3) * t195 + t192 * t224);
t65 = -pkin(3) * t137 + qJDD(4) - t90;
t142 = -pkin(1) - t270;
t135 = t318 * t266;
t132 = -qJ(3) * t267 + t167;
t131 = t230 * qJD(1);
t115 = Ifges(3,6) * qJD(2) + qJD(1) * t228;
t108 = pkin(4) * t276 + t151;
t105 = -qJ(3) * t265 + t170 - t264;
t103 = mrSges(4,1) * t137 - qJDD(2) * mrSges(4,3);
t99 = t217 * t195;
t96 = t173 * t196 - t192 * t273;
t95 = t172 * t196 + t192 * t272;
t94 = t172 * t274 + t272;
t93 = t173 * t274 - t273;
t92 = (-pkin(6) - t162) * t266;
t91 = -t162 * t268 - t166;
t84 = mrSges(5,1) * t268 + mrSges(5,3) * t124;
t83 = -mrSges(5,2) * t268 + mrSges(5,3) * t125;
t72 = -pkin(4) * t125 + t107;
t70 = -t121 * t187 + t130;
t60 = pkin(2) * t137 - qJD(1) * t264 + t220;
t53 = t110 * t195 + t266 * t348;
t52 = qJD(2) * t207 - qJD(5) * t98;
t48 = mrSges(6,1) * t157 + mrSges(6,3) * t66;
t47 = -mrSges(6,2) * t157 + mrSges(6,3) * t240;
t37 = pkin(7) * t192 * t263 + t44;
t35 = -pkin(4) * t87 + t65;
t33 = t88 * Ifges(5,1) + t87 * Ifges(5,4) + t138 * Ifges(5,5);
t31 = qJD(2) * t214 + t43;
t26 = Ifges(6,2) * t240 + t157 * Ifges(6,6) - t312;
t15 = -mrSges(6,2) * t126 + mrSges(6,3) * t24;
t14 = mrSges(6,1) * t126 - mrSges(6,3) * t23;
t4 = -qJD(5) * t19 - t191 * t37 + t194 * t31;
t3 = qJD(5) * t18 + t191 * t31 + t194 * t37;
t5 = [(Ifges(5,5) * t192 - t195 * t229) * t320 + (-t96 * mrSges(6,1) + t95 * mrSges(6,2) + (-m(4) * t213 + t368 * t195 + (m(5) * qJ(3) + m(6) * t159 + t231) * t192 + (m(3) + t362) * pkin(1) - t340) * t193 + ((-m(3) - t262) * pkin(6) + t364) * t196) * g(1) + (-m(3) * t269 - t94 * mrSges(6,1) - t93 * mrSges(6,2) - t231 * t274 + (-m(5) * qJ(4) - mrSges(5,3)) * t271 - t262 * (pkin(2) * t271 + qJ(3) * t274 + t269) + t364 * t193 + (-m(6) * (pkin(4) * t279 - t275) - t289 + t340) * t196) * g(2) + t195 * (Ifges(3,4) * t138 - Ifges(3,2) * t137 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t195 * (Ifges(4,5) * qJDD(2) - Ifges(4,6) * t138 + Ifges(4,3) * t137) / 0.2e1 - t192 * (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t138 + Ifges(4,6) * t137) / 0.2e1 + t151 * t45 - pkin(1) * (mrSges(3,1) * t137 + mrSges(3,2) * t138) + t142 * (-mrSges(4,2) * t137 - mrSges(4,3) * t138) + t105 * t131 - t135 * t74 + t108 * t7 + t234 * t285 + (-t1 * t98 + t2 * t99 - t52 * t8 + t53 * t9) * mrSges(6,3) + (-Ifges(6,5) * t99 - Ifges(6,6) * t98 + Ifges(6,3) * t192) * t317 + t35 * (mrSges(6,1) * t98 - mrSges(6,2) * t99) + (-Ifges(6,4) * t99 - Ifges(6,2) * t98 + Ifges(6,6) * t192) * t329 + (-Ifges(6,1) * t99 - Ifges(6,4) * t98 + Ifges(6,5) * t192) * t330 + t92 * t28 + t44 * t83 + t43 * t84 + t71 * t58 + t72 * (-mrSges(6,1) * t53 + mrSges(6,2) * t52) + t70 * t59 + t53 * t26 / 0.2e1 + t3 * t47 + t4 * t48 + t19 * t15 + t18 * t14 + t192 * t370 + (t210 + t199 + t195 * (-Ifges(3,2) * t192 + t299)) * t261 / 0.2e1 - t137 * t228 / 0.2e1 + t60 * t230 + t137 * t222 / 0.2e1 + (Ifges(6,4) * t52 + Ifges(6,2) * t53) * t324 + t350 * t242 + m(5) * (-t107 * t135 + t151 * t65 + t16 * t70 + t17 * t71 + t41 * t43 + t42 * t44) - t192 * t371 + m(4) * (pkin(6) * t346 + t105 * t119 + t142 * t60) + (-t137 * t177 + t138 * t310 + t347) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t347) + (Ifges(6,5) * t52 + Ifges(6,6) * t53) * t314 + (-qJDD(2) * mrSges(3,2) - t103) * t177 + (Ifges(6,3) * t314 + t342 / 0.2e1 + t351 * pkin(6) + t141 * mrSges(4,1) + Ifges(6,5) * t322 + Ifges(6,6) * t324 - t367) * t265 + (t349 * t373 + t334) * qJD(2) + (-qJDD(2) * mrSges(3,1) + t104) * t310 + (t192 * t359 - t195 * t358) * qJDD(2) / 0.2e1 + (-Ifges(3,4) * t137 + Ifges(3,5) * qJDD(2) + Ifges(5,5) * t88 + Ifges(5,6) * t87 + t138 * t360 + t256) * t192 / 0.2e1 + (t192 * t360 - t195 * t224 + t299) * t316 + t65 * t232 * t195 + t346 * mrSges(4,1) + (t336 / 0.2e1 + t352 * pkin(6) + t145 * mrSges(4,1) - t115 / 0.2e1) * t266 + (Ifges(5,6) * t192 - t195 * t226) * t321 + t276 * t326 + t52 * t327 - t99 * t331 - t98 * t332 + t223 * t361 + (Ifges(6,1) * t52 + Ifges(6,4) * t53) * t322 + m(6) * (t1 * t19 + t108 * t35 + t18 * t2 + t3 * t9 + t4 * t8 + t72 * t92) - t212 * t261 + t17 * (-mrSges(5,2) * t192 - mrSges(5,3) * t276) - t33 * t278 / 0.2e1 + t16 * (mrSges(5,1) * t192 + mrSges(5,3) * t278) + Ifges(2,3) * qJDD(1); (Ifges(5,5) * t188 - Ifges(5,6) * t187) * t316 + (-qJ(3) * t262 * t271 + t196 * t339) * g(1) + (-t262 * t286 + t339) * t366 + t188 * t33 / 0.2e1 + t159 * t7 - t132 * t131 - t133 * t74 + t122 * mrSges(3,2) - t123 * mrSges(3,1) - t355 * t26 / 0.2e1 - t217 * t332 + (qJ(3) * t65 - t221 * t302 + (-t187 * t42 - t188 * t41) * qJD(4) - t107 * t133 - t41 * t56 - t42 * t57) * m(5) - t345 * t302 + (-t334 + (t350 / 0.2e1 - t199 / 0.2e1 - t210 / 0.2e1 + t212) * qJD(1)) * qJD(1) + t356 * t48 + (-t103 + t45) * qJ(3) + (Ifges(6,5) * t86 + Ifges(6,6) * t85) * t315 - t336 * t268 / 0.2e1 + t101 * mrSges(4,2) - pkin(2) * t104 - t90 * mrSges(4,3) - t91 * t28 - t57 * t83 - t56 * t84 - t86 * t27 / 0.2e1 + t75 * t14 + t76 * t15 + (Ifges(6,4) * t86 + Ifges(6,2) * t85) * t325 + (-m(5) * t235 - t195 * mrSges(5,3) - m(4) * t270 - m(6) * (t270 - t275) - t289 + t369 * t192 + t365) * g(3) + (-m(4) * t145 + m(5) * t107 + m(6) * t72 - t341) * qJD(3) - (-Ifges(3,2) * t268 + t165 + t342) * t267 / 0.2e1 + t348 * t331 + (Ifges(6,5) * t348 - Ifges(6,6) * t217) * t317 + t35 * (mrSges(6,1) * t217 + mrSges(6,2) * t348) + (Ifges(6,4) * t348 - Ifges(6,2) * t217) * t329 + (Ifges(6,1) * t348 - Ifges(6,4) * t217) * t330 + t363 * mrSges(6,3) + t65 * t231 + (mrSges(6,1) * t355 - mrSges(6,2) * t372) * t72 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) - t351 * t168 - t352 * t166 + (-Ifges(6,5) * t110 - Ifges(6,6) * t109) * t314 + (-Ifges(6,1) * t110 - Ifges(6,4) * t109) * t322 + (-Ifges(6,4) * t110 - Ifges(6,2) * t109) * t324 + (-t187 * t83 - t188 * t84) * qJD(4) + t344 * t233 + t349 * t242 + (-pkin(2) * t101 - qJ(3) * t90 - t119 * t132) * m(4) + (Ifges(6,5) * t323 + Ifges(6,6) * t325 + Ifges(6,3) * t315 + t367) * t267 + t357 * t47 + (t1 * t76 + t159 * t35 + t2 * t75 + t356 * t8 + t357 * t9 - t72 * t91) * m(6) + t358 * t137 + t359 * t138 + (Ifges(6,1) * t86 + Ifges(6,4) * t85) * t323 + (Ifges(5,1) * t188 - t298) * t320 + (-Ifges(5,2) * t187 + t297) * t321 + t187 * t326 - t110 * t327 - t145 * t252 - t141 * t253 + t115 * t268 / 0.2e1 - t16 * t292 - mrSges(5,3) * t294; t217 * t15 + t348 * t14 - t372 * t48 + t355 * t47 + t341 * qJD(2) + t262 * t307 + ((-t187 * t84 + t188 * t83 + t131) * qJD(1) - t344 * t262) * t192 + t104 + t345 + (-qJD(2) * t72 - t363) * m(6) + (-qJD(2) * t107 - (t41 * t187 - t42 * t188) * t268 + t221) * m(5) + (qJD(2) * t145 + t119 * t268 + t101) * m(4); -t362 * t192 * g(3) - t124 * t84 - t125 * t83 - t240 * t47 - t66 * t48 + t45 + t7 + (-t240 * t9 - t66 * t8 + t35 - t353) * m(6) + (-t124 * t41 - t125 * t42 - t353 + t65) * m(5); -t371 + t370 - t72 * (-mrSges(6,1) * t66 + mrSges(6,2) * t240) + (Ifges(6,1) * t240 + t312) * t323 + t26 * t322 + (Ifges(6,5) * t240 + Ifges(6,6) * t66) * t315 - t8 * t47 + t9 * t48 - g(1) * (mrSges(6,1) * t93 - mrSges(6,2) * t94) - g(2) * (mrSges(6,1) * t95 + mrSges(6,2) * t96) - (-mrSges(6,1) * t173 + mrSges(6,2) * t172) * t307 + (t240 * t8 - t66 * t9) * mrSges(6,3) + t256 + (Ifges(6,2) * t66 + t27 + t64) * t325;];
tau = t5;
