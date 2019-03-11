% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:48:32
% EndTime: 2019-03-09 01:48:48
% DurationCPUTime: 11.70s
% Computational Cost: add. (4112->537), mult. (7706->726), div. (0->0), fcn. (4588->10), ass. (0->244)
t174 = sin(qJ(4));
t177 = cos(qJ(4));
t244 = qJD(1) * qJD(4);
t134 = qJDD(1) * t174 + t177 * t244;
t292 = t134 / 0.2e1;
t175 = sin(qJ(1));
t352 = g(1) * t175;
t178 = cos(qJ(1));
t351 = g(2) * t178;
t168 = sin(pkin(9));
t169 = cos(pkin(9));
t133 = qJDD(1) * t177 - t174 * t244;
t172 = pkin(1) + qJ(3);
t198 = qJDD(1) * t172 - qJDD(2);
t219 = -qJD(5) * t177 + qJD(3);
t41 = pkin(4) * t134 - qJ(5) * t133 + qJD(1) * t219 + t198;
t164 = qJD(1) * qJD(2);
t321 = qJDD(1) * qJ(2) + t164;
t144 = qJDD(3) + t321;
t129 = -qJDD(1) * pkin(7) + t144;
t154 = qJD(1) * qJ(2) + qJD(3);
t146 = -qJD(1) * pkin(7) + t154;
t247 = qJD(4) * t177;
t74 = t174 * t129 + t146 * t247;
t60 = qJDD(4) * qJ(5) + qJD(4) * qJD(5) + t74;
t16 = -t168 * t60 + t169 * t41;
t90 = qJDD(4) * t168 + t133 * t169;
t10 = pkin(5) * t134 - pkin(8) * t90 + t16;
t17 = t168 * t41 + t169 * t60;
t89 = qJDD(4) * t169 - t133 * t168;
t11 = pkin(8) * t89 + t17;
t173 = sin(qJ(6));
t176 = cos(qJ(6));
t251 = qJD(1) * t177;
t123 = qJD(4) * t168 + t169 * t251;
t246 = t174 * qJD(1);
t135 = t174 * t146;
t113 = qJD(4) * qJ(5) + t135;
t218 = qJ(5) * t177 - t172;
t136 = pkin(4) * t174 - t218;
t92 = qJD(1) * t136 - qJD(2);
t45 = -t113 * t168 + t169 * t92;
t29 = pkin(5) * t246 - pkin(8) * t123 + t45;
t122 = qJD(4) * t169 - t168 * t251;
t46 = t169 * t113 + t168 * t92;
t32 = pkin(8) * t122 + t46;
t8 = -t173 * t32 + t176 * t29;
t1 = qJD(6) * t8 + t10 * t173 + t11 * t176;
t9 = t173 * t29 + t176 * t32;
t2 = -qJD(6) * t9 + t10 * t176 - t11 * t173;
t350 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t349 = t174 / 0.2e1;
t151 = qJD(6) + t246;
t290 = t151 / 0.2e1;
t67 = t122 * t173 + t123 * t176;
t298 = t67 / 0.2e1;
t220 = t176 * t122 - t123 * t173;
t300 = t220 / 0.2e1;
t348 = Ifges(7,5) * t298 + Ifges(7,6) * t300 + Ifges(7,3) * t290;
t279 = Ifges(5,4) * t177;
t347 = Ifges(5,2) * t349 - t279 / 0.2e1;
t228 = t246 / 0.2e1;
t213 = -mrSges(6,1) * t169 + mrSges(6,2) * t168;
t191 = m(6) * pkin(4) - t213;
t214 = t174 * mrSges(5,1) + t177 * mrSges(5,2);
t346 = t191 * t174 + t214;
t233 = m(6) * qJ(5) + mrSges(6,3);
t281 = pkin(8) + qJ(5);
t345 = -m(7) * t281 - mrSges(7,3) - t233;
t288 = t177 / 0.2e1;
t344 = 0.2e1 * t288;
t323 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t122 - mrSges(6,2) * t123 - mrSges(5,3) * t251;
t153 = pkin(5) * t169 + pkin(4);
t333 = mrSges(3,2) - mrSges(4,3);
t341 = m(7) * (t153 * t174 - t177 * t281) + mrSges(2,1) - t177 * mrSges(7,3) - t333 + t346;
t170 = qJ(2) - pkin(7);
t250 = qJD(2) * t174;
t317 = t170 * t247 + t250;
t217 = g(1) * t178 + g(2) * t175;
t253 = t176 * t169;
t125 = t168 * t173 - t253;
t340 = qJD(6) * t125;
t296 = t90 / 0.2e1;
t339 = Ifges(6,1) * t296 + Ifges(6,5) * t292;
t338 = t8 * mrSges(7,1) - Ifges(5,6) * qJD(4) / 0.2e1 + qJD(1) * t347 + t123 * Ifges(6,5) / 0.2e1 + t122 * Ifges(6,6) / 0.2e1 + Ifges(6,3) * t228 - t9 * mrSges(7,2) + t348;
t124 = qJDD(6) + t134;
t21 = qJD(6) * t220 + t173 * t89 + t176 * t90;
t22 = -qJD(6) * t67 - t173 * t90 + t176 * t89;
t337 = -t21 * Ifges(7,4) / 0.2e1 - t22 * Ifges(7,2) / 0.2e1 - t124 * Ifges(7,6) / 0.2e1;
t306 = t21 / 0.2e1;
t305 = t22 / 0.2e1;
t297 = t89 / 0.2e1;
t336 = -m(6) - m(4);
t335 = -m(7) - m(5);
t293 = t124 / 0.2e1;
t334 = qJD(4) / 0.2e1;
t332 = -mrSges(4,2) - mrSges(3,3);
t259 = t169 * t174;
t242 = pkin(8) * t259;
t205 = pkin(4) * t177 + qJ(5) * t174;
t131 = t205 * qJD(1);
t260 = t168 * t177;
t70 = t169 * t131 - t146 * t260;
t42 = (pkin(5) * t177 + t242) * qJD(1) + t70;
t237 = t168 * t246;
t258 = t169 * t177;
t71 = t168 * t131 + t146 * t258;
t54 = pkin(8) * t237 + t71;
t139 = t281 * t168;
t140 = t281 * t169;
t75 = -t139 * t176 - t140 * t173;
t331 = -qJD(5) * t125 + qJD(6) * t75 - t173 * t42 - t176 * t54;
t126 = t168 * t176 + t169 * t173;
t76 = -t139 * t173 + t140 * t176;
t330 = -qJD(5) * t126 - qJD(6) * t76 + t173 * t54 - t176 * t42;
t38 = -t89 * mrSges(6,1) + t90 * mrSges(6,2);
t328 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t133 - t38;
t99 = t126 * t177;
t327 = t125 * qJD(1) - qJD(4) * t99 + t174 * t340;
t101 = t125 * t177;
t108 = t126 * qJD(1);
t313 = t126 * qJD(6);
t326 = -qJD(4) * t101 - t174 * t313 - t108;
t88 = -t173 * t237 + t246 * t253;
t325 = t340 - t88;
t87 = t174 * t108;
t324 = t313 + t87;
t322 = t177 * t217;
t277 = Ifges(6,4) * t169;
t208 = -Ifges(6,2) * t168 + t277;
t278 = Ifges(6,4) * t168;
t210 = Ifges(6,1) * t169 - t278;
t320 = t122 * (Ifges(6,6) * t177 - t174 * t208) + t123 * (Ifges(6,5) * t177 - t174 * t210);
t206 = Ifges(6,5) * t169 - Ifges(6,6) * t168;
t319 = t177 * (-Ifges(5,1) * t174 - t279) + t174 * (Ifges(6,3) * t177 - t174 * t206);
t248 = qJD(4) * t174;
t73 = t129 * t177 - t146 * t248;
t316 = -t174 * t74 - t177 * t73;
t204 = -t16 * t168 + t169 * t17;
t84 = -mrSges(6,2) * t246 + mrSges(6,3) * t122;
t85 = mrSges(6,1) * t246 - mrSges(6,3) * t123;
t315 = -t168 * t85 + t169 * t84;
t55 = -mrSges(6,2) * t134 + mrSges(6,3) * t89;
t56 = mrSges(6,1) * t134 - mrSges(6,3) * t90;
t314 = -t168 * t56 + t169 * t55;
t147 = qJD(1) * t172 - qJD(2);
t212 = t168 * mrSges(6,1) + t169 * mrSges(6,2);
t215 = mrSges(5,1) * t177 - mrSges(5,2) * t174;
t261 = t168 * t174;
t263 = t146 * t177;
t105 = -qJD(4) * pkin(4) + qJD(5) - t263;
t265 = t105 * t174;
t311 = t212 * t265 - t46 * (-mrSges(6,2) * t177 + mrSges(6,3) * t261) - t45 * (mrSges(6,1) * t177 + mrSges(6,3) * t259) - t147 * t215;
t286 = pkin(5) * t168;
t310 = m(6) * pkin(7) + m(7) * t286 + mrSges(2,2) + mrSges(5,3) + t212 + t332;
t179 = qJD(1) ^ 2;
t308 = -m(5) / 0.2e1;
t307 = Ifges(7,1) * t306 + Ifges(7,4) * t305 + Ifges(7,5) * t293;
t62 = Ifges(7,4) * t220;
t25 = Ifges(7,1) * t67 + Ifges(7,5) * t151 + t62;
t303 = t25 / 0.2e1;
t302 = Ifges(6,4) * t297 + t339;
t301 = -t220 / 0.2e1;
t299 = -t67 / 0.2e1;
t295 = -m(4) - m(5);
t294 = -m(6) - m(7);
t291 = -t151 / 0.2e1;
t287 = Ifges(7,4) * t67;
t280 = Ifges(5,4) * t174;
t63 = -qJDD(4) * pkin(4) + qJDD(5) - t73;
t268 = t177 * t63;
t266 = qJDD(1) * pkin(1);
t257 = t170 * t174;
t255 = t174 * t175;
t254 = t174 * t178;
t78 = t168 * t136 + t169 * t257;
t252 = t178 * pkin(1) + t175 * qJ(2);
t249 = qJD(2) * t177;
t245 = qJD(1) * qJD(3);
t28 = -mrSges(7,1) * t220 + mrSges(7,2) * t67;
t241 = t28 - t323;
t240 = Ifges(7,5) * t21 + Ifges(7,6) * t22 + Ifges(7,3) * t124;
t102 = qJD(4) * t205 + t219;
t52 = t168 * t102 + t169 * t317;
t239 = t294 + t295;
t238 = t178 * qJ(3) + t252;
t236 = t168 * t248;
t7 = -t22 * mrSges(7,1) + t21 * mrSges(7,2);
t227 = -t168 * t170 + pkin(5);
t226 = -t170 + t286;
t224 = -t244 / 0.2e1;
t222 = (t174 ^ 2 + t177 ^ 2) * t146;
t216 = t134 * mrSges(5,1) + t133 * mrSges(5,2);
t211 = t177 * Ifges(5,1) - t280;
t207 = -Ifges(5,5) * t174 - Ifges(5,6) * t177;
t203 = -t168 * t45 + t169 * t46;
t118 = t169 * t136;
t61 = -pkin(8) * t258 + t174 * t227 + t118;
t68 = -pkin(8) * t260 + t78;
t26 = -t173 * t68 + t176 * t61;
t27 = t173 * t61 + t176 * t68;
t130 = t198 + t245;
t200 = qJD(3) * t147 + t130 * t172;
t141 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t246;
t197 = -t141 - t315;
t195 = t174 * (-Ifges(5,2) * t177 - t280);
t162 = pkin(9) + qJ(6);
t156 = sin(t162);
t157 = cos(t162);
t189 = m(7) * t153 + mrSges(7,1) * t157 - mrSges(7,2) * t156;
t160 = t178 * qJ(2);
t155 = qJDD(2) - t266;
t132 = t214 * qJD(1);
t121 = t226 * t177;
t116 = Ifges(5,5) * qJD(4) + qJD(1) * t211;
t104 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t134;
t100 = t125 * t174;
t98 = t126 * t174;
t96 = -t156 * t175 + t157 * t254;
t95 = -t156 * t254 - t157 * t175;
t94 = -t156 * t178 - t157 * t255;
t93 = t156 * t255 - t157 * t178;
t91 = -pkin(5) * t237 + t135;
t86 = -t226 * t248 - t249;
t83 = t169 * t102;
t77 = -t168 * t257 + t118;
t69 = -pkin(5) * t122 + t105;
t59 = t123 * Ifges(6,1) + t122 * Ifges(6,4) + Ifges(6,5) * t246;
t58 = t123 * Ifges(6,4) + t122 * Ifges(6,2) + Ifges(6,6) * t246;
t51 = -t168 * t317 + t83;
t50 = t126 * t248 + t177 * t340;
t48 = t125 * t248 - t177 * t313;
t44 = mrSges(7,1) * t151 - mrSges(7,3) * t67;
t43 = -mrSges(7,2) * t151 + mrSges(7,3) * t220;
t35 = pkin(8) * t236 + t52;
t34 = -t168 * t250 + t83 + (t177 * t227 + t242) * qJD(4);
t33 = -pkin(5) * t89 + t63;
t30 = t90 * Ifges(6,4) + t89 * Ifges(6,2) + t134 * Ifges(6,6);
t24 = Ifges(7,2) * t220 + Ifges(7,6) * t151 + t287;
t13 = -mrSges(7,2) * t124 + mrSges(7,3) * t22;
t12 = mrSges(7,1) * t124 - mrSges(7,3) * t21;
t6 = -qJD(6) * t27 - t173 * t35 + t176 * t34;
t5 = qJD(6) * t26 + t173 * t34 + t176 * t35;
t3 = [t33 * (mrSges(7,1) * t99 - mrSges(7,2) * t101) + (-t16 * t258 - t17 * t260) * mrSges(6,3) + (-t94 * mrSges(7,1) - t93 * mrSges(7,2) + t335 * (-pkin(7) * t178 + t160) + (-m(3) + t336) * t160 + t310 * t178 + (m(3) * pkin(1) - m(6) * t218 + (m(4) - t335) * t172 + t341) * t175) * g(1) + (t155 - t266) * mrSges(3,2) + (-Ifges(7,1) * t101 - Ifges(7,4) * t99) * t306 - (t169 * t59 + t116) * t248 / 0.2e1 + 0.2e1 * t321 * mrSges(3,3) + (Ifges(7,5) * t48 + Ifges(7,6) * t50) * t290 + t172 * t216 + m(7) * (t1 * t27 + t121 * t33 + t2 * t26 + t5 * t9 + t6 * t8 + t69 * t86) + t133 * t211 / 0.2e1 + m(4) * (qJ(2) * t144 + qJD(2) * t154 + t200) + (Ifges(7,1) * t48 + Ifges(7,4) * t50) * t298 + (Ifges(7,4) * t48 + Ifges(7,2) * t50) * t300 + t58 * t236 / 0.2e1 + t316 * mrSges(5,3) + t317 * t141 + t319 * t244 / 0.2e1 + (t144 + t321) * mrSges(4,2) + (t130 + t245) * mrSges(4,3) + m(3) * (-pkin(1) * t155 + (t321 + t164) * qJ(2)) + m(5) * (qJD(2) * t222 - t170 * t316 + t200) + (-Ifges(7,4) * t101 - Ifges(7,2) * t99) * t305 + qJD(3) * t132 + t121 * t7 + t52 * t84 + t51 * t85 + t86 * t28 + t77 * t56 + t78 * t55 + t69 * (-mrSges(7,1) * t50 + mrSges(7,2) * t48) + t50 * t24 / 0.2e1 + t5 * t43 + t6 * t44 + t26 * t12 + t27 * t13 + m(6) * (t16 * t77 + t17 * t78 + t45 * t51 + t46 * t52) + (-mrSges(6,3) * t352 + t233 * t351 + (-m(6) * t63 + t328) * t170 + t210 * t296 + t208 * t297 + t206 * t292) * t177 + (-Ifges(7,5) * t101 - Ifges(7,6) * t99) * t293 + Ifges(5,5) * t344 * qJDD(4) + (-t1 * t99 + t101 * t2 - t48 * t8 + t50 * t9) * mrSges(7,3) + (t338 + t348) * t247 + (Ifges(7,5) * t306 + Ifges(7,3) * t293 + Ifges(7,6) * t305 - Ifges(5,4) * t133 / 0.2e1 + Ifges(6,5) * t296 + Ifges(6,6) * t297 + t16 * mrSges(6,1) - t17 * mrSges(6,2) - Ifges(5,6) * qJDD(4) + t350 + (Ifges(5,2) + Ifges(6,3)) * t292) * t174 + (Ifges(5,1) * t133 - Ifges(5,4) * t134) * t288 + t130 * t214 + t134 * t347 + t258 * t302 + t48 * t303 - t101 * t307 + (Ifges(6,5) * t90 + Ifges(6,6) * t89 + Ifges(6,3) * t134 + t240) * t349 + t104 * t257 + t212 * t268 + t99 * t337 + t320 * t334 + t195 * t224 - t30 * t260 / 0.2e1 + (m(6) * t105 - t323) * (t170 * t248 - t249) + (-m(3) * t252 - t96 * mrSges(7,1) - t95 * mrSges(7,2) + t336 * t238 + t335 * (-pkin(7) * t175 + t238) + t310 * t175 - t341 * t178) * g(2) + (t207 * t334 - t311) * qJD(4) + (mrSges(4,3) * t172 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1); t125 * t12 - t126 * t13 - t168 * t55 - t169 * t56 + t324 * t44 + t325 * t43 + t333 * qJDD(1) + (-m(3) * qJ(2) + t332) * t179 + m(3) * t155 + m(6) * (-t16 * t169 - t168 * t17) + 0.2e1 * (-m(4) / 0.2e1 + t308) * t130 - t216 + (t351 - t352) * (m(3) - t239) + (-t1 * t126 + t125 * t2 + t324 * t8 + t325 * t9) * m(7) + (-m(4) * t154 + t197 * t174 + t241 * t177 + m(7) * t69 * t344 - m(6) * (-t105 * t177 + t259 * t46 - t261 * t45) + 0.2e1 * t222 * t308) * qJD(1); qJDD(1) * mrSges(4,2) - t179 * mrSges(4,3) - t100 * t13 - t98 * t12 + t327 * t44 + t326 * t43 + (-t7 + t328) * t177 + (t104 + t314) * t174 + (t147 * t295 - t168 * t84 - t169 * t85 - t132) * qJD(1) + (t174 * t241 - t177 * t197) * qJD(4) + m(4) * t144 - m(5) * t316 + t217 * t239 + (-t1 * t100 - t177 * t33 - t2 * t98 + t69 * t248 + t326 * t9 + t327 * t8) * m(7) + (-t268 + t204 * t174 + (t177 * t203 + t265) * qJD(4) - (t168 * t46 + t169 * t45) * qJD(1)) * m(6); (-t320 / 0.2e1 + t311) * qJD(1) + (-Ifges(7,1) * t88 + Ifges(7,4) * t87) * t299 + (-t319 / 0.2e1 + t195 / 0.2e1) * t179 + t323 * t135 - t324 * t24 / 0.2e1 + (-t1 * t125 - t126 * t2 - t324 * t9 + t325 * t8) * mrSges(7,3) + (mrSges(7,1) * t324 - mrSges(7,2) * t325) * t69 + (t302 + t339) * t168 + t63 * t213 + t277 * t296 - t58 * t237 / 0.2e1 + t314 * qJ(5) + t315 * qJD(5) + t204 * mrSges(6,3) - t340 * t25 / 0.2e1 + (-Ifges(7,1) * t340 - Ifges(7,4) * t313) * t298 + (-Ifges(7,4) * t340 - Ifges(7,2) * t313) * t300 + (-Ifges(7,5) * t340 - Ifges(7,6) * t313) * t290 + (-Ifges(7,4) * t88 + Ifges(7,2) * t87) * t301 - t153 * t7 - Ifges(5,6) * t134 + t33 * (mrSges(7,1) * t125 + mrSges(7,2) * t126) + Ifges(5,5) * t133 - t91 * t28 - t71 * t84 - t70 * t85 + t73 * mrSges(5,1) - t74 * mrSges(5,2) + t75 * t12 + t76 * t13 - pkin(4) * t38 + t278 * t297 + (Ifges(6,2) * t297 + Ifges(6,6) * t292 + t59 * t228 + t30 / 0.2e1) * t169 - t141 * t263 + (-pkin(4) * t63 + qJ(5) * t204 + qJD(5) * t203 - t105 * t135 - t45 * t70 - t46 * t71) * m(6) + (Ifges(7,5) * t299 + Ifges(7,6) * t301 + Ifges(7,3) * t291 - t338) * t251 + (-Ifges(7,5) * t88 + Ifges(7,6) * t87) * t291 + t217 * (-t215 + (-t189 - t191) * t177 + t345 * t174) + (t174 * t189 + t345 * t177 + t346) * g(3) + t88 * t303 + (Ifges(7,4) * t126 - Ifges(7,2) * t125) * t305 + (Ifges(7,1) * t126 - Ifges(7,4) * t125) * t306 + t126 * t307 + (Ifges(7,5) * t126 - Ifges(7,6) * t125) * t293 + t125 * t337 + t207 * t224 + t116 * t228 + Ifges(5,3) * qJDD(4) + t330 * t44 + t331 * t43 + (t1 * t76 - t153 * t33 + t2 * t75 + t330 * t8 + t331 * t9 - t69 * t91) * m(7); t294 * t174 * g(3) - t122 * t84 + t123 * t85 - t220 * t43 + t67 * t44 + t38 + t7 + (-t220 * t9 + t67 * t8 + t322 + t33) * m(7) + (-t122 * t46 + t123 * t45 + t322 + t63) * m(6); -t69 * (mrSges(7,1) * t67 + mrSges(7,2) * t220) + (Ifges(7,1) * t220 - t287) * t299 + t24 * t298 + (Ifges(7,5) * t220 - Ifges(7,6) * t67) * t291 - t8 * t43 + t9 * t44 - g(1) * (mrSges(7,1) * t95 - mrSges(7,2) * t96) - g(2) * (-mrSges(7,1) * t93 + mrSges(7,2) * t94) - g(3) * (-mrSges(7,1) * t156 - mrSges(7,2) * t157) * t177 + (t220 * t8 + t67 * t9) * mrSges(7,3) + t240 + (-Ifges(7,2) * t67 + t25 + t62) * t301 + t350;];
tau  = t3;
