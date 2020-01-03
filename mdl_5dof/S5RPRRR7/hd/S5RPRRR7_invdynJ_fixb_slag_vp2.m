% Calculate vector of inverse dynamics joint torques for
% S5RPRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:12
% EndTime: 2019-12-31 19:03:34
% DurationCPUTime: 12.26s
% Computational Cost: add. (4925->563), mult. (10431->776), div. (0->0), fcn. (6614->14), ass. (0->258)
t361 = mrSges(5,3) + mrSges(6,3);
t182 = sin(qJ(4));
t187 = cos(qJ(3));
t259 = qJD(1) * t187;
t238 = t182 * t259;
t189 = -pkin(8) - pkin(7);
t239 = qJD(4) * t189;
t179 = sin(pkin(9));
t164 = pkin(1) * t179 + pkin(6);
t152 = t164 * qJD(1);
t183 = sin(qJ(3));
t111 = qJD(2) * t187 - t183 * t152;
t224 = pkin(3) * t183 - pkin(7) * t187;
t144 = t224 * qJD(1);
t186 = cos(qJ(4));
t66 = t186 * t111 + t182 * t144;
t360 = pkin(8) * t238 + t182 * t239 - t66;
t262 = t186 * t187;
t206 = pkin(4) * t183 - pkin(8) * t262;
t65 = -t111 * t182 + t186 * t144;
t359 = -qJD(1) * t206 + t186 * t239 - t65;
t169 = pkin(4) * t186 + pkin(3);
t178 = qJ(4) + qJ(5);
t174 = sin(t178);
t175 = cos(t178);
t221 = -mrSges(5,1) * t186 + mrSges(5,2) * t182;
t358 = m(5) * pkin(3) + m(6) * t169 + mrSges(6,1) * t175 - mrSges(6,2) * t174 - t221;
t357 = -m(5) * pkin(7) + m(6) * t189 - t361;
t250 = qJD(1) * qJD(3);
t148 = qJDD(1) * t187 - t183 * t250;
t356 = t148 / 0.2e1;
t149 = qJDD(1) * t183 + t187 * t250;
t355 = t149 / 0.2e1;
t354 = t250 / 0.2e1;
t353 = qJD(2) * qJD(3) + t164 * qJDD(1);
t258 = qJD(2) * t183;
t112 = t152 * t187 + t258;
t102 = qJD(3) * pkin(7) + t112;
t225 = pkin(3) * t187 + pkin(7) * t183;
t207 = -pkin(2) - t225;
t180 = cos(pkin(9));
t304 = pkin(1) * t180;
t131 = t207 - t304;
t103 = t131 * qJD(1);
t251 = qJD(4) * t186;
t253 = qJD(4) * t182;
t256 = qJD(3) * t183;
t60 = t183 * qJDD(2) - t152 * t256 + t187 * t353;
t55 = qJDD(3) * pkin(7) + t60;
t165 = -pkin(2) - t304;
t151 = t165 * qJDD(1);
t70 = -pkin(3) * t148 - pkin(7) * t149 + t151;
t13 = -t102 * t253 + t103 * t251 + t182 * t70 + t186 * t55;
t50 = t102 * t186 + t103 * t182;
t14 = -qJD(4) * t50 - t182 * t55 + t186 * t70;
t210 = t13 * t186 - t14 * t182;
t49 = -t102 * t182 + t186 * t103;
t352 = -t49 * t251 - t50 * t253 + t210;
t257 = qJD(3) * t182;
t260 = qJD(1) * t183;
t140 = t186 * t260 + t257;
t73 = -qJD(4) * t140 + qJDD(3) * t186 - t149 * t182;
t10 = pkin(8) * t73 + t13;
t185 = cos(qJ(5));
t181 = sin(qJ(5));
t255 = qJD(3) * t186;
t139 = -t182 * t260 + t255;
t39 = pkin(8) * t139 + t50;
t281 = t181 * t39;
t160 = qJD(4) - t259;
t38 = -pkin(8) * t140 + t49;
t36 = pkin(4) * t160 + t38;
t11 = t185 * t36 - t281;
t135 = qJDD(4) - t148;
t72 = qJD(4) * t139 + qJDD(3) * t182 + t149 * t186;
t9 = pkin(4) * t135 - pkin(8) * t72 + t14;
t2 = qJD(5) * t11 + t10 * t185 + t181 * t9;
t276 = t185 * t39;
t12 = t181 * t36 + t276;
t3 = -qJD(5) * t12 - t10 * t181 + t185 * t9;
t351 = -t3 * mrSges(6,1) + t2 * mrSges(6,2);
t350 = -t14 * mrSges(5,1) + t13 * mrSges(5,2);
t289 = Ifges(4,4) * t183;
t216 = t187 * Ifges(4,2) + t289;
t349 = t12 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t216 / 0.2e1 - t11 * mrSges(6,1);
t327 = m(6) * pkin(4);
t226 = t185 * t139 - t140 * t181;
t21 = qJD(5) * t226 + t181 * t73 + t185 * t72;
t326 = t21 / 0.2e1;
t78 = t139 * t181 + t140 * t185;
t22 = -qJD(5) * t78 - t181 * t72 + t185 * t73;
t325 = t22 / 0.2e1;
t318 = t72 / 0.2e1;
t317 = t73 / 0.2e1;
t129 = qJDD(5) + t135;
t312 = t129 / 0.2e1;
t311 = t135 / 0.2e1;
t347 = -qJD(1) / 0.2e1;
t346 = mrSges(3,2) - mrSges(4,3);
t157 = t189 * t182;
t158 = t189 * t186;
t87 = t157 * t181 - t158 * t185;
t345 = -qJD(5) * t87 - t360 * t181 + t359 * t185;
t86 = t157 * t185 + t158 * t181;
t344 = qJD(5) * t86 + t359 * t181 + t360 * t185;
t159 = qJD(5) + t160;
t343 = t140 * Ifges(5,5) + t78 * Ifges(6,5) + t139 * Ifges(5,6) + Ifges(6,6) * t226 + t160 * Ifges(5,3) + t159 * Ifges(6,3);
t342 = mrSges(5,1) + t327;
t35 = -mrSges(5,1) * t73 + mrSges(5,2) * t72;
t341 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t149 + t35;
t244 = mrSges(4,3) * t260;
t340 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t139 + mrSges(5,2) * t140 + t244;
t301 = pkin(4) * t182;
t339 = pkin(4) * t253 - t258 - (qJD(1) * t301 + t152) * t187;
t208 = t181 * t182 - t185 * t186;
t114 = t208 * t183;
t170 = Ifges(4,4) * t259;
t132 = Ifges(5,4) * t139;
t64 = Ifges(5,1) * t140 + Ifges(5,5) * t160 + t132;
t338 = Ifges(4,1) * t260 + Ifges(4,5) * qJD(3) + t186 * t64 + t170;
t143 = t164 * t262;
t80 = t182 * t131 + t143;
t46 = mrSges(5,1) * t135 - mrSges(5,3) * t72;
t47 = -mrSges(5,2) * t135 + mrSges(5,3) * t73;
t336 = -t182 * t46 + t186 * t47;
t254 = qJD(3) * t187;
t61 = qJDD(2) * t187 - t152 * t254 - t183 * t353;
t335 = -t183 * t61 + t187 * t60;
t334 = -m(5) - m(6) - m(4);
t333 = qJD(4) + qJD(5);
t223 = mrSges(4,1) * t187 - mrSges(4,2) * t183;
t331 = t361 * t183 + mrSges(3,1) + t223;
t101 = -qJD(3) * pkin(3) - t111;
t330 = -m(5) * t101 - t340;
t329 = Ifges(6,4) * t326 + Ifges(6,2) * t325 + Ifges(6,6) * t312;
t328 = Ifges(6,1) * t326 + Ifges(6,4) * t325 + Ifges(6,5) * t312;
t324 = Ifges(5,1) * t318 + Ifges(5,4) * t317 + Ifges(5,5) * t311;
t305 = Ifges(6,4) * t78;
t31 = Ifges(6,2) * t226 + Ifges(6,6) * t159 + t305;
t323 = -t31 / 0.2e1;
t322 = t31 / 0.2e1;
t74 = Ifges(6,4) * t226;
t32 = Ifges(6,1) * t78 + Ifges(6,5) * t159 + t74;
t321 = -t32 / 0.2e1;
t320 = t32 / 0.2e1;
t287 = Ifges(5,4) * t140;
t63 = Ifges(5,2) * t139 + Ifges(5,6) * t160 + t287;
t319 = -t63 / 0.2e1;
t316 = -t226 / 0.2e1;
t315 = t226 / 0.2e1;
t314 = -t78 / 0.2e1;
t313 = t78 / 0.2e1;
t309 = t140 / 0.2e1;
t308 = -t159 / 0.2e1;
t307 = t159 / 0.2e1;
t184 = sin(qJ(1));
t303 = pkin(1) * t184;
t302 = pkin(4) * t140;
t298 = g(3) * t183;
t297 = t11 * mrSges(6,3);
t296 = t12 * mrSges(6,3);
t188 = cos(qJ(1));
t176 = t188 * pkin(1);
t177 = qJ(1) + pkin(9);
t171 = sin(t177);
t172 = cos(t177);
t267 = t174 * t187;
t96 = t171 * t267 + t172 * t175;
t266 = t175 * t187;
t97 = -t171 * t266 + t172 * t174;
t293 = -t96 * mrSges(6,1) + t97 * mrSges(6,2);
t98 = t171 * t175 - t172 * t267;
t99 = t171 * t174 + t172 * t266;
t292 = t98 * mrSges(6,1) - t99 * mrSges(6,2);
t291 = mrSges(5,3) * t139;
t290 = mrSges(5,3) * t140;
t288 = Ifges(4,4) * t187;
t286 = Ifges(5,4) * t182;
t285 = Ifges(5,4) * t186;
t269 = t171 * t182;
t268 = t172 * t182;
t265 = t182 * t183;
t264 = t182 * t187;
t263 = t183 * t186;
t147 = t224 * qJD(3);
t237 = t164 * t256;
t261 = t186 * t147 + t182 * t237;
t252 = qJD(4) * t183;
t248 = Ifges(6,5) * t21 + Ifges(6,6) * t22 + Ifges(6,3) * t129;
t247 = Ifges(5,5) * t72 + Ifges(5,6) * t73 + Ifges(5,3) * t135;
t243 = mrSges(4,3) * t259;
t236 = t164 * t254;
t235 = t182 * t254;
t222 = mrSges(4,1) * t183 + mrSges(4,2) * t187;
t220 = mrSges(5,1) * t182 + mrSges(5,2) * t186;
t219 = -mrSges(6,1) * t174 - mrSges(6,2) * t175;
t218 = Ifges(5,1) * t186 - t286;
t217 = Ifges(5,1) * t182 + t285;
t215 = -Ifges(5,2) * t182 + t285;
t214 = Ifges(5,2) * t186 + t286;
t213 = Ifges(4,5) * t187 - Ifges(4,6) * t183;
t212 = Ifges(5,5) * t186 - Ifges(5,6) * t182;
t211 = Ifges(5,5) * t182 + Ifges(5,6) * t186;
t116 = t186 * t131;
t57 = -pkin(8) * t263 + t116 + (-t164 * t182 - pkin(4)) * t187;
t67 = -pkin(8) * t265 + t80;
t27 = -t181 * t67 + t185 * t57;
t28 = t181 * t57 + t185 * t67;
t209 = t169 * t187 - t183 * t189;
t142 = t181 * t186 + t182 * t185;
t205 = t248 - t351;
t108 = t171 * t186 - t172 * t264;
t106 = t171 * t264 + t172 * t186;
t203 = t165 * qJD(1) * t222;
t202 = t183 * (Ifges(4,1) * t187 - t289);
t201 = t142 * t187;
t200 = t208 * t187;
t197 = -t182 * t252 + t186 * t254;
t196 = t183 * t251 + t235;
t56 = -qJDD(3) * pkin(3) - t61;
t195 = Ifges(5,5) * t183 + t187 * t218;
t194 = Ifges(5,6) * t183 + t187 * t215;
t193 = Ifges(5,3) * t183 + t187 * t212;
t40 = t131 * t251 + t182 * t147 + (-t183 * t255 - t187 * t253) * t164;
t82 = t333 * t142;
t155 = -qJD(3) * mrSges(4,2) + t243;
t130 = t220 * t183;
t119 = (t164 + t301) * t183;
t117 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t148;
t113 = t142 * t183;
t109 = t172 * t262 + t269;
t107 = -t171 * t262 + t268;
t105 = qJD(1) * t200;
t104 = qJD(1) * t201;
t95 = mrSges(5,1) * t160 - t290;
t94 = -mrSges(5,2) * t160 + t291;
t84 = pkin(4) * t196 + t236;
t81 = t333 * t208;
t79 = -t164 * t264 + t116;
t71 = -pkin(4) * t139 + t101;
t54 = mrSges(6,1) * t159 - mrSges(6,3) * t78;
t53 = -mrSges(6,2) * t159 + mrSges(6,3) * t226;
t43 = -qJD(3) * t201 + t114 * t333;
t42 = -qJD(3) * t200 - t183 * t82;
t41 = -qJD(4) * t80 + t261;
t37 = -mrSges(6,1) * t226 + mrSges(6,2) * t78;
t34 = -pkin(8) * t196 + t40;
t33 = -pkin(4) * t73 + t56;
t29 = t206 * qJD(3) + (-t143 + (pkin(8) * t183 - t131) * t182) * qJD(4) + t261;
t25 = t72 * Ifges(5,4) + t73 * Ifges(5,2) + t135 * Ifges(5,6);
t18 = -mrSges(6,2) * t129 + mrSges(6,3) * t22;
t17 = mrSges(6,1) * t129 - mrSges(6,3) * t21;
t16 = t185 * t38 - t281;
t15 = -t181 * t38 - t276;
t8 = -mrSges(6,1) * t22 + mrSges(6,2) * t21;
t5 = -qJD(5) * t28 - t181 * t34 + t185 * t29;
t4 = qJD(5) * t27 + t181 * t29 + t185 * t34;
t1 = [(-Ifges(6,4) * t114 - Ifges(6,2) * t113) * t325 + (Ifges(4,4) * t355 + Ifges(4,2) * t356 - t248 / 0.2e1 - t247 / 0.2e1 - Ifges(6,6) * t325 - Ifges(6,5) * t326 - Ifges(6,3) * t312 + (-Ifges(4,2) * t183 + t288) * t354 - Ifges(5,3) * t311 - Ifges(5,6) * t317 - Ifges(5,5) * t318 + Ifges(4,6) * qJDD(3) + t350 + t351) * t187 + (Ifges(6,5) * t42 + Ifges(6,6) * t43) * t307 + (Ifges(6,4) * t42 + Ifges(6,2) * t43) * t315 + (Ifges(6,1) * t42 + Ifges(6,4) * t43) * t313 + (-t13 * t265 - t14 * t263 - t196 * t50 - t197 * t49) * mrSges(5,3) + (Ifges(4,1) * t149 + Ifges(4,5) * qJDD(3) + t212 * t311 + t215 * t317 + t218 * t318) * t183 + (-t11 * t42 - t113 * t2 + t114 * t3 + t12 * t43) * mrSges(6,3) + t202 * t354 + t288 * t355 + (t289 + t216) * t356 + (-t111 * t254 + t335) * mrSges(4,3) + (-Ifges(6,1) * t114 - Ifges(6,4) * t113) * t326 + (-Ifges(6,5) * t114 - Ifges(6,6) * t113) * t312 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t180 - 0.2e1 * mrSges(3,2) * t179 + m(3) * (t179 ^ 2 + t180 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (m(4) * t165 - t223) * t151 + (-t50 * mrSges(5,2) + t49 * mrSges(5,1) + t343 / 0.2e1 + Ifges(6,3) * t307 + Ifges(6,5) * t313 + Ifges(6,6) * t315 - t112 * mrSges(4,3) - t349) * t256 + t33 * (mrSges(6,1) * t113 - mrSges(6,2) * t114) - (t182 * t64 + t186 * t63) * t252 / 0.2e1 + (-t269 * t327 - m(3) * t176 - mrSges(2,1) * t188 - t109 * mrSges(5,1) - t99 * mrSges(6,1) + mrSges(2,2) * t184 - t108 * mrSges(5,2) - t98 * mrSges(6,2) + t334 * (t172 * pkin(2) + t171 * pkin(6) + t176) + t346 * t171 + (-m(5) * t225 - m(6) * t209 - t331) * t172) * g(2) + (-t268 * t327 + m(3) * t303 + mrSges(2,1) * t184 - t107 * mrSges(5,1) - t97 * mrSges(6,1) + mrSges(2,2) * t188 - t106 * mrSges(5,2) - t96 * mrSges(6,2) + t334 * (t172 * pkin(6) - t303) + t346 * t172 + (-m(5) * t207 - m(6) * (-pkin(2) - t209) + m(4) * pkin(2) + t331) * t171) * g(1) - t25 * t265 / 0.2e1 + t139 * (qJD(3) * t194 - t214 * t252) / 0.2e1 + t160 * (qJD(3) * t193 - t211 * t252) / 0.2e1 - t155 * t237 + t235 * t319 + t42 * t320 + t43 * t322 + t263 * t324 - t114 * t328 - t113 * t329 + (qJD(3) * t195 - t217 * t252) * t309 + m(6) * (t11 * t5 + t119 * t33 + t12 * t4 + t2 * t28 + t27 * t3 + t71 * t84) + t338 * t254 / 0.2e1 + t340 * t236 + t119 * t8 + t56 * t130 + t165 * (-mrSges(4,1) * t148 + mrSges(4,2) * t149) + qJD(3) * t203 + t101 * (mrSges(5,1) * t196 + mrSges(5,2) * t197) + m(5) * (t13 * t80 + t14 * t79 + t40 * t50 + t41 * t49) + (m(5) * (t101 * t254 + t183 * t56) + t341 * t183 + m(4) * ((-t111 * t187 - t112 * t183) * qJD(3) + t335) + t187 * t117) * t164 + qJD(3) ^ 2 * t213 / 0.2e1 + t27 * t17 + t28 * t18 + t4 * t53 + t5 * t54 + t71 * (-mrSges(6,1) * t43 + mrSges(6,2) * t42) + t79 * t46 + t80 * t47 + t84 * t37 + t40 * t94 + t41 * t95; -t113 * t17 - t114 * t18 + m(6) * (t11 * t43 - t113 * t3 - t114 * t2 + t12 * t42) + m(3) * qJDD(2) + t42 * t53 + t43 * t54 + (-m(3) + t334) * g(3) + (-t8 + (-t182 * t95 + t186 * t94 + t155) * qJD(3) + m(4) * (qJD(3) * t112 + t61) - m(6) * t33 + m(5) * (t255 * t50 - t257 * t49 - t56) - t341) * t187 + (t117 + (-t182 * t94 - t186 * t95) * qJD(4) + m(4) * t60 + m(5) * t352 + (-m(4) * t111 + m(6) * t71 - t330 + t37) * qJD(3) + t336) * t183; t352 * mrSges(5,3) + (g(1) * t172 + g(2) * t171) * (t358 * t183 + t357 * t187 + t222) + (t357 * t183 - t358 * t187 - t223) * g(3) + (t243 - t155) * t111 - (Ifges(6,4) * t313 + Ifges(6,2) * t315 + Ifges(6,6) * t307 + t296 + t322) * t82 - (Ifges(6,1) * t313 + Ifges(6,4) * t315 + Ifges(6,5) * t307 - t297 + t320) * t81 + (-pkin(3) * t56 - t49 * t65 - t50 * t66) * m(5) + (-Ifges(6,4) * t105 - Ifges(6,2) * t104) * t316 + (-Ifges(6,1) * t105 - Ifges(6,4) * t104) * t314 + (-Ifges(6,5) * t105 - Ifges(6,6) * t104) * t308 + (t12 * t104 - t11 * t105 - t3 * t142 - t2 * t208) * mrSges(6,3) + ((t105 - t81) * mrSges(6,2) + (-t104 + t82) * mrSges(6,1)) * t71 + (Ifges(6,4) * t142 - Ifges(6,2) * t208) * t325 + (Ifges(6,1) * t142 - Ifges(6,4) * t208) * t326 + (Ifges(6,5) * t142 - Ifges(6,6) * t208) * t312 + t33 * (mrSges(6,1) * t208 + mrSges(6,2) * t142) - t208 * t329 + (Ifges(6,5) * t314 + Ifges(6,6) * t316 + Ifges(6,3) * t308 + t349) * t260 + (t139 * t215 + t140 * t218 + t160 * t212) * qJD(4) / 0.2e1 + (t244 + t330) * t112 + (t139 * t194 + t140 * t195 + t160 * t193) * t347 + t64 * t251 / 0.2e1 - t213 * t250 / 0.2e1 + t63 * t238 / 0.2e1 - t343 * t260 / 0.2e1 + t344 * t53 + t345 * t54 + (t11 * t345 + t12 * t344 - t169 * t33 + t2 * t87 + t3 * t86 + t339 * t71) * m(6) + t253 * t319 - t105 * t321 - t104 * t323 + t182 * t324 + t142 * t328 + t211 * t311 + t214 * t317 + t217 * t318 + t160 * t101 * t220 + (-t95 * t251 - t94 * t253 + m(5) * ((-t182 * t50 - t186 * t49) * qJD(4) + t210) + t336) * pkin(7) - (-Ifges(4,2) * t260 + t170 + t338) * t259 / 0.2e1 + t339 * t37 + (-t49 * (mrSges(5,1) * t183 - mrSges(5,3) * t262) - t50 * (-mrSges(5,2) * t183 - mrSges(5,3) * t264) - t203 + t202 * t347) * qJD(1) + Ifges(4,6) * t148 + Ifges(4,5) * t149 - t169 * t8 + t186 * t25 / 0.2e1 + Ifges(4,3) * qJDD(3) + t56 * t221 - pkin(3) * t35 - t60 * mrSges(4,2) + t61 * mrSges(4,1) + t86 * t17 + t87 * t18 - t66 * t94 - t65 * t95; -t350 - (-Ifges(5,2) * t140 + t132 + t64) * t139 / 0.2e1 - (t71 * mrSges(6,1) + Ifges(6,4) * t314 + Ifges(6,2) * t316 + Ifges(6,6) * t308 - t296 + t323) * t78 + (-t71 * mrSges(6,2) + Ifges(6,1) * t314 + Ifges(6,4) * t316 + Ifges(6,5) * t308 + t297 + t321) * t226 - t37 * t302 - m(6) * (t11 * t15 + t12 * t16 + t302 * t71) - (-m(6) * t301 + t219) * t298 + ((-t181 * t54 + t185 * t53) * qJD(5) + t185 * t17 + t181 * t18) * pkin(4) - t140 * (Ifges(5,1) * t139 - t287) / 0.2e1 + t205 + (t181 * t2 + t185 * t3 + (-t11 * t181 + t12 * t185) * qJD(5)) * t327 + t63 * t309 + t247 + (-mrSges(5,2) * t107 + t106 * t342 - t293) * g(2) + (mrSges(5,2) * t109 - t108 * t342 - t292) * g(1) + (t290 + t95) * t50 + (t291 - t94) * t49 + g(3) * t130 - t101 * (mrSges(5,1) * t140 + mrSges(5,2) * t139) - t160 * (Ifges(5,5) * t139 - Ifges(5,6) * t140) / 0.2e1 - t16 * t53 - t15 * t54; -t71 * (mrSges(6,1) * t78 + mrSges(6,2) * t226) + (Ifges(6,1) * t226 - t305) * t314 + t31 * t313 + (Ifges(6,5) * t226 - Ifges(6,6) * t78) * t308 - t11 * t53 + t12 * t54 - g(1) * t292 - g(2) * t293 - t219 * t298 + (t11 * t226 + t12 * t78) * mrSges(6,3) + t205 + (-Ifges(6,2) * t78 + t32 + t74) * t316;];
tau = t1;
