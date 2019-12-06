% Calculate vector of inverse dynamics joint torques for
% S5PRRRR7
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:48
% EndTime: 2019-12-05 17:12:06
% DurationCPUTime: 8.03s
% Computational Cost: add. (4653->462), mult. (10361->644), div. (0->0), fcn. (7251->14), ass. (0->220)
t195 = sin(qJ(3));
t201 = -pkin(7) - pkin(6);
t236 = qJD(3) * t201;
t155 = t195 * t236;
t199 = cos(qJ(3));
t156 = t199 * t236;
t165 = t201 * t195;
t166 = t201 * t199;
t194 = sin(qJ(4));
t198 = cos(qJ(4));
t221 = t194 * t195 - t198 * t199;
t245 = qJD(4) * t198;
t246 = qJD(4) * t194;
t200 = cos(qJ(2));
t252 = qJD(1) * t200;
t308 = t155 * t198 + t156 * t194 + t165 * t245 + t166 * t246 + t221 * t252;
t105 = t165 * t194 - t166 * t198;
t149 = t194 * t199 + t195 * t198;
t213 = t149 * t200;
t307 = qJD(1) * t213 - qJD(4) * t105 - t155 * t194 + t156 * t198;
t211 = t221 * qJD(4);
t99 = -qJD(3) * t221 - t211;
t323 = -pkin(8) * t99 + t307;
t212 = t149 * qJD(4);
t100 = -qJD(3) * t149 - t212;
t322 = pkin(8) * t100 + t308;
t190 = qJ(3) + qJ(4);
t183 = cos(t190);
t287 = pkin(3) * t199;
t162 = pkin(4) * t183 + t287;
t184 = qJ(5) + t190;
t173 = sin(t184);
t174 = cos(t184);
t176 = pkin(2) + t287;
t182 = sin(t190);
t228 = -mrSges(4,1) * t199 + mrSges(4,2) * t195;
t321 = mrSges(3,1) + m(6) * (pkin(2) + t162) + t174 * mrSges(6,1) - t173 * mrSges(6,2) + m(5) * t176 + t183 * mrSges(5,1) - t182 * mrSges(5,2) + m(4) * pkin(2) - t228;
t320 = mrSges(3,2) + m(6) * (-pkin(8) + t201) - mrSges(6,3) + m(5) * t201 - mrSges(5,3) - m(4) * pkin(6) - mrSges(4,3);
t185 = qJDD(3) + qJDD(4);
t179 = qJDD(5) + t185;
t193 = sin(qJ(5));
t197 = cos(qJ(5));
t139 = t221 * qJD(2);
t140 = t149 * qJD(2);
t230 = -t139 * t197 - t140 * t193;
t240 = qJD(2) * qJD(3);
t157 = qJDD(2) * t199 - t195 * t240;
t158 = qJDD(2) * t195 + t199 * t240;
t68 = -qJD(2) * t211 + t157 * t194 + t158 * t198;
t69 = -qJD(2) * t212 + t157 * t198 - t158 * t194;
t19 = qJD(5) * t230 + t193 * t69 + t197 * t68;
t284 = pkin(8) * t139;
t196 = sin(qJ(2));
t242 = t196 * qJD(1);
t167 = qJD(2) * pkin(6) + t242;
t231 = pkin(7) * qJD(2) + t167;
t128 = t231 * t199;
t120 = t198 * t128;
t127 = t231 * t195;
t121 = qJD(3) * pkin(3) - t127;
t73 = t121 * t194 + t120;
t47 = t73 - t284;
t273 = t193 * t47;
t186 = qJD(3) + qJD(4);
t134 = t140 * pkin(8);
t118 = t194 * t128;
t72 = t121 * t198 - t118;
t46 = -t134 + t72;
t44 = pkin(4) * t186 + t46;
t21 = t197 * t44 - t273;
t241 = qJD(1) * qJD(2);
t171 = t200 * t241;
t160 = qJDD(1) * t196 + t171;
t145 = qJDD(2) * pkin(6) + t160;
t247 = qJD(3) * t199;
t97 = -t145 * t195 - t167 * t247;
t67 = qJDD(3) * pkin(3) - pkin(7) * t158 + t97;
t248 = qJD(3) * t195;
t96 = t145 * t199 - t167 * t248;
t74 = pkin(7) * t157 + t96;
t16 = -qJD(4) * t73 - t194 * t74 + t198 * t67;
t7 = pkin(4) * t185 - pkin(8) * t68 + t16;
t15 = t121 * t245 - t128 * t246 + t194 * t67 + t198 * t74;
t8 = pkin(8) * t69 + t15;
t2 = qJD(5) * t21 + t193 * t7 + t197 * t8;
t88 = -t139 * t193 + t140 * t197;
t20 = -qJD(5) * t88 - t193 * t68 + t197 * t69;
t271 = t197 * t47;
t22 = t193 * t44 + t271;
t3 = -qJD(5) * t22 - t193 * t8 + t197 * t7;
t180 = qJD(5) + t186;
t82 = Ifges(6,4) * t230;
t38 = Ifges(6,1) * t88 + Ifges(6,5) * t180 + t82;
t319 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t179 - (-Ifges(6,2) * t88 + t38 + t82) * t230 / 0.2e1;
t290 = Ifges(6,4) * t88;
t37 = Ifges(6,2) * t230 + Ifges(6,6) * t180 + t290;
t318 = t37 / 0.2e1;
t317 = t157 / 0.2e1;
t104 = t165 * t198 + t166 * t194;
t83 = -pkin(8) * t149 + t104;
t84 = -pkin(8) * t221 + t105;
t39 = -t193 * t84 + t197 * t83;
t316 = qJD(5) * t39 + t193 * t323 + t197 * t322;
t40 = t193 * t83 + t197 * t84;
t315 = -qJD(5) * t40 - t193 * t322 + t197 * t323;
t288 = pkin(3) * t198;
t175 = pkin(4) + t288;
t243 = qJD(5) * t197;
t244 = qJD(5) * t193;
t263 = t194 * t197;
t75 = t127 * t194 - t120;
t50 = t75 + t284;
t76 = -t127 * t198 - t118;
t51 = -t134 + t76;
t314 = t193 * t51 - t197 * t50 - t175 * t244 + (-t194 * t243 + (-t193 * t198 - t263) * qJD(4)) * pkin(3);
t264 = t193 * t194;
t313 = -t193 * t50 - t197 * t51 + t175 * t243 + (-t194 * t244 + (t197 * t198 - t264) * qJD(4)) * pkin(3);
t312 = mrSges(6,2) * t230;
t311 = Ifges(6,1) * t230;
t310 = Ifges(6,5) * t230;
t298 = m(5) * pkin(3);
t309 = -mrSges(4,1) - t298;
t130 = t221 * t196;
t225 = -mrSges(6,1) * t173 - mrSges(6,2) * t174;
t306 = mrSges(5,1) * t182 + mrSges(5,2) * t183 - t225;
t191 = sin(pkin(9));
t192 = cos(pkin(9));
t265 = t192 * t200;
t219 = -t182 * t265 + t183 * t191;
t255 = (-t173 * t265 + t174 * t191) * mrSges(6,1) + (-t173 * t191 - t174 * t265) * mrSges(6,2);
t305 = -t219 * mrSges(5,1) - (-t182 * t191 - t183 * t265) * mrSges(5,2) - t255;
t266 = t191 * t200;
t218 = -t182 * t266 - t183 * t192;
t256 = (-t173 * t266 - t174 * t192) * mrSges(6,1) + (t173 * t192 - t174 * t266) * mrSges(6,2);
t304 = -t218 * mrSges(5,1) - (t182 * t192 - t183 * t266) * mrSges(5,2) - t256;
t251 = qJD(2) * t195;
t163 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t251;
t249 = qJD(2) * t199;
t164 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t249;
t303 = t163 * t199 + t164 * t195;
t302 = t199 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t157) - t195 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t158);
t301 = t163 * t195 - t164 * t199;
t222 = -t195 * t97 + t199 * t96;
t42 = -mrSges(6,1) * t230 + mrSges(6,2) * t88;
t93 = mrSges(5,1) * t139 + mrSges(5,2) * t140;
t300 = qJD(2) * t228 + t42 + t93;
t202 = qJD(2) ^ 2;
t297 = m(6) * pkin(4);
t295 = -t88 / 0.2e1;
t294 = t88 / 0.2e1;
t292 = t140 / 0.2e1;
t291 = -t180 / 0.2e1;
t289 = pkin(3) * t195;
t286 = pkin(4) * t140;
t285 = pkin(4) * t182;
t281 = g(3) * t196;
t280 = t21 * mrSges(6,3);
t279 = t22 * mrSges(6,3);
t278 = mrSges(5,3) * t139;
t277 = mrSges(5,3) * t140;
t276 = Ifges(4,4) * t195;
t275 = Ifges(4,4) * t199;
t274 = t140 * Ifges(5,4);
t260 = t195 * t200;
t257 = t199 * t200;
t250 = qJD(2) * t196;
t239 = pkin(3) * t251;
t238 = pkin(3) * t248;
t170 = t196 * t241;
t159 = qJDD(1) * t200 - t170;
t227 = mrSges(4,1) * t195 + mrSges(4,2) * t199;
t224 = t199 * Ifges(4,2) + t276;
t223 = Ifges(4,5) * t199 - Ifges(4,6) * t195;
t129 = t149 * t196;
t77 = -t129 * t197 + t130 * t193;
t78 = -t129 * t193 - t130 * t197;
t94 = -t149 * t193 - t197 * t221;
t95 = t149 * t197 - t193 * t221;
t168 = -qJD(2) * pkin(2) - t252;
t215 = t168 * t227;
t214 = t195 * (Ifges(4,1) * t199 - t276);
t144 = -qJDD(2) * pkin(2) - t159;
t207 = t168 * t196 + (t195 ^ 2 + t199 ^ 2) * t200 * t167;
t146 = -qJD(2) * t176 - t252;
t101 = -pkin(3) * t157 + t144;
t133 = Ifges(5,4) * t139;
t80 = -t139 * Ifges(5,2) + t186 * Ifges(5,6) + t274;
t81 = t140 * Ifges(5,1) + t186 * Ifges(5,5) - t133;
t98 = pkin(4) * t139 + t146;
t203 = t16 * mrSges(5,1) - t15 * mrSges(5,2) - t146 * (mrSges(5,1) * t140 - mrSges(5,2) * t139) + t277 * t73 + t230 * t280 - t98 * t312 + Ifges(5,3) * t185 + t311 * t295 + t310 * t291 - t72 * t278 + t80 * t292 - t140 * (-Ifges(5,1) * t139 - t274) / 0.2e1 + Ifges(5,6) * t69 + Ifges(5,5) * t68 - t186 * (-Ifges(5,5) * t139 - Ifges(5,6) * t140) / 0.2e1 + (-t98 * mrSges(6,1) - Ifges(6,4) * t295 - Ifges(6,6) * t291 + t279 + t318) * t88 + (-Ifges(5,2) * t140 - t133 + t81) * t139 / 0.2e1 + t319;
t177 = Ifges(4,4) * t249;
t161 = -t285 - t289;
t138 = Ifges(4,1) * t251 + Ifges(4,5) * qJD(3) + t177;
t137 = Ifges(4,6) * qJD(3) + qJD(2) * t224;
t136 = pkin(3) * t263 + t175 * t193;
t135 = -pkin(3) * t264 + t175 * t197;
t124 = pkin(4) * t221 - t176;
t113 = mrSges(5,1) * t186 - t277;
t112 = -mrSges(5,2) * t186 - t278;
t111 = t239 + t286;
t102 = -mrSges(4,1) * t157 + mrSges(4,2) * t158;
t89 = -pkin(4) * t100 + t238;
t71 = mrSges(6,1) * t180 - mrSges(6,3) * t88;
t70 = -mrSges(6,2) * t180 + mrSges(6,3) * t230;
t57 = -mrSges(5,2) * t185 + mrSges(5,3) * t69;
t56 = mrSges(5,1) * t185 - mrSges(5,3) * t68;
t49 = -qJD(2) * t213 + t130 * t186;
t48 = t100 * t196 - t139 * t200;
t43 = -pkin(4) * t69 + t101;
t33 = -mrSges(5,1) * t69 + mrSges(5,2) * t68;
t32 = -qJD(5) * t95 + t100 * t197 - t193 * t99;
t31 = qJD(5) * t94 + t100 * t193 + t197 * t99;
t24 = t197 * t46 - t273;
t23 = -t193 * t46 - t271;
t13 = -mrSges(6,2) * t179 + mrSges(6,3) * t20;
t12 = mrSges(6,1) * t179 - mrSges(6,3) * t19;
t11 = -qJD(5) * t78 - t193 * t48 + t197 * t49;
t10 = qJD(5) * t77 + t193 * t49 + t197 * t48;
t4 = -mrSges(6,1) * t20 + mrSges(6,2) * t19;
t1 = [m(2) * qJDD(1) + t10 * t70 + t11 * t71 + t48 * t112 + t49 * t113 + t77 * t12 - t129 * t56 + t78 * t13 - t130 * t57 + (-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3) + (qJDD(2) * mrSges(3,1) - t202 * mrSges(3,2) - qJD(2) * t301 - t102 - t33 - t4) * t200 + (-t202 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + qJD(2) * t300 - qJD(3) * t303 + t302) * t196 + m(6) * (t10 * t22 + t11 * t21 + t2 * t78 - t200 * t43 + t250 * t98 + t3 * t77) + m(5) * (-t101 * t200 - t129 * t16 - t130 * t15 + t146 * t250 + t48 * t73 + t49 * t72) + m(3) * (t159 * t200 + t160 * t196) + m(4) * (qJD(2) * t207 - t144 * t200 + t196 * t222); t158 * t275 / 0.2e1 - (-mrSges(5,1) * t101 + mrSges(5,3) * t15 + Ifges(5,4) * t68 + Ifges(5,2) * t69 + Ifges(5,6) * t185) * t221 + (t199 * (-Ifges(4,2) * t195 + t275) + t214) * t240 / 0.2e1 + t230 * (Ifges(6,4) * t31 + Ifges(6,2) * t32) / 0.2e1 + (t170 + t159) * mrSges(3,1) + (g(1) * t192 + g(2) * t191) * (t196 * t321 + t200 * t320) + (t196 * t320 - t200 * t321) * g(3) + (-pkin(2) * t144 - qJD(1) * t207) * m(4) + (t100 * t73 - t72 * t99) * mrSges(5,3) + t144 * t228 + (t215 + t223 * qJD(3) / 0.2e1) * qJD(3) + t40 * t13 + t31 * t38 / 0.2e1 + t39 * t12 + t93 * t238 + (mrSges(6,2) * t43 - mrSges(6,3) * t3 + Ifges(6,1) * t19 + Ifges(6,4) * t20 + Ifges(6,5) * t179) * t95 + t199 * (Ifges(4,4) * t158 + Ifges(4,2) * t157) / 0.2e1 + (Ifges(5,1) * t99 + Ifges(5,4) * t100) * t292 + (Ifges(6,1) * t31 + Ifges(6,4) * t32) * t294 + t32 * t279 + t32 * t318 + t224 * t317 + (t171 - t160) * mrSges(3,2) + (-mrSges(6,1) * t43 + mrSges(6,3) * t2 + Ifges(6,4) * t19 + Ifges(6,2) * t20 + Ifges(6,6) * t179) * t94 + t89 * t42 + t98 * (-mrSges(6,1) * t32 + mrSges(6,2) * t31) + t99 * t81 / 0.2e1 + t100 * t80 / 0.2e1 - pkin(2) * t102 + t104 * t56 + t105 * t57 + (-m(5) * t146 - m(6) * t98 - t300) * t242 + t222 * mrSges(4,3) + t301 * t252 + (m(4) * t222 - t163 * t247 - t164 * t248 + t302) * pkin(6) + Ifges(3,3) * qJDD(2) + t124 * t4 - t139 * (Ifges(5,4) * t99 + Ifges(5,2) * t100) / 0.2e1 + t307 * t113 + t308 * t112 + (-t101 * t176 + t104 * t16 + t105 * t15 + t146 * t238 + t307 * t72 + t308 * t73) * m(5) + t146 * (-mrSges(5,1) * t100 + mrSges(5,2) * t99) + t315 * t71 + t316 * t70 + (t124 * t43 + t2 * t40 + t21 * t315 + t22 * t316 + t3 * t39 + t89 * t98) * m(6) + (Ifges(4,1) * t158 + Ifges(4,4) * t317) * t195 - t176 * t33 + t180 * (Ifges(6,5) * t31 + Ifges(6,6) * t32) / 0.2e1 + t186 * (Ifges(5,5) * t99 + Ifges(5,6) * t100) / 0.2e1 + qJDD(3) * (Ifges(4,5) * t195 + Ifges(4,6) * t199) + t138 * t247 / 0.2e1 - t137 * t248 / 0.2e1 - t31 * t280 + (mrSges(5,2) * t101 - mrSges(5,3) * t16 + Ifges(5,1) * t68 + Ifges(5,4) * t69 + Ifges(5,5) * t185) * t149; -(-Ifges(4,2) * t251 + t138 + t177) * t249 / 0.2e1 - qJD(2) * t215 + t203 + (t15 * t194 + t16 * t198 + (-t194 * t72 + t198 * t73) * qJD(4)) * t298 + t56 * t288 + (t112 * t245 - t113 * t246 + t194 * t57) * pkin(3) - t96 * mrSges(4,2) + t97 * mrSges(4,1) - t111 * t42 - t76 * t112 - t75 * t113 + t303 * t167 + (m(5) * t289 - m(6) * t161 + t227 + t306) * t281 + t135 * t12 + t136 * t13 + (-(-t191 * t257 + t192 * t195) * mrSges(4,2) - m(6) * (t161 * t266 - t162 * t192) + t309 * (-t191 * t260 - t192 * t199) + t304) * g(2) + (-(-t191 * t195 - t192 * t257) * mrSges(4,2) - m(6) * (t161 * t265 + t162 * t191) + t309 * (t191 * t199 - t192 * t260) + t305) * g(1) + Ifges(4,6) * t157 + Ifges(4,5) * t158 + t313 * t70 + t314 * t71 + (-t111 * t98 + t135 * t3 + t136 * t2 + t314 * t21 + t313 * t22) * m(6) - t93 * t239 - m(5) * (t146 * t239 + t72 * t75 + t73 * t76) - t223 * t240 / 0.2e1 + t137 * t251 / 0.2e1 - t202 * t214 / 0.2e1 + Ifges(4,3) * qJDD(3); t203 - m(6) * (t21 * t23 + t22 * t24 + t286 * t98) + (t193 * t2 + t197 * t3 + (-t193 * t21 + t197 * t22) * qJD(5)) * t297 - t24 * t70 - t23 * t71 - t72 * t112 + t73 * t113 - t42 * t286 + (m(6) * t285 + t306) * t281 + (-t218 * t297 + t304) * g(2) + (-t219 * t297 + t305) * g(1) + (t197 * t12 + t193 * t13 + t243 * t70 - t244 * t71) * pkin(4); -t98 * (mrSges(6,1) * t88 + t312) + (-t290 + t311) * t295 + t37 * t294 + (-Ifges(6,6) * t88 + t310) * t291 - t21 * t70 + t22 * t71 - g(1) * t255 - g(2) * t256 - t225 * t281 + (t21 * t230 + t22 * t88) * mrSges(6,3) + t319;];
tau = t1;
