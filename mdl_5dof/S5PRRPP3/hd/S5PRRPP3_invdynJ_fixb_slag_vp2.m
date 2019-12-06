% Calculate vector of inverse dynamics joint torques for
% S5PRRPP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:00
% EndTime: 2019-12-05 16:12:24
% DurationCPUTime: 9.40s
% Computational Cost: add. (1828->420), mult. (4190->556), div. (0->0), fcn. (2635->8), ass. (0->193)
t302 = Ifges(6,1) + Ifges(5,1);
t140 = sin(qJ(3));
t142 = cos(qJ(3));
t206 = qJD(2) * qJD(3);
t101 = -t142 * qJDD(2) + t140 * t206;
t254 = t101 / 0.2e1;
t102 = qJDD(2) * t140 + t142 * t206;
t136 = sin(pkin(8));
t138 = cos(pkin(8));
t66 = qJDD(3) * t136 + t102 * t138;
t255 = t66 / 0.2e1;
t284 = Ifges(5,5) + Ifges(6,4);
t304 = t284 * t254 + t302 * t255;
t65 = -t138 * qJDD(3) + t102 * t136;
t256 = t65 / 0.2e1;
t291 = -m(5) - m(6);
t303 = mrSges(3,2) - mrSges(4,3);
t285 = -mrSges(5,3) - mrSges(6,2);
t301 = -Ifges(5,4) + Ifges(6,5);
t283 = -Ifges(6,6) + Ifges(5,6);
t282 = -Ifges(5,3) - Ifges(6,2);
t183 = mrSges(4,1) * t142 - mrSges(4,2) * t140;
t300 = -mrSges(3,1) - t183;
t239 = Ifges(6,5) * t136;
t241 = Ifges(5,4) * t136;
t299 = t302 * t138 + t239 - t241;
t288 = -t101 / 0.2e1;
t169 = pkin(3) * t140 - qJ(4) * t142;
t209 = qJD(4) * t140;
t79 = qJD(3) * t169 - t209;
t141 = sin(qJ(2));
t143 = cos(qJ(2));
t220 = t142 * t143;
t88 = t136 * t141 + t138 * t220;
t298 = -t88 * qJD(1) + t136 * t79;
t217 = qJD(1) * t143;
t200 = t142 * t217;
t218 = qJD(1) * t141;
t297 = t136 * t200 + (-t218 + t79) * t138;
t215 = qJD(2) * t142;
t127 = Ifges(4,4) * t215;
t216 = qJD(2) * t140;
t96 = -t138 * qJD(3) + t136 * t216;
t97 = qJD(3) * t136 + t138 * t216;
t296 = Ifges(4,1) * t216 + Ifges(4,5) * qJD(3) + t127 + t138 * (-t215 * t284 + t301 * t96 + t302 * t97);
t179 = -t138 * mrSges(6,1) - t136 * mrSges(6,3);
t181 = -mrSges(5,1) * t138 + mrSges(5,2) * t136;
t295 = -t179 - t181;
t207 = qJD(1) * qJD(2);
t125 = t141 * t207;
t103 = qJDD(1) * t143 - t125;
t93 = -qJDD(2) * pkin(2) - t103;
t13 = pkin(3) * t101 - qJ(4) * t102 - qJD(2) * t209 + t93;
t112 = qJD(2) * pkin(6) + t218;
t105 = t140 * t112;
t126 = t143 * t207;
t104 = t141 * qJDD(1) + t126;
t94 = qJDD(2) * pkin(6) + t104;
t74 = t142 * t94;
t21 = qJDD(3) * qJ(4) + t74 + (qJD(4) - t105) * qJD(3);
t4 = t13 * t138 - t136 * t21;
t2 = -pkin(4) * t101 + qJDD(5) - t4;
t294 = t2 * mrSges(6,2) + t256 * t301 + t304;
t243 = Ifges(4,4) * t140;
t273 = Ifges(4,2) * t142;
t175 = t243 + t273;
t293 = t284 * t97 / 0.2e1 - t283 * t96 / 0.2e1 + t282 * t215 / 0.2e1 - Ifges(4,6) * qJD(3) / 0.2e1 - qJD(2) * t175 / 0.2e1;
t287 = t102 / 0.2e1;
t286 = qJD(3) / 0.2e1;
t16 = t65 * mrSges(5,1) + t66 * mrSges(5,2);
t281 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t102 + t16;
t23 = -mrSges(6,2) * t65 + mrSges(6,3) * t101;
t24 = -mrSges(5,2) * t101 - mrSges(5,3) * t65;
t280 = t23 + t24;
t25 = mrSges(5,1) * t101 - mrSges(5,3) * t66;
t26 = -t101 * mrSges(6,1) + t66 * mrSges(6,2);
t279 = t26 - t25;
t202 = pkin(6) * t136 + pkin(4);
t213 = qJD(3) * t140;
t277 = -t202 * t213 - t297;
t61 = -mrSges(6,2) * t96 - mrSges(6,3) * t215;
t62 = mrSges(5,2) * t215 - mrSges(5,3) * t96;
t245 = t61 + t62;
t63 = -mrSges(5,1) * t215 - mrSges(5,3) * t97;
t64 = mrSges(6,1) * t215 + mrSges(6,2) * t97;
t244 = t64 - t63;
t204 = pkin(6) * t213;
t276 = t136 * t204 + t297;
t208 = qJD(5) * t142;
t275 = -t208 + (-pkin(6) * t138 + qJ(5)) * t213 + t298;
t274 = -t138 * t204 + t298;
t272 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t96 + mrSges(5,2) * t97 + mrSges(4,3) * t216;
t36 = -t112 * t213 + t74;
t212 = qJD(3) * t142;
t37 = -t112 * t212 - t140 * t94;
t166 = -t140 * t37 + t142 * t36;
t5 = t136 * t13 + t138 * t21;
t271 = -t136 * t4 + t138 * t5;
t270 = -m(4) + t291;
t269 = mrSges(4,2) + t285;
t41 = mrSges(6,1) * t96 - mrSges(6,3) * t97;
t205 = t41 + t272;
t168 = -pkin(4) * t138 - qJ(5) * t136;
t268 = -m(6) * t168 + mrSges(4,1) + t295;
t238 = Ifges(6,5) * t138;
t170 = Ifges(6,3) * t136 + t238;
t267 = t96 * (Ifges(6,6) * t140 + t142 * t170) + (t140 * t284 + t142 * t299) * t97;
t264 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t262 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t1 = qJ(5) * t101 - qJD(2) * t208 + t5;
t261 = Ifges(6,5) * t255 + Ifges(6,6) * t254 - t66 * Ifges(5,4) / 0.2e1 + Ifges(5,6) * t288 - t1 * mrSges(6,2) + (Ifges(6,3) + Ifges(5,2)) * t256;
t108 = -pkin(3) * t142 - qJ(4) * t140 - pkin(2);
t68 = qJD(2) * t108 - t217;
t106 = t142 * t112;
t90 = qJD(3) * qJ(4) + t106;
t17 = -t136 * t90 + t138 * t68;
t11 = pkin(4) * t215 + qJD(5) - t17;
t113 = -qJD(2) * pkin(2) - t217;
t18 = t136 * t68 + t138 * t90;
t12 = -qJ(5) * t215 + t18;
t80 = -qJD(3) * pkin(3) + qJD(4) + t105;
t14 = pkin(4) * t96 - qJ(5) * t97 + t80;
t240 = Ifges(5,4) * t138;
t174 = -Ifges(5,2) * t136 + t240;
t178 = mrSges(6,1) * t136 - mrSges(6,3) * t138;
t180 = mrSges(5,1) * t136 + mrSges(5,2) * t138;
t182 = mrSges(4,1) * t140 + mrSges(4,2) * t142;
t225 = t138 * t142;
t227 = t136 * t142;
t258 = -t18 * (-mrSges(5,2) * t140 - mrSges(5,3) * t227) - t17 * (mrSges(5,1) * t140 - mrSges(5,3) * t225) - t12 * (-mrSges(6,2) * t227 + mrSges(6,3) * t140) - t11 * (-mrSges(6,1) * t140 + mrSges(6,2) * t225) - t142 * t80 * t180 - t14 * t142 * t178 + t96 * (Ifges(5,6) * t140 + t142 * t174) / 0.2e1 - t113 * t182;
t144 = qJD(2) ^ 2;
t257 = -t65 / 0.2e1;
t242 = Ifges(4,4) * t142;
t34 = -qJDD(3) * pkin(3) + qJDD(4) - t37;
t236 = t140 * t34;
t233 = t142 * (-qJDD(3) * mrSges(4,2) - mrSges(4,3) * t101);
t60 = pkin(6) * t225 + t136 * t108;
t100 = t169 * qJD(2);
t230 = t100 * t138;
t229 = t108 * t138;
t226 = t138 * t140;
t224 = t140 * t141;
t223 = t140 * t143;
t222 = t141 * t138;
t221 = t141 * t142;
t219 = t143 * pkin(2) + t141 * pkin(6);
t214 = qJD(2) * t143;
t199 = t136 * t215;
t198 = t141 * t213;
t197 = t141 * t212;
t15 = t65 * mrSges(6,1) - t66 * mrSges(6,3);
t194 = -t215 / 0.2e1;
t189 = -t206 / 0.2e1;
t188 = t206 / 0.2e1;
t173 = Ifges(6,4) * t138 + Ifges(6,6) * t136;
t172 = Ifges(4,5) * t142 - Ifges(4,6) * t140;
t171 = Ifges(5,5) * t138 - Ifges(5,6) * t136;
t167 = pkin(4) * t136 - qJ(5) * t138;
t107 = -pkin(3) + t168;
t165 = pkin(6) + t167;
t84 = -t136 * t143 + t138 * t221;
t83 = t136 * t221 + t138 * t143;
t163 = t140 * (Ifges(4,1) * t142 - t243);
t137 = sin(pkin(7));
t139 = cos(pkin(7));
t81 = t137 * t223 + t139 * t142;
t85 = -t137 * t142 + t139 * t223;
t156 = -g(1) * t85 - g(2) * t81 - g(3) * t224;
t155 = t113 * t141 + (t140 ^ 2 + t142 ^ 2) * t143 * t112;
t150 = t142 * (Ifges(6,2) * t140 + t142 * t173);
t149 = t142 * (Ifges(5,3) * t140 + t142 * t171);
t111 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t215;
t99 = t183 * qJD(2);
t89 = t136 * t100;
t86 = t137 * t140 + t139 * t220;
t82 = t137 * t220 - t139 * t140;
t69 = t165 * t140;
t59 = -pkin(6) * t227 + t229;
t50 = t142 * t202 - t229;
t49 = -qJ(5) * t142 + t60;
t48 = mrSges(4,1) * t101 + mrSges(4,2) * t102;
t47 = qJD(2) * t88 - t138 * t198;
t46 = -qJD(2) * t222 - t136 * t198 + t143 * t199;
t43 = t167 * t215 + t106;
t40 = -t105 * t138 + t89;
t39 = t105 * t136 + t230;
t38 = -qJD(5) * t226 + t165 * t212;
t31 = t97 * Ifges(5,4) - t96 * Ifges(5,2) - Ifges(5,6) * t215;
t28 = t97 * Ifges(6,5) - Ifges(6,6) * t215 + t96 * Ifges(6,3);
t27 = -t230 + (-pkin(4) * qJD(2) - t112 * t136) * t140;
t22 = t89 + (qJ(5) * qJD(2) - t112 * t138) * t140;
t3 = pkin(4) * t65 - qJ(5) * t66 - qJD(5) * t97 + t34;
t6 = [m(2) * qJDD(1) + t280 * t84 + t279 * t83 + t245 * t47 + t244 * t46 + (-m(2) - m(3) + t270) * g(3) + (qJDD(2) * mrSges(3,1) - t144 * mrSges(3,2) - t48 + (t111 * t142 + t140 * t205) * qJD(2)) * t143 + (-t144 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - qJD(2) * t99 + t233 + (t15 + t281) * t140 + (-t111 * t140 + t142 * t205) * qJD(3)) * t141 + m(4) * (qJD(2) * t155 + t141 * t166 - t143 * t93) + m(6) * (t3 * t224 + t1 * t84 + t11 * t46 + t12 * t47 + t2 * t83 + (t140 * t214 + t197) * t14) + m(5) * (t80 * t197 - t17 * t46 + t18 * t47 - t4 * t83 + t5 * t84 + (t141 * t34 + t214 * t80) * t140) + m(3) * (t103 * t143 + t104 * t141); t163 * t188 + (g(1) * t139 + g(2) * t137) * (t83 * t262 - t84 * t264 + (m(4) * pkin(2) + t108 * t291 - t140 * t285 - t300) * t141 + (pkin(6) * t270 + t303) * t143) + qJDD(3) * Ifges(4,6) * t142 + t166 * mrSges(4,3) - t93 * t183 + (-t4 * mrSges(5,1) + t2 * mrSges(6,1) + t5 * mrSges(5,2) - t1 * mrSges(6,3) + Ifges(4,4) * t287 + Ifges(4,2) * t288 - Ifges(5,6) * t257 - Ifges(6,6) * t256 + t242 * t188 + t254 * t282 - t255 * t284) * t142 + (t172 * t286 - t258) * qJD(3) + (t1 * t49 + t11 * t277 + t12 * t275 + t14 * t38 + t2 * t50 + t3 * t69) * m(6) + t277 * t64 - (-t101 * t282 - t283 * t65 + t284 * t66) * t142 / 0.2e1 + t180 * t236 + (t150 + t149) * t189 + t274 * t62 + t275 * t61 + t276 * t63 + (t281 * pkin(6) + t170 * t256 + t174 * t257 + Ifges(4,1) * t102 + Ifges(4,4) * t288 - t188 * t273 + t3 * t178 + t299 * t255 + (t171 + t173) * t254 + Ifges(4,5) * qJDD(3) + (-m(5) * t80 - m(6) * t14 - t205) * t217 + (-t5 * mrSges(5,3) + t261) * t136) * t140 + (t276 * t17 + t274 * t18 + t4 * t59 + t5 * t60) * m(5) + (-pkin(2) * t93 - qJD(1) * t155) * m(4) + t69 * t15 + t59 * t25 + t60 * t24 + t38 * t41 - pkin(2) * t48 + t49 * t23 + t50 * t26 + (-m(4) * t219 + t291 * (pkin(3) * t220 + qJ(4) * t223 + t219) + t264 * t88 - t262 * (t136 * t220 - t222) + t285 * t223 + t300 * t143 + t303 * t141) * g(3) + (t136 * t28 + t296) * t212 / 0.2e1 - t136 * t31 * t212 / 0.2e1 + (t233 + t272 * t212 + m(5) * (t212 * t80 + t236) + m(4) * t166) * pkin(6) + t293 * t213 + (-t4 * mrSges(5,3) + t294) * t226 + t267 * t286 + t99 * t218 + t242 * t287 + t175 * t288 + (t103 + t125) * mrSges(3,1) + (-t104 + t126) * mrSges(3,2) + (-t204 - t200) * t111 + Ifges(3,3) * qJDD(2); t172 * t189 + (t28 * t194 + m(6) * (qJ(4) * t2 - qJD(5) * t14) - qJD(5) * t41 + t279 * qJ(4) + (-t17 * m(5) + t11 * m(6) + t244) * qJD(4) + t294 + t304) * t136 + t271 * mrSges(5,3) + t34 * t181 + t3 * t179 + (t291 * (-t85 * pkin(3) + qJ(4) * t86) + t269 * t86 + t268 * t85) * g(1) + (t291 * (-t81 * pkin(3) + qJ(4) * t82) + t269 * t82 + t268 * t81) * g(2) + (t149 / 0.2e1 - t163 / 0.2e1 + t150 / 0.2e1) * t144 + t111 * t105 - Ifges(4,6) * t101 + Ifges(4,5) * t102 + t107 * t15 - t272 * t106 + (t291 * qJ(4) * t221 + (t182 + t285 * t142 + (m(5) * pkin(3) - m(6) * t107 + t295) * t140) * t141) * g(3) + (t258 - t267 / 0.2e1) * qJD(2) - t27 * t64 - t22 * t61 - t40 * t62 - t39 * t63 - t43 * t41 - t36 * mrSges(4,2) + t37 * mrSges(4,1) - pkin(3) * t16 + (-Ifges(4,2) * t216 + t127 + t296) * t194 + (-t238 + t240) * t255 + (Ifges(5,2) * t257 - Ifges(6,3) * t256 + t283 * t254 - t261 + (t1 * m(6) + t280) * qJ(4) + (t18 * m(5) + t12 * m(6) + t245) * qJD(4)) * t138 + t239 * t256 + t241 * t257 + t31 * t199 / 0.2e1 + (t107 * t3 - t11 * t27 - t12 * t22 - t14 * t43) * m(6) + (-pkin(3) * t34 + t271 * qJ(4) - t80 * t106 - t17 * t39 - t18 * t40) * m(5) - t293 * t216 + Ifges(4,3) * qJDD(3); -t244 * t97 + t245 * t96 + t15 + t16 + (-t11 * t97 + t12 * t96 + t156 + t3) * m(6) + (t17 * t97 + t18 * t96 + t156 + t34) * m(5); t61 * t215 + t97 * t41 + (-g(1) * (t136 * t86 - t139 * t222) - g(2) * (t136 * t82 - t137 * t222) - g(3) * t83 + t12 * t215 + t14 * t97 + t2) * m(6) + t26;];
tau = t6;
