% Calculate vector of inverse dynamics joint torques for
% S5RRRPP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:08
% EndTime: 2019-12-31 20:53:18
% DurationCPUTime: 4.96s
% Computational Cost: add. (2074->381), mult. (3089->444), div. (0->0), fcn. (1360->8), ass. (0->184)
t329 = mrSges(5,2) - mrSges(4,1);
t322 = Ifges(4,5) + Ifges(6,5);
t317 = -Ifges(5,4) + t322;
t323 = Ifges(6,4) + Ifges(5,5);
t316 = -Ifges(4,6) + t323;
t328 = Ifges(6,2) + Ifges(5,3);
t164 = qJD(1) + qJD(2);
t174 = cos(qJ(3));
t257 = t164 * t174;
t107 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t257;
t242 = mrSges(5,1) * t257;
t109 = -qJD(3) * mrSges(5,3) - t242;
t326 = t109 - t107;
t171 = sin(qJ(3));
t258 = t164 * t171;
t247 = mrSges(5,1) * t258;
t325 = mrSges(4,3) * t258 + t329 * qJD(3) + t247;
t324 = -mrSges(6,1) - mrSges(5,1) - mrSges(4,3) + mrSges(3,2);
t130 = Ifges(4,4) * t257;
t277 = Ifges(6,6) * t174;
t200 = t171 * Ifges(6,3) - t277;
t321 = Ifges(4,1) * t258 + t322 * qJD(3) + t164 * t200 + t130;
t129 = Ifges(6,6) * t258;
t280 = Ifges(5,6) * t171;
t201 = -t174 * Ifges(5,3) - t280;
t320 = -Ifges(6,2) * t257 + t323 * qJD(3) + t164 * t201 + t129;
t172 = sin(qJ(2));
t276 = pkin(1) * qJD(1);
t246 = t172 * t276;
t112 = pkin(7) * t164 + t246;
t252 = qJD(3) * t171;
t175 = cos(qJ(2));
t275 = pkin(1) * qJD(2);
t240 = qJD(1) * t275;
t264 = pkin(1) * qJDD(1);
t100 = t172 * t264 + t175 * t240;
t163 = qJDD(1) + qJDD(2);
t78 = pkin(7) * t163 + t100;
t56 = t174 * t78;
t21 = -t112 * t252 + t56;
t251 = qJD(3) * t174;
t22 = -t112 * t251 - t171 * t78;
t319 = -t171 * t22 + t174 * t21;
t94 = t171 * t112;
t11 = qJD(3) * (-qJD(4) + t94) - qJDD(3) * qJ(4) - t56;
t239 = qJDD(4) - t22;
t18 = -qJDD(3) * pkin(3) + t239;
t318 = -t11 * t174 + t171 * t18;
t169 = qJ(1) + qJ(2);
t156 = sin(t169);
t157 = cos(t169);
t304 = g(1) * t157 + g(2) * t156;
t279 = Ifges(5,6) * t174;
t315 = t171 * t280 + (t277 - t279 + (-Ifges(5,2) + t328) * t171) * t174;
t245 = t175 * t276;
t113 = -pkin(2) * t164 - t245;
t209 = mrSges(4,1) * t171 + mrSges(4,2) * t174;
t267 = t174 * mrSges(5,3);
t268 = t174 * mrSges(6,2);
t158 = t171 * qJ(4);
t227 = -pkin(2) - t158;
t283 = pkin(3) + qJ(5);
t27 = -t245 + (-t174 * t283 + t227) * t164;
t161 = t174 * pkin(3);
t194 = t227 - t161;
t39 = t164 * t194 - t245;
t314 = t113 * t209 + t39 * (-t171 * mrSges(5,2) - t267) + t27 * (t171 * mrSges(6,3) - t268);
t313 = t316 * t171 + t317 * t174;
t312 = m(5) + m(6);
t86 = -t174 * t163 + t164 * t252;
t60 = mrSges(5,1) * t86 - qJDD(3) * mrSges(5,3);
t311 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t86 - t60;
t87 = t163 * t171 + t164 * t251;
t62 = t87 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t310 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t87 + t62;
t271 = t171 * mrSges(5,3);
t208 = t174 * mrSges(5,2) - t271;
t82 = t208 * t164;
t272 = t171 * mrSges(6,2);
t207 = -t174 * mrSges(6,3) - t272;
t85 = t207 * t164;
t309 = t82 + t85;
t241 = mrSges(6,1) * t257;
t110 = qJD(3) * mrSges(6,2) + t241;
t307 = t109 - t110;
t52 = -pkin(4) * t258 - t94;
t303 = -t52 + qJD(4);
t95 = t174 * t112;
t53 = pkin(4) * t257 + t95;
t302 = -t53 - qJD(5);
t99 = -t172 * t240 + t175 * t264;
t77 = -pkin(2) * t163 - t99;
t301 = m(4) * t77 + mrSges(4,1) * t86 + mrSges(4,2) * t87;
t213 = -m(6) * t283 - mrSges(6,3);
t261 = t156 * t174;
t273 = t171 * mrSges(4,2);
t300 = -t329 * t261 + t324 * t157 + (-m(5) * t194 - m(6) * t227 - t174 * t213 + mrSges(3,1) + t271 + t272 - t273) * t156;
t259 = t157 * t174;
t299 = (-mrSges(6,3) + t329) * t259 + t324 * t156 + (-mrSges(3,1) + (mrSges(4,2) - mrSges(6,2) - mrSges(5,3)) * t171) * t157;
t288 = pkin(1) * t175;
t290 = pkin(1) * t172;
t298 = m(4) * pkin(1) * (t113 * t172 + (t171 ^ 2 + t174 ^ 2) * t175 * t112) + (-mrSges(3,1) * t290 - mrSges(3,2) * t288) * t164;
t297 = (-g(1) * t259 - g(2) * t261) * qJ(4);
t69 = -qJD(3) * pkin(3) + qJD(4) + t94;
t254 = qJD(3) * qJ(4);
t79 = -t95 - t254;
t198 = t171 * t79 + t174 * t69;
t296 = m(4) * t319 + m(5) * (qJD(3) * t198 + t318) + t326 * t252 + t325 * t251;
t173 = sin(qJ(1));
t289 = pkin(1) * t173;
t287 = pkin(7) * t171;
t286 = pkin(7) * t174;
t176 = cos(qJ(1));
t162 = t176 * pkin(1);
t282 = Ifges(4,4) * t171;
t281 = Ifges(4,4) * t174;
t278 = Ifges(6,6) * t171;
t149 = pkin(7) + t290;
t263 = t149 * t171;
t262 = t149 * t174;
t256 = t157 * pkin(2) + t156 * pkin(7);
t255 = t161 + t158;
t253 = qJD(3) * t164;
t250 = qJD(4) * t171;
t244 = t175 * t275;
t154 = t172 * t275;
t243 = mrSges(6,1) * t258;
t150 = -pkin(2) - t288;
t231 = -t253 / 0.2e1;
t61 = -t86 * mrSges(6,1) + qJDD(3) * mrSges(6,2);
t145 = t157 * pkin(7);
t228 = t145 - t289;
t226 = -pkin(2) * t156 + t145;
t59 = t87 * mrSges(6,1) - qJDD(3) * mrSges(6,3);
t220 = pkin(3) * t252 - t250;
t217 = t171 * t244;
t216 = t174 * t244;
t215 = -qJD(5) - t254;
t214 = pkin(3) * t259 + t157 * t158 + t256;
t211 = qJ(5) * t174 + t255;
t115 = -mrSges(4,1) * t174 + t273;
t206 = t174 * Ifges(4,2) + t282;
t202 = -t171 * Ifges(5,2) - t279;
t93 = -pkin(2) - t211;
t193 = t156 * pkin(4) + qJ(5) * t259 + t214;
t189 = t171 * (Ifges(4,1) * t174 - t282);
t185 = (t171 * t69 - t174 * t79) * t175;
t81 = -qJ(4) * t251 + t220;
t181 = -qJ(4) * t87 + t77;
t32 = qJ(5) * t252 + t174 * t215 + t220;
t2 = t283 * t86 + (-qJD(5) * t174 - t250) * t164 + t181;
t29 = -qJD(3) * t283 + t303;
t31 = -t215 + t53;
t6 = pkin(4) * t87 - qJD(3) * qJD(5) - qJDD(3) * t283 + t239;
t7 = -pkin(4) * t86 + qJDD(5) - t11;
t70 = Ifges(4,6) * qJD(3) + t164 * t206;
t75 = Ifges(5,4) * qJD(3) + t164 * t202;
t8 = pkin(3) * t86 - t164 * t250 + t181;
t177 = (t171 * Ifges(4,1) + t200 + t281) * t87 / 0.2e1 + (-t174 * Ifges(6,2) + t201 + t278) * t86 / 0.2e1 + (t171 * t6 + t174 * t7 + t251 * t29 - t252 * t31) * mrSges(6,1) - ((-Ifges(5,6) + Ifges(6,6)) * t87 + t328 * t86 + t323 * qJDD(3)) * t174 / 0.2e1 + (t164 * (Ifges(6,3) * t174 + t278) + t320) * t252 / 0.2e1 + (t164 * (-Ifges(4,2) * t171 + t281) + t321) * t251 / 0.2e1 + ((Ifges(4,1) + Ifges(6,3)) * t87 + (-Ifges(4,4) + Ifges(6,6)) * t86 + t322 * qJDD(3)) * t171 / 0.2e1 + (t251 * t69 + t252 * t79 + t318) * mrSges(5,1) + t319 * mrSges(4,3) + t315 * t231 + (t317 * t171 - t316 * t174) * qJDD(3) / 0.2e1 + t8 * t208 - t87 * t202 / 0.2e1 - t86 * t206 / 0.2e1 + t2 * t207 + (t314 + t313 * qJD(3) / 0.2e1) * qJD(3) + t174 * (Ifges(4,4) * t87 - Ifges(4,2) * t86 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t171 * (Ifges(5,4) * qJDD(3) - Ifges(5,2) * t87 + Ifges(5,6) * t86) / 0.2e1 - t75 * t251 / 0.2e1 - t70 * t252 / 0.2e1 + t189 * t253 / 0.2e1 + Ifges(3,3) * t163 + t77 * t115 + t99 * mrSges(3,1) - t100 * mrSges(3,2);
t160 = t174 * pkin(4);
t159 = t171 * pkin(4);
t153 = pkin(4) * t251;
t146 = t157 * pkin(4);
t131 = pkin(3) * t258;
t117 = t160 + t286;
t116 = t159 + t287;
t114 = -pkin(2) - t255;
t108 = -qJD(3) * mrSges(6,3) + t243;
t105 = pkin(7) * t251 + t153;
t104 = (-pkin(4) - pkin(7)) * t252;
t98 = t160 + t262;
t97 = t159 + t263;
t96 = t150 - t255;
t84 = -qJ(4) * t257 + t131;
t83 = t115 * t164;
t76 = t93 - t288;
t54 = t154 + t81;
t49 = t131 + (-qJ(4) * t174 + qJ(5) * t171) * t164;
t48 = t149 * t251 + t153 + t217;
t47 = t216 + (-pkin(4) - t149) * t252;
t28 = t154 + t32;
t26 = -mrSges(5,2) * t86 - mrSges(5,3) * t87;
t24 = -mrSges(6,2) * t87 + mrSges(6,3) * t86;
t1 = [m(3) * (t100 * t172 + t175 * t99) * pkin(1) + t177 + m(6) * (t2 * t76 + t27 * t28 + t29 * t48 + t31 * t47 + t6 * t97 + t7 * t98) + t48 * t108 + t47 * t110 + t96 * t26 + t97 * t59 + t98 * t61 + t54 * t82 + t28 * t85 + t76 * t24 + t83 * t154 + Ifges(2,3) * qJDD(1) + m(5) * (t185 * t275 + t39 * t54 + t8 * t96) + t310 * t263 + t311 * t262 + t325 * t217 - t326 * t216 + (mrSges(3,1) * t288 - mrSges(3,2) * t290) * t163 + t301 * t150 + t298 * qJD(2) + t296 * t149 + (-m(6) * (t162 + t193) - m(3) * t162 - mrSges(2,1) * t176 + t173 * mrSges(2,2) - m(4) * (t162 + t256) - m(5) * (t162 + t214) + t299) * g(2) + (m(3) * t289 + t173 * mrSges(2,1) + mrSges(2,2) * t176 - m(4) * (t226 - t289) - m(6) * (t146 + t228) - m(5) * t228 + t300) * g(1); t104 * t110 + t105 * t108 + t114 * t26 + t116 * t59 + t117 * t61 + t93 * t24 + t32 * t85 + t81 * t82 + t177 + t310 * t287 + t311 * t286 - t301 * pkin(2) + (t104 * t31 + t105 * t29 + t116 * t6 + t117 * t7 + t2 * t93 + t27 * t32 - (t172 * t27 + (t171 * t29 + t174 * t31) * t175) * t276) * m(6) + (-(t172 * t39 + t185) * t276 + t114 * t8 + t39 * t81) * m(5) + (-t83 - t309) * t246 - t298 * qJD(1) + (-m(4) * t256 - m(5) * t214 - m(6) * t193 + t299) * g(2) + (-m(6) * (t145 + t146) - m(5) * t145 - m(4) * t226 + t300) * g(1) + t296 * pkin(7) + ((-t108 - t325) * t171 + (-t107 + t307) * t174) * t245; (qJ(4) * t7 - t27 * t49 - t283 * t6 + t29 * t302 + t303 * t31 + t297) * m(6) - t283 * t59 - (Ifges(6,3) * t257 + t129 + t320) * t258 / 0.2e1 - (-Ifges(4,2) * t258 + t130 + t321) * t257 / 0.2e1 - t325 * t95 - t326 * t94 + (-t267 - t268 + t209 + (m(5) * pkin(3) - mrSges(5,2) - t213) * t171) * t304 + t313 * t231 + t316 * t86 + t317 * t87 + (-t314 + (t315 / 0.2e1 - t189 / 0.2e1) * t164) * t164 + t75 * t257 / 0.2e1 + t70 * t258 / 0.2e1 - t29 * t241 - t69 * t242 - t307 * qJD(4) + t302 * t108 + (Ifges(5,1) + Ifges(4,3) + Ifges(6,1)) * qJDD(3) + (-t60 + t61) * qJ(4) + (-m(5) * t255 - m(6) * t211 + t115 + t207 + t208) * g(3) - t79 * t247 - t52 * t110 - t84 * t82 - t49 * t85 - pkin(3) * t62 + t18 * mrSges(5,2) - t21 * mrSges(4,2) + t22 * mrSges(4,1) - t6 * mrSges(6,3) + t7 * mrSges(6,2) - t11 * mrSges(5,3) + t31 * t243 + (-pkin(3) * t18 - qJ(4) * t11 - qJD(4) * t79 - t112 * t198 - t39 * t84 + t297) * m(5); t307 * qJD(3) + t312 * t174 * g(3) + (t309 * t164 - t304 * t312) * t171 + t59 + t62 + (-qJD(3) * t31 + t258 * t27 + t6) * m(6) + (t79 * qJD(3) + t258 * t39 + t18) * m(5); t85 * t257 + qJD(3) * t108 + (-g(3) * t171 + t29 * qJD(3) - t174 * t304 + t257 * t27 + t7) * m(6) + t61;];
tau = t1;
