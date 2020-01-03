% Calculate vector of inverse dynamics joint torques for
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:32
% EndTime: 2020-01-03 12:03:40
% DurationCPUTime: 4.20s
% Computational Cost: add. (5951->392), mult. (9273->515), div. (0->0), fcn. (6441->16), ass. (0->190)
t227 = cos(qJ(2));
t275 = qJD(1) * pkin(1);
t257 = t227 * t275;
t241 = qJD(3) - t257;
t218 = sin(pkin(9));
t222 = sin(qJ(4));
t219 = cos(pkin(9));
t226 = cos(qJ(4));
t265 = t219 * t226;
t157 = -t218 * t222 + t265;
t220 = -pkin(7) - qJ(3);
t171 = t220 * t218;
t205 = t219 * pkin(7);
t172 = qJ(3) * t219 + t205;
t259 = qJD(4) * t226;
t298 = t171 * t259 + qJD(3) * t265 + (-qJD(3) * t218 - qJD(4) * t172) * t222 - t157 * t257;
t113 = t222 * t171 + t226 * t172;
t158 = t218 * t226 + t219 * t222;
t297 = -t113 * qJD(4) - t241 * t158;
t144 = t157 * qJD(4);
t287 = pkin(8) * t144;
t312 = -t287 + t297;
t145 = t158 * qJD(4);
t142 = t145 * pkin(8);
t311 = t142 - t298;
t214 = pkin(9) + qJ(4);
t200 = sin(t214);
t201 = cos(t214);
t202 = qJ(5) + t214;
t190 = sin(t202);
t279 = mrSges(6,2) * t190;
t310 = -mrSges(5,1) * t201 + mrSges(5,2) * t200 + t279;
t240 = -mrSges(4,1) * t219 + mrSges(4,2) * t218;
t209 = qJDD(4) + qJDD(5);
t215 = qJD(4) + qJD(5);
t225 = cos(qJ(5));
t221 = sin(qJ(5));
t216 = qJD(1) + qJD(2);
t131 = t157 * t216;
t223 = sin(qJ(2));
t255 = t223 * t275;
t165 = qJ(3) * t216 + t255;
t252 = pkin(7) * t216 + t165;
t119 = t252 * t218;
t120 = t252 * t219;
t71 = -t119 * t222 + t120 * t226;
t53 = pkin(8) * t131 + t71;
t274 = t221 * t53;
t132 = t158 * t216;
t270 = t120 * t222;
t70 = -t226 * t119 - t270;
t52 = -pkin(8) * t132 + t70;
t50 = qJD(4) * pkin(4) + t52;
t22 = t225 * t50 - t274;
t273 = t225 * t53;
t23 = t221 * t50 + t273;
t246 = t225 * t131 - t132 * t221;
t86 = t131 * t221 + t132 * t225;
t290 = Ifges(6,4) * t86;
t254 = qJD(2) * t275;
t271 = qJDD(1) * pkin(1);
t153 = t223 * t271 + t227 * t254;
t210 = qJDD(1) + qJDD(2);
t118 = qJ(3) * t210 + qJD(3) * t216 + t153;
t253 = pkin(7) * t210 + t118;
t100 = t253 * t218;
t101 = t253 * t219;
t37 = -t71 * qJD(4) - t226 * t100 - t101 * t222;
t91 = t216 * t144 + t158 * t210;
t18 = qJDD(4) * pkin(4) - pkin(8) * t91 + t37;
t36 = -qJD(4) * t270 - t222 * t100 + t226 * t101 - t119 * t259;
t92 = -t216 * t145 + t157 * t210;
t19 = pkin(8) * t92 + t36;
t3 = t22 * qJD(5) + t18 * t221 + t19 * t225;
t34 = t246 * qJD(5) + t221 * t92 + t225 * t91;
t35 = -t86 * qJD(5) - t221 * t91 + t225 * t92;
t4 = -t23 * qJD(5) + t18 * t225 - t19 * t221;
t78 = Ifges(6,4) * t246;
t43 = Ifges(6,1) * t86 + Ifges(6,5) * t215 + t78;
t193 = pkin(3) * t219 + pkin(2);
t129 = -t193 * t216 + t241;
t93 = -pkin(4) * t131 + t129;
t309 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t34 + Ifges(6,6) * t35 + Ifges(6,3) * t209 - (Ifges(6,5) * t246 - Ifges(6,6) * t86) * t215 / 0.2e1 + (t22 * t246 + t23 * t86) * mrSges(6,3) - (-Ifges(6,2) * t86 + t43 + t78) * t246 / 0.2e1 - t93 * (mrSges(6,1) * t86 + mrSges(6,2) * t246) - (Ifges(6,1) * t246 - t290) * t86 / 0.2e1;
t308 = -mrSges(5,3) - mrSges(6,3);
t260 = t218 ^ 2 + t219 ^ 2;
t250 = t260 * t118;
t307 = -mrSges(3,1) + t240 + t310;
t42 = Ifges(6,2) * t246 + Ifges(6,6) * t215 + t290;
t305 = t42 / 0.2e1;
t112 = t226 * t171 - t172 * t222;
t286 = pkin(8) * t158;
t94 = t112 - t286;
t150 = t157 * pkin(8);
t95 = t150 + t113;
t45 = -t221 * t95 + t225 * t94;
t304 = t45 * qJD(5) + t312 * t221 - t311 * t225;
t46 = t221 * t94 + t225 * t95;
t303 = -t46 * qJD(5) + t311 * t221 + t312 * t225;
t192 = pkin(1) * t223 + qJ(3);
t147 = (-pkin(7) - t192) * t218;
t148 = t192 * t219 + t205;
t103 = t222 * t147 + t226 * t148;
t294 = t86 / 0.2e1;
t292 = t132 / 0.2e1;
t289 = pkin(1) * t227;
t288 = pkin(4) * t145;
t217 = qJ(1) + qJ(2);
t203 = sin(t217);
t285 = g(2) * t203;
t204 = cos(t217);
t284 = g(3) * t204;
t281 = mrSges(6,1) * t190;
t277 = Ifges(5,4) * t132;
t276 = pkin(1) * qJD(2);
t191 = cos(t202);
t178 = t191 * mrSges(6,1);
t272 = mrSges(5,1) * t131 - mrSges(5,2) * t132 - t240 * t216;
t268 = t191 * t204;
t267 = t210 * t218;
t266 = t210 * t219;
t160 = pkin(4) * t201 + t193;
t211 = -pkin(8) + t220;
t264 = t203 * t160 + t204 * t211;
t263 = mrSges(6,2) * t268 + t204 * t281;
t262 = t203 * t193 + t204 * t220;
t261 = t204 * pkin(2) + t203 * qJ(3);
t256 = t223 * t276;
t47 = -t92 * mrSges(5,1) + t91 * mrSges(5,2);
t9 = -t35 * mrSges(6,1) + t34 * mrSges(6,2);
t251 = t260 * mrSges(4,3);
t249 = t260 * t165;
t248 = t260 * t210;
t247 = t260 * t216;
t102 = t226 * t147 - t148 * t222;
t245 = t204 * t160 - t203 * t211;
t244 = t204 * t193 - t203 * t220;
t243 = qJD(3) * t260;
t242 = m(4) * qJ(3) - t308;
t139 = -mrSges(4,1) * t266 + mrSges(4,2) * t267;
t152 = -t223 * t254 + t227 * t271;
t239 = -mrSges(5,1) * t200 - mrSges(5,2) * t201;
t238 = -mrSges(6,2) * t191 - t281;
t76 = t102 - t286;
t77 = t150 + t103;
t38 = -t221 * t77 + t225 * t76;
t39 = t221 * t76 + t225 * t77;
t105 = t157 * t225 - t158 * t221;
t106 = t157 * t221 + t158 * t225;
t128 = -pkin(4) * t157 - t193;
t235 = qJDD(3) - t152;
t179 = t227 * t276 + qJD(3);
t60 = t147 * t259 + t179 * t265 + (-qJD(4) * t148 - t179 * t218) * t222;
t110 = -t193 * t210 + t235;
t231 = -t204 * mrSges(3,2) + (-t178 + t307) * t203;
t61 = -t103 * qJD(4) - t158 * t179;
t230 = -mrSges(6,1) * t268 + t307 * t204 + (mrSges(3,2) - mrSges(4,3) + t308) * t203;
t130 = -pkin(2) * t210 + t235;
t56 = t105 * qJD(5) + t144 * t225 - t145 * t221;
t57 = -t106 * qJD(5) - t144 * t221 - t145 * t225;
t58 = -pkin(4) * t92 + t110;
t81 = t131 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t277;
t127 = Ifges(5,4) * t131;
t82 = Ifges(5,1) * t132 + Ifges(5,5) * qJD(4) + t127;
t229 = mrSges(4,3) * t250 + (t110 * mrSges(5,2) - t37 * mrSges(5,3) + Ifges(5,1) * t91 + Ifges(5,4) * t92 + Ifges(5,5) * qJDD(4)) * t158 + (-t110 * mrSges(5,1) + t36 * mrSges(5,3) + Ifges(5,4) * t91 + Ifges(5,2) * t92 + Ifges(5,6) * qJDD(4)) * t157 + (t58 * mrSges(6,2) - t4 * mrSges(6,3) + Ifges(6,1) * t34 + Ifges(6,4) * t35 + Ifges(6,5) * t209) * t106 + (-t58 * mrSges(6,1) + t3 * mrSges(6,3) + Ifges(6,4) * t34 + Ifges(6,2) * t35 + Ifges(6,6) * t209) * t105 + (Ifges(4,4) * t218 + Ifges(4,2) * t219) * t266 + (Ifges(4,1) * t218 + Ifges(4,4) * t219) * t267 + t130 * t240 + (-t144 * t70 - t71 * t145) * mrSges(5,3) + (Ifges(5,1) * t144 - Ifges(5,4) * t145) * t292 + qJD(4) * (Ifges(5,5) * t144 - Ifges(5,6) * t145) / 0.2e1 + t129 * (mrSges(5,1) * t145 + mrSges(5,2) * t144) + t131 * (Ifges(5,4) * t144 - Ifges(5,2) * t145) / 0.2e1 + t246 * (Ifges(6,4) * t56 + Ifges(6,2) * t57) / 0.2e1 + (-t22 * t56 + t23 * t57) * mrSges(6,3) + (Ifges(6,1) * t56 + Ifges(6,4) * t57) * t294 + t57 * t305 + t56 * t43 / 0.2e1 + t93 * (-mrSges(6,1) * t57 + mrSges(6,2) * t56) + t144 * t82 / 0.2e1 - t145 * t81 / 0.2e1 + t152 * mrSges(3,1) - t153 * mrSges(3,2) + Ifges(3,3) * t210 + t215 * (Ifges(6,5) * t56 + Ifges(6,6) * t57) / 0.2e1;
t228 = cos(qJ(1));
t224 = sin(qJ(1));
t207 = t228 * pkin(1);
t206 = t224 * pkin(1);
t196 = -pkin(2) - t289;
t194 = t203 * pkin(2);
t170 = -t193 - t289;
t159 = -pkin(2) * t216 + t241;
t121 = t256 + t288;
t117 = t128 - t289;
t116 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t132;
t115 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t131;
t80 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t92;
t79 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t91;
t73 = mrSges(6,1) * t215 - mrSges(6,3) * t86;
t72 = -mrSges(6,2) * t215 + mrSges(6,3) * t246;
t49 = t61 - t287;
t48 = -t142 + t60;
t44 = -mrSges(6,1) * t246 + mrSges(6,2) * t86;
t29 = -mrSges(6,2) * t209 + mrSges(6,3) * t35;
t28 = mrSges(6,1) * t209 - mrSges(6,3) * t34;
t25 = t225 * t52 - t274;
t24 = -t221 * t52 - t273;
t8 = -t39 * qJD(5) - t221 * t48 + t225 * t49;
t7 = t38 * qJD(5) + t221 * t49 + t225 * t48;
t1 = [t229 + m(6) * (t117 * t58 + t121 * t93 + t22 * t8 + t23 * t7 + t3 * t39 + t38 * t4) + m(5) * (t102 * t37 + t103 * t36 + t110 * t170 + t129 * t256 + t60 * t71 + t61 * t70) + t38 * t28 + t39 * t29 + m(4) * (t130 * t196 + t159 * t256 + t179 * t249 + t192 * t250) + t7 * t72 + t8 * t73 + t102 * t79 + t103 * t80 + t60 * t115 + t61 * t116 + t117 * t9 + t121 * t44 + t170 * t47 + t196 * t139 + Ifges(2,3) * qJDD(1) + (t179 * t247 + t192 * t248) * mrSges(4,3) + (-m(6) * (t206 + t264) - m(5) * (t206 + t262) - m(4) * (t194 + t206) - t224 * mrSges(2,1) - mrSges(2,2) * t228 + (mrSges(4,3) + t242) * t204 + t231) * g(3) + (-m(6) * (t207 + t245) - m(4) * (t207 + t261) - m(5) * (t207 + t244) - mrSges(2,1) * t228 + t224 * mrSges(2,2) + t230) * g(2) + ((mrSges(3,1) * t227 - mrSges(3,2) * t223) * t210 + (-mrSges(3,2) * t216 * t227 + (-mrSges(3,1) * t216 - t272) * t223) * qJD(2) + (-g(2) * t228 - g(3) * t224 + t152 * t227 + t153 * t223) * m(3)) * pkin(1); t229 + t297 * t116 + t298 * t115 + ((-t44 + t272) * t223 + (mrSges(3,1) * t223 + (mrSges(3,2) - t251) * t227) * t216) * t275 + t303 * t73 + t304 * t72 + t45 * t28 + t46 * t29 + t44 * t288 + (t242 * t204 + t231) * g(3) + (qJ(3) * t248 + t216 * t243 + t284) * mrSges(4,3) + t230 * g(2) + t112 * t79 + t113 * t80 + t128 * t9 - pkin(2) * t139 - t193 * t47 + (-t245 * g(2) - t264 * g(3) + t128 * t58 + t3 * t46 + t4 * t45 + (-t255 + t288) * t93 + t304 * t23 + t303 * t22) * m(6) + (-t244 * g(2) - t262 * g(3) - t110 * t193 + t112 * t37 + t113 * t36 - t129 * t255 + t297 * t70 + t298 * t71) * m(5) + (-pkin(2) * t130 + qJ(3) * t250 + t165 * t243 - (t159 * t223 + t227 * t249) * t275 - t194 * g(3) - t261 * g(2)) * m(4); -t216 ^ 2 * t251 - t131 * t115 + t132 * t116 - t246 * t72 + t86 * t73 + t139 + t47 + t9 + (g(2) * t204 + g(3) * t203) * (m(4) + m(5) + m(6)) + (t22 * t86 - t23 * t246 + t58) * m(6) + (-t131 * t71 + t132 * t70 + t110) * m(5) + (-t165 * t247 + t130) * m(4); (-t178 + t310) * g(1) - t132 * (Ifges(5,1) * t131 - t277) / 0.2e1 - (-Ifges(5,2) * t132 + t127 + t82) * t131 / 0.2e1 + t86 * t305 + (t239 * t204 - t263) * g(3) + (-t238 - t239) * t285 + t81 * t292 + (-t132 * t44 + t221 * t29 + t225 * t28 + (-t132 * t93 + t221 * t3 + t225 * t4 - g(1) * t201 + (-t284 + t285) * t200) * m(6) + (-t221 * t73 + t225 * t72 + (-t22 * t221 + t225 * t23) * m(6)) * qJD(5)) * pkin(4) + t37 * mrSges(5,1) - t36 * mrSges(5,2) - m(6) * (t22 * t24 + t23 * t25) + t309 - t25 * t72 - t24 * t73 + Ifges(5,5) * t91 + Ifges(5,6) * t92 - t70 * t115 + t71 * t116 - qJD(4) * (Ifges(5,5) * t131 - Ifges(5,6) * t132) / 0.2e1 - t129 * (mrSges(5,1) * t132 + mrSges(5,2) * t131) + Ifges(5,3) * qJDD(4) + (t131 * t70 + t132 * t71) * mrSges(5,3); t42 * t294 - t22 * t72 + t23 * t73 - g(1) * (t178 - t279) - g(3) * t263 - t238 * t285 + t309;];
tau = t1;
