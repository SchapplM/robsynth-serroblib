% Calculate vector of inverse dynamics joint torques for
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:03
% EndTime: 2020-01-03 12:13:10
% DurationCPUTime: 2.96s
% Computational Cost: add. (5576->402), mult. (7695->542), div. (0->0), fcn. (4376->16), ass. (0->210)
t200 = sin(qJ(4));
t289 = t200 / 0.2e1;
t205 = cos(qJ(4));
t280 = mrSges(5,1) * t205;
t299 = -mrSges(4,1) - t280;
t202 = sin(qJ(2));
t272 = qJD(1) * pkin(1);
t247 = t202 * t272;
t207 = cos(qJ(2));
t287 = pkin(1) * t207;
t126 = -qJD(2) * t247 + qJDD(1) * t287;
t192 = qJDD(1) + qJDD(2);
t298 = pkin(2) * t192 - qJD(3) * t247 + t126;
t297 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t230 = mrSges(5,1) * t200 + mrSges(5,2) * t205;
t194 = qJD(1) + qJD(2);
t182 = qJD(3) + t194;
t246 = t207 * t272;
t142 = pkin(2) * t194 + t246;
t201 = sin(qJ(3));
t153 = t201 * t247;
t206 = cos(qJ(3));
t98 = t142 * t206 - t153;
t83 = -pkin(3) * t182 - t98;
t296 = t83 * t230 + qJD(4) * (Ifges(5,5) * t205 - Ifges(5,6) * t200) / 0.2e1;
t199 = sin(qJ(5));
t204 = cos(qJ(5));
t131 = -t199 * t200 + t204 * t205;
t209 = -pkin(9) - pkin(8);
t242 = qJD(4) * t209;
t139 = t200 * t242;
t140 = t205 * t242;
t149 = t209 * t200;
t189 = t205 * pkin(9);
t150 = pkin(8) * t205 + t189;
t91 = t149 * t204 - t150 * t199;
t295 = qJD(5) * t91 - t131 * t98 + t139 * t204 + t140 * t199;
t132 = t199 * t205 + t200 * t204;
t92 = t149 * t199 + t150 * t204;
t294 = -qJD(5) * t92 + t132 * t98 - t139 * t199 + t140 * t204;
t115 = t206 * t246 - t153;
t173 = pkin(2) * t201 + pkin(8);
t281 = -pkin(9) - t173;
t124 = t281 * t200;
t262 = t173 * t205;
t125 = t189 + t262;
t77 = t124 * t199 + t125 * t204;
t238 = qJD(4) * t281;
t250 = qJD(3) * t206;
t245 = pkin(2) * t250;
t93 = t200 * t238 + t205 * t245;
t94 = -t200 * t245 + t205 * t238;
t293 = -qJD(5) * t77 + t132 * t115 - t199 * t93 + t204 * t94;
t76 = t124 * t204 - t125 * t199;
t292 = qJD(5) * t76 - t131 * t115 + t199 * t94 + t204 * t93;
t99 = t142 * t201 + t206 * t247;
t84 = pkin(8) * t182 + t99;
t237 = (t200 ^ 2 + t205 ^ 2) * t84;
t103 = t132 * t182;
t290 = t103 / 0.2e1;
t286 = pkin(2) * t206;
t285 = pkin(4) * t205;
t198 = qJ(1) + qJ(2);
t187 = qJ(3) + t198;
t171 = sin(t187);
t284 = g(2) * t171;
t248 = qJD(4) * t205;
t181 = qJDD(3) + t192;
t252 = qJD(2) * t207;
t127 = (qJD(1) * t252 + qJDD(1) * t202) * pkin(1);
t50 = t206 * t127 + t142 * t250 + t201 * t298;
t40 = pkin(8) * t181 + t50;
t17 = -t200 * t40 - t248 * t84;
t283 = t17 * mrSges(5,3);
t176 = pkin(2) + t287;
t258 = t202 * t206;
t119 = pkin(1) * t258 + t201 * t176;
t113 = pkin(8) + t119;
t282 = -pkin(9) - t113;
t279 = mrSges(5,2) * t200;
t197 = qJ(4) + qJ(5);
t183 = sin(t197);
t278 = mrSges(6,2) * t183;
t102 = t131 * t182;
t277 = mrSges(6,3) * t102;
t276 = Ifges(5,4) * t200;
t275 = Ifges(5,4) * t205;
t274 = Ifges(6,4) * t103;
t273 = Ifges(5,2) * t205;
t271 = qJD(4) * pkin(4);
t270 = t103 * mrSges(6,3);
t249 = qJD(4) * t200;
t16 = t205 * t40 - t249 * t84;
t269 = t16 * t205;
t268 = t17 * t200;
t185 = cos(t197);
t166 = t185 * mrSges(6,1);
t241 = pkin(9) * t182 + t84;
t69 = t241 * t205;
t267 = t199 * t69;
t266 = t204 * t69;
t172 = cos(t187);
t264 = t172 * t183;
t263 = t172 * t185;
t261 = t182 * t200;
t260 = t182 * t205;
t259 = t201 * t202;
t129 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t260;
t257 = t205 * t129;
t174 = pkin(3) + t285;
t256 = t171 * t174 + t172 * t209;
t255 = mrSges(6,1) * t264 + mrSges(6,2) * t263;
t254 = t172 * pkin(3) + t171 * pkin(8);
t184 = sin(t198);
t169 = pkin(2) * t184;
t203 = sin(qJ(1));
t188 = t203 * pkin(1);
t253 = t169 + t188;
t251 = qJD(3) * t201;
t178 = pkin(4) * t249;
t244 = t169 + t256;
t186 = cos(t198);
t170 = pkin(2) * t186;
t243 = t170 + t254;
t239 = qJD(4) * t282;
t147 = t279 - t280;
t108 = t147 * t182;
t236 = -t182 * mrSges(4,1) + t108;
t235 = t166 - t278;
t234 = -t171 * t209 + t172 * t174;
t118 = -pkin(1) * t259 + t176 * t206;
t64 = -mrSges(6,1) * t102 + mrSges(6,2) * t103;
t232 = -t236 - t64;
t112 = -pkin(3) - t118;
t231 = t170 + t234;
t68 = t241 * t200;
t229 = -mrSges(6,1) * t183 - mrSges(6,2) * t185;
t228 = t278 + t279;
t227 = t273 + t276;
t225 = -t268 + t269;
t67 = -t68 + t271;
t21 = t204 * t67 - t267;
t22 = t199 * t67 + t266;
t89 = t282 * t200;
t90 = t113 * t205 + t189;
t56 = -t199 * t90 + t204 * t89;
t57 = t199 * t89 + t204 * t90;
t128 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t261;
t224 = t128 * t205 + t129 * t200;
t223 = t201 * t207 + t258;
t51 = -t201 * t127 - t142 * t251 + t206 * t298;
t221 = t200 * (Ifges(5,1) * t205 - t276);
t220 = t132 * qJD(5);
t219 = t131 * qJD(5);
t218 = -mrSges(6,1) * t263 + t171 * t297 + t172 * t299;
t217 = t182 * mrSges(4,2) + t200 * t128 - t257;
t41 = -pkin(3) * t181 - t51;
t216 = (m(5) * pkin(8) - t297) * t172 + (-t166 + t299) * t171;
t109 = t181 * t205 - t182 * t249;
t95 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t109;
t110 = t181 * t200 + t182 * t248;
t96 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t110;
t215 = -qJD(4) * t224 - t200 * t96 + t205 * t95;
t191 = qJDD(4) + qJDD(5);
t193 = qJD(4) + qJD(5);
t10 = qJDD(4) * pkin(4) - pkin(9) * t110 + t17;
t11 = pkin(9) * t109 + t16;
t3 = qJD(5) * t21 + t10 * t199 + t11 * t204;
t4 = -qJD(5) * t22 + t10 * t204 - t11 * t199;
t44 = t109 * t199 + t110 * t204 + t182 * t219;
t45 = t109 * t204 - t110 * t199 - t182 * t220;
t59 = Ifges(6,2) * t102 + Ifges(6,6) * t193 + t274;
t97 = Ifges(6,4) * t102;
t60 = Ifges(6,1) * t103 + Ifges(6,5) * t193 + t97;
t73 = -t174 * t182 - t98;
t214 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t21 * t277 + t59 * t290 - t73 * (mrSges(6,1) * t103 + mrSges(6,2) * t102) + Ifges(6,3) * t191 - t103 * (Ifges(6,1) * t102 - t274) / 0.2e1 + Ifges(6,6) * t45 + Ifges(6,5) * t44 - t193 * (Ifges(6,5) * t102 - Ifges(6,6) * t103) / 0.2e1 - (-Ifges(6,2) * t103 + t60 + t97) * t102 / 0.2e1;
t75 = t176 * t251 + (qJD(2) * t223 + t202 * t250) * pkin(1);
t213 = -t186 * mrSges(3,1) + mrSges(3,2) * t184 + t172 * t228 + t218;
t212 = -t184 * mrSges(3,1) - t186 * mrSges(3,2) + t171 * t228 + t216;
t100 = Ifges(5,6) * qJD(4) + t182 * t227;
t148 = Ifges(5,4) * t260;
t101 = Ifges(5,1) * t261 + Ifges(5,5) * qJD(4) + t148;
t27 = -pkin(4) * t109 + t41;
t81 = qJD(4) * t131 + t219;
t82 = -qJD(4) * t132 - t220;
t211 = t205 * (Ifges(5,4) * t110 + Ifges(5,2) * t109) / 0.2e1 + t193 * (Ifges(6,5) * t81 + Ifges(6,6) * t82) / 0.2e1 + Ifges(4,3) * t181 + t41 * t147 + t102 * (Ifges(6,4) * t81 + Ifges(6,2) * t82) / 0.2e1 + t109 * t227 / 0.2e1 + t110 * (t200 * Ifges(5,1) + t275) / 0.2e1 + t81 * t60 / 0.2e1 + t82 * t59 / 0.2e1 + t73 * (-mrSges(6,1) * t82 + mrSges(6,2) * t81) + t51 * mrSges(4,1) - t100 * t249 / 0.2e1 + (Ifges(5,1) * t110 + Ifges(5,4) * t109) * t289 + (Ifges(6,1) * t81 + Ifges(6,4) * t82) * t290 + mrSges(5,3) * t269 + (t101 + t182 * (-Ifges(5,2) * t200 + t275)) * t248 / 0.2e1 + (-t21 * t81 + t22 * t82) * mrSges(6,3) + (0.2e1 * Ifges(5,5) * t289 + Ifges(5,6) * t205) * qJDD(4) + (t221 * t182 / 0.2e1 + t296) * qJD(4) + (mrSges(6,2) * t27 - mrSges(6,3) * t4 + Ifges(6,1) * t44 + Ifges(6,4) * t45 + Ifges(6,5) * t191) * t132 + (-mrSges(6,1) * t27 + mrSges(6,3) * t3 + Ifges(6,4) * t44 + Ifges(6,2) * t45 + Ifges(6,6) * t191) * t131;
t210 = t126 * mrSges(3,1) + Ifges(3,3) * t192 + t211;
t208 = cos(qJ(1));
t190 = t208 * pkin(1);
t175 = -pkin(3) - t286;
t160 = t171 * pkin(3);
t145 = -t174 - t286;
t141 = pkin(2) * t251 + t178;
t114 = t223 * t272;
t107 = t112 - t285;
t86 = mrSges(6,1) * t193 - t270;
t85 = -mrSges(6,2) * t193 + t277;
t74 = t176 * t250 + (-t202 * t251 + (t206 * t207 - t259) * qJD(2)) * pkin(1);
t72 = t178 + t75;
t71 = -mrSges(5,1) * t109 + mrSges(5,2) * t110;
t48 = -t200 * t74 + t205 * t239;
t47 = t200 * t239 + t205 * t74;
t36 = -mrSges(6,2) * t191 + mrSges(6,3) * t45;
t35 = mrSges(6,1) * t191 - mrSges(6,3) * t44;
t26 = -t204 * t68 - t267;
t25 = t199 * t68 - t266;
t9 = -mrSges(6,1) * t45 + mrSges(6,2) * t44;
t6 = -qJD(5) * t57 - t199 * t47 + t204 * t48;
t5 = qJD(5) * t56 + t199 * t48 + t204 * t47;
t1 = [t74 * t257 - t127 * mrSges(3,2) + t75 * t108 + t112 * t71 + t107 * t9 + (-mrSges(2,1) * t208 + mrSges(2,2) * t203 - m(6) * (t190 + t231) - m(5) * (t190 + t243) - m(4) * (t170 + t190) + t213) * g(2) + (-mrSges(2,1) * t203 - mrSges(2,2) * t208 - m(5) * (t160 + t253) - m(6) * (t188 + t244) - m(4) * t253 + t212) * g(3) + (t118 * t181 - t182 * t75) * mrSges(4,1) + t215 * t113 + t5 * t85 + t6 * t86 + t72 * t64 + t56 * t35 + t57 * t36 + m(4) * (t118 * t51 + t119 * t50 + t74 * t99 - t75 * t98) + m(6) * (t107 * t27 + t21 * t6 + t22 * t5 + t3 * t57 + t4 * t56 + t72 * t73) + (-t74 * t128 - t283) * t200 + t210 + (-t119 * t181 - t182 * t74 - t50) * mrSges(4,2) + Ifges(2,3) * qJDD(1) + m(5) * (t112 * t41 + t113 * t225 + t237 * t74 + t75 * t83) + ((-t192 * t202 - t194 * t252) * mrSges(3,2) + (-qJD(2) * t194 * t202 + t192 * t207) * mrSges(3,1) + (-g(2) * t208 - g(3) * t203 + t126 * t207 + t127 * t202) * m(3)) * pkin(1); -mrSges(5,3) * t268 + t175 * t71 + t141 * t64 + t145 * t9 - m(5) * (t114 * t83 + t115 * t237) + (-m(5) * t243 + t213) * g(2) + (t194 * t246 - t127) * mrSges(3,2) + t76 * t35 + t77 * t36 - t50 * mrSges(4,2) + t215 * t173 + (-m(5) * (t160 + t169) + t212) * g(3) + t217 * t115 + t232 * t114 - m(4) * (-t114 * t98 + t115 * t99) + t210 + t293 * t86 + m(5) * (t16 * t262 - t173 * t268 + t175 * t41) + t292 * t85 + t194 * mrSges(3,1) * t247 + ((mrSges(4,1) * t206 - mrSges(4,2) * t201) * t181 + (-g(2) * t186 - g(3) * t184 + t201 * t50 + t206 * t51) * m(4) + (t236 * t201 - t217 * t206 + m(5) * (t201 * t83 + t206 * t237) + m(4) * (-t201 * t98 + t206 * t99)) * qJD(3)) * pkin(2) + (-t231 * g(2) - t244 * g(3) + t145 * t27 + t3 * t77 + t4 * t76 + (t141 - t114) * t73 + t292 * t22 + t293 * t21) * m(6); -t174 * t9 + t91 * t35 + t92 * t36 + (t64 * t271 - t283 + t98 * t128 + (-qJD(4) * t129 - t96) * pkin(8) + (g(2) * t172 + g(3) * t171) * mrSges(5,2)) * t200 + (t171 * t278 + t216) * g(3) + (-t98 * t129 + (-qJD(4) * t128 + t95) * pkin(8)) * t205 - pkin(3) * t71 + t294 * t86 + t295 * t85 + (t182 * t98 - t50) * mrSges(4,2) + t211 + (mrSges(6,2) * t264 + t218) * g(2) + t232 * t99 + (-t234 * g(2) - t256 * g(3) - t174 * t27 + t3 * t92 + t4 * t91 + (t178 - t99) * t73 + t295 * t22 + t294 * t21) * m(6) + (-pkin(3) * t41 + pkin(8) * t225 - t254 * g(2) - t160 * g(3) - t237 * t98 - t83 * t99) * m(5); (t100 * t289 + (t273 * t289 - t221 / 0.2e1) * t182 - (t148 + t101) * t205 / 0.2e1 - t296) * t182 + Ifges(5,6) * t109 + Ifges(5,5) * t110 - t26 * t85 - t25 * t86 - t16 * mrSges(5,2) + t17 * mrSges(5,1) + (t147 - t235) * g(1) + t22 * t270 + t214 + (-t229 + t230) * t284 + (-t172 * t230 - t255) * g(3) + Ifges(5,3) * qJDD(4) + t224 * t84 + (-t64 * t261 + t199 * t36 + t204 * t35 + (t199 * t3 + t204 * t4 - g(1) * t205 + (-g(3) * t172 - t182 * t73 + t284) * t200) * m(6) + (-t199 * t86 + t204 * t85 + (-t199 * t21 + t204 * t22) * m(6)) * qJD(5)) * pkin(4) - m(6) * (t21 * t25 + t22 * t26); (t86 + t270) * t22 - t229 * t284 - g(1) * t235 - g(3) * t255 - t21 * t85 + t214;];
tau = t1;
