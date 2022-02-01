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
% m [6x1]
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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 12:01:41
% EndTime: 2022-01-20 12:01:50
% DurationCPUTime: 3.02s
% Computational Cost: add. (5576->390), mult. (7695->525), div. (0->0), fcn. (4376->16), ass. (0->205)
t192 = sin(qJ(4));
t278 = t192 / 0.2e1;
t194 = sin(qJ(2));
t258 = pkin(1) * qJD(1);
t238 = t194 * t258;
t199 = cos(qJ(2));
t276 = pkin(1) * t199;
t126 = -qJD(2) * t238 + qJDD(1) * t276;
t184 = qJDD(1) + qJDD(2);
t292 = pkin(2) * t184 - qJD(3) * t238 + t126;
t189 = qJ(4) + qJ(5);
t176 = sin(t189);
t265 = mrSges(6,2) * t176;
t266 = mrSges(5,2) * t192;
t291 = t265 + t266;
t190 = qJ(1) + qJ(2);
t180 = qJ(3) + t190;
t164 = sin(t180);
t165 = cos(t180);
t283 = g(1) * t165 + g(2) * t164;
t290 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t178 = cos(t189);
t161 = t178 * mrSges(6,1);
t197 = cos(qJ(4));
t267 = mrSges(5,1) * t197;
t289 = -mrSges(4,1) - t161 - t267;
t222 = mrSges(5,1) * t192 + mrSges(5,2) * t197;
t186 = qJD(1) + qJD(2);
t175 = qJD(3) + t186;
t237 = t199 * t258;
t139 = pkin(2) * t186 + t237;
t193 = sin(qJ(3));
t149 = t193 * t238;
t198 = cos(qJ(3));
t98 = t139 * t198 - t149;
t83 = -pkin(3) * t175 - t98;
t288 = t83 * t222 + qJD(4) * (Ifges(5,5) * t197 - Ifges(5,6) * t192) / 0.2e1;
t191 = sin(qJ(5));
t196 = cos(qJ(5));
t132 = t191 * t197 + t192 * t196;
t201 = -pkin(9) - pkin(8);
t234 = qJD(4) * t201;
t136 = t192 * t234;
t137 = t197 * t234;
t146 = t201 * t192;
t181 = t197 * pkin(9);
t147 = pkin(8) * t197 + t181;
t92 = t146 * t191 + t147 * t196;
t287 = -qJD(5) * t92 + t132 * t98 - t136 * t191 + t137 * t196;
t131 = -t191 * t192 + t196 * t197;
t91 = t146 * t196 - t147 * t191;
t286 = qJD(5) * t91 - t131 * t98 + t136 * t196 + t137 * t191;
t115 = t198 * t237 - t149;
t166 = pkin(2) * t193 + pkin(8);
t268 = -pkin(9) - t166;
t124 = t268 * t192;
t251 = t166 * t197;
t125 = t181 + t251;
t77 = t124 * t191 + t125 * t196;
t230 = qJD(4) * t268;
t242 = qJD(3) * t198;
t236 = pkin(2) * t242;
t93 = t192 * t230 + t197 * t236;
t94 = -t192 * t236 + t197 * t230;
t285 = -qJD(5) * t77 + t132 * t115 - t191 * t93 + t196 * t94;
t76 = t124 * t196 - t125 * t191;
t284 = qJD(5) * t76 - t131 * t115 + t191 * t94 + t196 * t93;
t99 = t139 * t193 + t198 * t238;
t84 = pkin(8) * t175 + t99;
t229 = (t192 ^ 2 + t197 ^ 2) * t84;
t282 = m(3) * pkin(1);
t281 = m(5) + m(6);
t103 = t132 * t175;
t279 = t103 / 0.2e1;
t275 = pkin(2) * t198;
t274 = pkin(4) * t197;
t177 = sin(t190);
t272 = g(1) * t177;
t240 = qJD(4) * t197;
t174 = qJDD(3) + t184;
t244 = qJD(2) * t199;
t127 = (qJD(1) * t244 + qJDD(1) * t194) * pkin(1);
t50 = t198 * t127 + t139 * t242 + t292 * t193;
t40 = pkin(8) * t174 + t50;
t17 = -t192 * t40 - t240 * t84;
t270 = t17 * mrSges(5,3);
t169 = pkin(2) + t276;
t247 = t194 * t198;
t119 = pkin(1) * t247 + t193 * t169;
t113 = pkin(8) + t119;
t269 = -pkin(9) - t113;
t102 = t131 * t175;
t264 = mrSges(6,3) * t102;
t263 = mrSges(6,3) * t103;
t262 = Ifges(5,4) * t192;
t261 = Ifges(5,4) * t197;
t260 = Ifges(6,4) * t103;
t259 = Ifges(5,2) * t197;
t257 = qJD(4) * pkin(4);
t241 = qJD(4) * t192;
t16 = t197 * t40 - t241 * t84;
t256 = t16 * t197;
t255 = t17 * t192;
t233 = pkin(9) * t175 + t84;
t69 = t233 * t197;
t254 = t191 * t69;
t253 = t196 * t69;
t250 = t175 * t192;
t249 = t175 * t197;
t248 = t193 * t194;
t129 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t249;
t246 = t197 * t129;
t245 = t165 * pkin(3) + t164 * pkin(8);
t243 = qJD(3) * t193;
t239 = m(4) + t281;
t171 = pkin(4) * t241;
t179 = cos(t190);
t163 = pkin(2) * t179;
t235 = t163 + t245;
t167 = pkin(3) + t274;
t231 = qJD(4) * t269;
t144 = t266 - t267;
t108 = t144 * t175;
t228 = -t175 * mrSges(4,1) + t108;
t227 = t161 - t265;
t226 = -t164 * t201 + t165 * t167;
t118 = -pkin(1) * t248 + t169 * t198;
t64 = -mrSges(6,1) * t102 + mrSges(6,2) * t103;
t224 = -t228 - t64;
t112 = -pkin(3) - t118;
t223 = t163 + t226;
t68 = t233 * t192;
t221 = mrSges(6,1) * t176 + mrSges(6,2) * t178;
t220 = t259 + t262;
t218 = -t255 + t256;
t67 = -t68 + t257;
t21 = t196 * t67 - t254;
t22 = t191 * t67 + t253;
t89 = t269 * t192;
t90 = t113 * t197 + t181;
t56 = -t191 * t90 + t196 * t89;
t57 = t191 * t89 + t196 * t90;
t128 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t250;
t217 = t128 * t197 + t129 * t192;
t216 = t193 * t199 + t247;
t51 = -t193 * t127 - t139 * t243 + t292 * t198;
t214 = t192 * (Ifges(5,1) * t197 - t262);
t213 = t132 * qJD(5);
t212 = t131 * qJD(5);
t211 = t290 * t164 + t289 * t165;
t210 = t175 * mrSges(4,2) + t192 * t128 - t246;
t41 = -pkin(3) * t174 - t51;
t109 = t174 * t197 - t175 * t241;
t95 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t109;
t110 = t174 * t192 + t175 * t240;
t96 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t110;
t208 = -qJD(4) * t217 - t192 * t96 + t197 * t95;
t183 = qJDD(4) + qJDD(5);
t185 = qJD(4) + qJD(5);
t10 = qJDD(4) * pkin(4) - pkin(9) * t110 + t17;
t11 = pkin(9) * t109 + t16;
t3 = qJD(5) * t21 + t10 * t191 + t11 * t196;
t4 = -qJD(5) * t22 + t10 * t196 - t11 * t191;
t44 = t109 * t191 + t110 * t196 + t175 * t212;
t45 = t109 * t196 - t110 * t191 - t175 * t213;
t59 = Ifges(6,2) * t102 + Ifges(6,6) * t185 + t260;
t97 = Ifges(6,4) * t102;
t60 = Ifges(6,1) * t103 + Ifges(6,5) * t185 + t97;
t73 = -t167 * t175 - t98;
t207 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t21 * t264 + t59 * t279 - t73 * (mrSges(6,1) * t103 + mrSges(6,2) * t102) + Ifges(6,3) * t183 - t103 * (Ifges(6,1) * t102 - t260) / 0.2e1 + Ifges(6,6) * t45 + Ifges(6,5) * t44 - t185 * (Ifges(6,5) * t102 - Ifges(6,6) * t103) / 0.2e1 - (-Ifges(6,2) * t103 + t60 + t97) * t102 / 0.2e1;
t75 = t169 * t243 + (qJD(2) * t216 + t194 * t242) * pkin(1);
t206 = -t179 * mrSges(3,1) + mrSges(3,2) * t177 + t291 * t165 + t211;
t205 = (m(5) * pkin(3) + m(6) * t167 - t289 - t291) * t164 + (-m(5) * pkin(8) + m(6) * t201 + t290) * t165;
t204 = mrSges(3,2) * t179 + t205;
t100 = Ifges(5,6) * qJD(4) + t175 * t220;
t145 = Ifges(5,4) * t249;
t101 = Ifges(5,1) * t250 + Ifges(5,5) * qJD(4) + t145;
t27 = -pkin(4) * t109 + t41;
t81 = qJD(4) * t131 + t212;
t82 = -qJD(4) * t132 - t213;
t203 = t197 * (Ifges(5,4) * t110 + Ifges(5,2) * t109) / 0.2e1 + t185 * (Ifges(6,5) * t81 + Ifges(6,6) * t82) / 0.2e1 + Ifges(4,3) * t174 + t41 * t144 + t102 * (Ifges(6,4) * t81 + Ifges(6,2) * t82) / 0.2e1 + t82 * t59 / 0.2e1 + t73 * (-mrSges(6,1) * t82 + mrSges(6,2) * t81) + t81 * t60 / 0.2e1 + t51 * mrSges(4,1) + t110 * (Ifges(5,1) * t192 + t261) / 0.2e1 - t100 * t241 / 0.2e1 + (Ifges(5,1) * t110 + Ifges(5,4) * t109) * t278 + (Ifges(6,1) * t81 + Ifges(6,4) * t82) * t279 + t109 * t220 / 0.2e1 + mrSges(5,3) * t256 + (t101 + t175 * (-Ifges(5,2) * t192 + t261)) * t240 / 0.2e1 + (-t21 * t81 + t22 * t82) * mrSges(6,3) + (0.2e1 * Ifges(5,5) * t278 + Ifges(5,6) * t197) * qJDD(4) + (t214 * t175 / 0.2e1 + t288) * qJD(4) + (t27 * mrSges(6,2) - t4 * mrSges(6,3) + Ifges(6,1) * t44 + Ifges(6,4) * t45 + Ifges(6,5) * t183) * t132 + (-t27 * mrSges(6,1) + t3 * mrSges(6,3) + Ifges(6,4) * t44 + Ifges(6,2) * t45 + Ifges(6,6) * t183) * t131;
t202 = t126 * mrSges(3,1) + Ifges(3,3) * t184 + t203;
t200 = cos(qJ(1));
t195 = sin(qJ(1));
t182 = t200 * pkin(1);
t168 = -pkin(3) - t275;
t142 = -t167 - t275;
t138 = pkin(2) * t243 + t171;
t114 = t216 * t258;
t107 = t112 - t274;
t86 = mrSges(6,1) * t185 - t263;
t85 = -mrSges(6,2) * t185 + t264;
t74 = t169 * t242 + (-t194 * t243 + (t198 * t199 - t248) * qJD(2)) * pkin(1);
t72 = t171 + t75;
t71 = -mrSges(5,1) * t109 + mrSges(5,2) * t110;
t48 = -t192 * t74 + t197 * t231;
t47 = t192 * t231 + t197 * t74;
t36 = -mrSges(6,2) * t183 + mrSges(6,3) * t45;
t35 = mrSges(6,1) * t183 - mrSges(6,3) * t44;
t26 = -t196 * t68 - t254;
t25 = t191 * t68 - t253;
t9 = -mrSges(6,1) * t45 + mrSges(6,2) * t44;
t6 = -qJD(5) * t57 - t191 * t47 + t196 * t48;
t5 = qJD(5) * t56 + t191 * t48 + t196 * t47;
t1 = [t202 + t74 * t246 - t127 * mrSges(3,2) + t112 * t71 + t107 * t9 + t75 * t108 + (t126 * t199 + t127 * t194) * t282 + t5 * t85 + t6 * t86 + t72 * t64 + t56 * t35 + t57 * t36 + (-t119 * t174 - t175 * t74 - t50) * mrSges(4,2) + (mrSges(2,2) * t200 + (pkin(2) * t239 + mrSges(3,1)) * t177 + (mrSges(2,1) + (m(3) + t239) * pkin(1)) * t195 + t204) * g(1) + (mrSges(2,2) * t195 - m(6) * (t182 + t223) - m(5) * (t182 + t235) - m(4) * (t163 + t182) + (-mrSges(2,1) - t282) * t200 + t206) * g(2) + ((-t184 * t194 - t186 * t244) * mrSges(3,2) + (-qJD(2) * t186 * t194 + t184 * t199) * mrSges(3,1)) * pkin(1) + t208 * t113 + m(4) * (t118 * t51 + t119 * t50 + t74 * t99 - t75 * t98) + m(6) * (t107 * t27 + t21 * t6 + t22 * t5 + t3 * t57 + t4 * t56 + t72 * t73) + (t118 * t174 - t175 * t75) * mrSges(4,1) + (-t128 * t74 - t270) * t192 + Ifges(2,3) * qJDD(1) + m(5) * (t112 * t41 + t113 * t218 + t229 * t74 + t75 * t83); t202 - mrSges(5,3) * t255 + t168 * t71 + t142 * t9 + t138 * t64 + (mrSges(3,1) * t177 + t204) * g(1) + t208 * t166 + (-m(5) * t235 + t206) * g(2) - m(5) * (t114 * t83 + t115 * t229) + (t186 * t237 - t127) * mrSges(3,2) + t76 * t35 + t77 * t36 - t50 * mrSges(4,2) + t210 * t115 + t284 * t85 - m(4) * (-t114 * t98 + t115 * t99) + t224 * t114 + t285 * t86 + m(5) * (t16 * t251 - t166 * t255 + t168 * t41) + t186 * mrSges(3,1) * t238 + ((mrSges(4,1) * t198 - mrSges(4,2) * t193) * t174 + t281 * t272 + (-g(2) * t179 + t193 * t50 + t198 * t51 + t272) * m(4) + (t228 * t193 - t210 * t198 + m(5) * (t193 * t83 + t198 * t229) + m(4) * (-t193 * t98 + t198 * t99)) * qJD(3)) * pkin(2) + (-t223 * g(2) + t142 * t27 + t3 * t77 + t4 * t76 + (t138 - t114) * t73 + t284 * t22 + t285 * t21) * m(6); t203 - t167 * t9 + t91 * t35 + t92 * t36 + (g(2) * mrSges(5,2) * t165 + t64 * t257 - t270 + t98 * t128 + (-qJD(4) * t129 - t96) * pkin(8)) * t192 + t205 * g(1) + (-t98 * t129 + (-qJD(4) * t128 + t95) * pkin(8)) * t197 - pkin(3) * t71 + (t175 * t98 - t50) * mrSges(4,2) + t287 * t86 + t286 * t85 + (t165 * t265 + t211) * g(2) + t224 * t99 + (-t226 * g(2) - t167 * t27 + t3 * t92 + t4 * t91 + (t171 - t99) * t73 + t286 * t22 + t287 * t21) * m(6) + (-pkin(3) * t41 + pkin(8) * t218 - g(2) * t245 - t229 * t98 - t83 * t99) * m(5); t217 * t84 + t207 + t22 * t263 + Ifges(5,5) * t110 + Ifges(5,6) * t109 - t26 * t85 - t25 * t86 - t16 * mrSges(5,2) + t17 * mrSges(5,1) + (t100 * t278 + (t259 * t278 - t214 / 0.2e1) * t175 - (t145 + t101) * t197 / 0.2e1 - t288) * t175 + (t144 - t227) * g(3) + Ifges(5,3) * qJDD(4) + (-t64 * t250 + t191 * t36 + t196 * t35 + (t191 * t3 + t196 * t4 - g(3) * t197 + (-t175 * t73 + t283) * t192) * m(6) + (-t191 * t86 + t196 * t85 + (-t191 * t21 + t196 * t22) * m(6)) * qJD(5)) * pkin(4) - m(6) * (t21 * t25 + t22 * t26) + t283 * (t221 + t222); (t86 + t263) * t22 + t207 - g(3) * t227 - t21 * t85 + t283 * t221;];
tau = t1;
