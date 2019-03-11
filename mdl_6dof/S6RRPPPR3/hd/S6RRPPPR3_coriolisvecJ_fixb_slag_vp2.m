% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:53
% EndTime: 2019-03-09 08:14:03
% DurationCPUTime: 4.73s
% Computational Cost: add. (3922->504), mult. (8920->675), div. (0->0), fcn. (4971->6), ass. (0->219)
t300 = -mrSges(3,1) - mrSges(4,1);
t184 = sin(pkin(9));
t185 = cos(pkin(9));
t192 = cos(qJ(2));
t193 = -pkin(2) - pkin(3);
t181 = -qJ(5) + t193;
t190 = sin(qJ(2));
t201 = pkin(4) * t192 + t181 * t190;
t196 = qJD(2) * t201 + qJD(5) * t192;
t176 = t190 * qJD(3);
t231 = qJD(1) * qJD(2);
t223 = t192 * t231;
t244 = qJ(3) * t223 + qJD(1) * t176;
t46 = qJD(1) * t196 + t244;
t237 = qJD(4) * t190;
t238 = qJD(2) * t192;
t95 = pkin(7) * t223 + (-qJ(4) * t238 - t237) * qJD(1);
t85 = -qJD(2) * qJD(5) + t95;
t20 = -t184 * t85 + t185 * t46;
t247 = t185 * t190;
t204 = pkin(5) * t192 - pkin(8) * t247;
t200 = t204 * qJD(2);
t14 = qJD(1) * t200 + t20;
t21 = t184 * t46 + t185 * t85;
t224 = t190 * t231;
t219 = t184 * t224;
t15 = -pkin(8) * t219 + t21;
t189 = sin(qJ(6));
t191 = cos(qJ(6));
t233 = t184 * qJD(2);
t240 = qJD(1) * t192;
t124 = -t185 * t240 - t233;
t241 = qJD(1) * t190;
t208 = pkin(4) * t190 + qJ(5) * t192;
t121 = -qJD(1) * pkin(1) - pkin(2) * t240 - qJ(3) * t241;
t89 = pkin(3) * t240 + qJD(4) - t121;
t69 = qJD(1) * t208 + t89;
t164 = qJ(4) * t241;
t171 = pkin(7) * t241;
t232 = qJD(3) + t171;
t225 = -t164 + t232;
t87 = qJD(2) * t181 + t225;
t28 = -t184 * t87 + t185 * t69;
t16 = pkin(5) * t241 - pkin(8) * t124 + t28;
t123 = -t185 * qJD(2) + t184 * t240;
t29 = t184 * t69 + t185 * t87;
t22 = pkin(8) * t123 + t29;
t5 = t16 * t191 - t189 * t22;
t1 = qJD(6) * t5 + t14 * t189 + t15 * t191;
t246 = t191 * t185;
t125 = t184 * t189 - t246;
t126 = t184 * t191 + t185 * t189;
t6 = t16 * t189 + t191 * t22;
t2 = -qJD(6) * t6 + t14 * t191 - t15 * t189;
t234 = qJD(6) * t191;
t235 = qJD(6) * t189;
t112 = -t184 * t234 - t185 * t235;
t92 = t126 * t241;
t251 = t112 - t92;
t111 = t184 * t235 - t185 * t234;
t226 = t184 * t241;
t93 = -t189 * t226 + t241 * t246;
t252 = t111 - t93;
t299 = t1 * t125 + t126 * t2 - t251 * t6 - t252 * t5;
t172 = pkin(7) * t240;
t133 = -qJ(4) * t240 + t172;
t166 = qJ(3) * t240;
t75 = qJD(1) * t201 + t166;
t44 = -t133 * t184 + t185 * t75;
t30 = qJD(1) * t204 + t44;
t45 = t185 * t133 + t184 * t75;
t38 = -pkin(8) * t226 + t45;
t263 = pkin(8) - t181;
t128 = t263 * t184;
t129 = t263 * t185;
t73 = t128 * t189 - t129 * t191;
t297 = qJD(5) * t126 - qJD(6) * t73 + t189 * t38 - t191 * t30;
t294 = -qJD(1) / 0.2e1;
t293 = qJD(1) / 0.2e1;
t291 = qJD(2) / 0.2e1;
t72 = t128 * t191 + t129 * t189;
t290 = qJD(5) * t125 + qJD(6) * t72 - t189 * t30 - t191 * t38;
t102 = t125 * t192;
t222 = t191 * t123 - t124 * t189;
t130 = t171 - t164;
t289 = t130 + qJD(3);
t239 = qJD(2) * t190;
t198 = t125 * t239;
t36 = -qJD(1) * t198 + qJD(6) * t222;
t199 = t126 * t239;
t68 = t123 * t189 + t124 * t191;
t37 = -qJD(1) * t199 - qJD(6) * t68;
t288 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t36 + Ifges(7,6) * t37;
t227 = t193 * qJD(2);
t100 = t227 + t225;
t136 = -qJD(2) * pkin(2) + t232;
t161 = qJD(6) + t241;
t170 = Ifges(3,4) * t240;
t287 = (m(4) * t136 + (mrSges(4,2) + mrSges(3,3)) * t241 + t300 * qJD(2)) * pkin(7) + t136 * mrSges(4,2) + t28 * mrSges(6,1) + t5 * mrSges(7,1) + t89 * mrSges(5,1) + (-Ifges(5,4) * t192 - t190 * Ifges(5,2)) * t294 + (t190 * Ifges(4,1) - Ifges(4,5) * t192) * t293 + t170 / 0.2e1 + t161 * Ifges(7,3) + t68 * Ifges(7,5) + t222 * Ifges(7,6) + t124 * Ifges(6,5) + t123 * Ifges(6,6) - t100 * mrSges(5,3) - t121 * mrSges(4,3) - t29 * mrSges(6,2) - t6 * mrSges(7,2) + (Ifges(3,1) + Ifges(6,3)) * t241 / 0.2e1 + (Ifges(5,6) + Ifges(4,4) + Ifges(3,5)) * t291;
t183 = qJD(2) * qJ(3);
t114 = -t133 - t183;
t144 = t172 + t183;
t147 = mrSges(4,2) * t240 + qJD(2) * mrSges(4,3);
t169 = Ifges(4,5) * t241;
t205 = t184 * t29 + t185 * t28;
t258 = Ifges(6,2) * t184;
t260 = Ifges(6,4) * t185;
t209 = -t258 + t260;
t261 = Ifges(6,4) * t184;
t210 = Ifges(6,1) * t185 - t261;
t262 = mrSges(6,2) * t185;
t211 = mrSges(6,1) * t184 + t262;
t269 = t185 / 0.2e1;
t271 = -t184 / 0.2e1;
t99 = qJD(2) * pkin(4) + qJD(5) - t114;
t286 = -t205 * mrSges(6,3) - (m(4) * t144 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t240 + t147) * pkin(7) + t121 * mrSges(4,1) + t89 * mrSges(5,2) + t99 * t211 + Ifges(4,6) * t291 - Ifges(4,3) * t240 / 0.2e1 + t169 / 0.2e1 + (Ifges(3,4) * t190 + t192 * Ifges(3,2)) * t294 + (-t192 * Ifges(5,1) - Ifges(5,4) * t190) * t293 - t114 * mrSges(5,3) + t123 * t209 / 0.2e1 + t124 * t210 / 0.2e1 - t144 * mrSges(4,2) + (t124 * Ifges(6,4) + t123 * Ifges(6,2) + Ifges(6,6) * t241) * t271 + (t124 * Ifges(6,1) + t123 * Ifges(6,4) + Ifges(6,5) * t241) * t269 - (Ifges(3,6) + Ifges(5,5)) * qJD(2) / 0.2e1;
t285 = t36 / 0.2e1;
t284 = t37 / 0.2e1;
t283 = -t222 / 0.2e1;
t282 = t222 / 0.2e1;
t281 = -t68 / 0.2e1;
t280 = t68 / 0.2e1;
t279 = pkin(1) * mrSges(3,1);
t278 = pkin(1) * mrSges(3,2);
t101 = t126 * t192;
t277 = t101 / 0.2e1;
t276 = t102 / 0.2e1;
t275 = t125 / 0.2e1;
t274 = -t126 / 0.2e1;
t273 = -t161 / 0.2e1;
t272 = t161 / 0.2e1;
t270 = -t185 / 0.2e1;
t268 = Ifges(7,4) * t68;
t267 = pkin(8) * t192;
t264 = pkin(7) - qJ(4);
t188 = qJ(3) + pkin(4);
t108 = t238 * t264 - t237;
t242 = qJ(3) * t238 + t176;
t56 = t196 + t242;
t32 = t185 * t108 + t184 * t56;
t259 = Ifges(6,5) * t185;
t257 = Ifges(6,6) * t184;
t178 = t192 * qJ(4);
t149 = pkin(7) * t192 - t178;
t236 = qJD(4) * t192;
t202 = pkin(7) * t239 + t236;
t182 = qJD(2) * qJD(3);
t245 = -qJ(4) * t224 - t182;
t84 = qJD(1) * t202 + t245;
t254 = t149 * t84;
t148 = t264 * t190;
t137 = -t192 * pkin(2) - t190 * qJ(3) - pkin(1);
t122 = t192 * pkin(3) - t137;
t86 = t208 + t122;
t55 = t185 * t148 + t184 * t86;
t250 = -qJD(2) * mrSges(5,1) + mrSges(6,1) * t123 - mrSges(6,2) * t124 + mrSges(5,3) * t240;
t248 = t184 * t190;
t96 = mrSges(6,1) * t219 + t224 * t262;
t243 = mrSges(5,1) * t223 + mrSges(5,2) * t224;
t24 = -mrSges(7,1) * t222 + mrSges(7,2) * t68;
t230 = t24 - t250;
t229 = m(4) * pkin(7) + mrSges(4,2);
t228 = pkin(5) * t184 - pkin(7);
t11 = -t37 * mrSges(7,1) + t36 * mrSges(7,2);
t31 = -t108 * t184 + t185 * t56;
t54 = -t148 * t184 + t185 * t86;
t221 = -t147 - t230;
t220 = t190 * t227;
t218 = 0.3e1 / 0.2e1 * Ifges(3,4) + 0.3e1 / 0.2e1 * Ifges(5,4) - 0.3e1 / 0.2e1 * Ifges(4,5);
t217 = -Ifges(5,5) / 0.2e1 + Ifges(4,6) / 0.2e1 - Ifges(3,6) / 0.2e1;
t216 = Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t207 = t184 * t21 + t185 * t20;
t206 = -t184 * t20 + t185 * t21;
t41 = pkin(5) * t190 + t185 * t267 + t54;
t47 = t184 * t267 + t55;
t12 = -t189 * t47 + t191 * t41;
t13 = t189 * t41 + t191 * t47;
t203 = -t259 / 0.2e1 + t257 / 0.2e1;
t197 = t228 * t239 - t236;
t167 = qJ(4) * t239;
t159 = Ifges(7,3) * t223;
t154 = pkin(5) * t185 + t188;
t141 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t241;
t135 = -pkin(7) * t224 + t182;
t134 = (mrSges(5,1) * t190 - mrSges(5,2) * t192) * qJD(1);
t131 = (-mrSges(4,1) * t192 - mrSges(4,3) * t190) * qJD(1);
t110 = -t192 * t228 - t178;
t107 = pkin(2) * t239 - t242;
t106 = -t167 + t202;
t105 = t193 * t241 + t166;
t104 = (mrSges(6,1) * t192 - mrSges(6,3) * t247) * t231;
t103 = (-mrSges(6,2) * t192 - mrSges(6,3) * t248) * t231;
t97 = pkin(5) * t226 - t130;
t94 = pkin(2) * t224 - t244;
t91 = mrSges(6,1) * t241 - mrSges(6,3) * t124;
t90 = -mrSges(6,2) * t241 + mrSges(6,3) * t123;
t88 = t220 + t242;
t83 = t167 + t197;
t79 = (Ifges(6,5) * t192 + t190 * t210) * t231;
t78 = (Ifges(6,6) * t192 + t190 * t209) * t231;
t77 = qJD(1) * t220 + t244;
t64 = qJD(1) * t197 - t245;
t63 = -pkin(5) * t123 + t99;
t62 = Ifges(7,4) * t222;
t53 = -qJD(6) * t102 - t199;
t52 = qJD(6) * t101 - t198;
t49 = mrSges(7,1) * t161 - mrSges(7,3) * t68;
t48 = -mrSges(7,2) * t161 + mrSges(7,3) * t222;
t27 = -mrSges(7,2) * t223 + mrSges(7,3) * t37;
t26 = mrSges(7,1) * t223 - mrSges(7,3) * t36;
t25 = -pkin(8) * t190 * t233 + t32;
t23 = t200 + t31;
t19 = Ifges(7,1) * t68 + Ifges(7,5) * t161 + t62;
t18 = Ifges(7,2) * t222 + Ifges(7,6) * t161 + t268;
t8 = t36 * Ifges(7,1) + t37 * Ifges(7,4) + Ifges(7,5) * t223;
t7 = Ifges(7,4) * t36 + Ifges(7,2) * t37 + Ifges(7,6) * t223;
t4 = -qJD(6) * t13 - t189 * t25 + t191 * t23;
t3 = qJD(6) * t12 + t189 * t23 + t191 * t25;
t9 = [t108 * t141 + t149 * t96 + t107 * t131 + t88 * t134 + t110 * t11 + t64 * (-mrSges(7,1) * t101 + mrSges(7,2) * t102) + t55 * t103 + t54 * t104 + t83 * t24 + t32 * t90 + t31 * t91 + t63 * (-mrSges(7,1) * t53 + mrSges(7,2) * t52) + t52 * t19 / 0.2e1 + t53 * t18 / 0.2e1 + t3 * t48 + t4 * t49 + t12 * t26 + t13 * t27 + m(7) * (t1 * t13 + t110 * t64 + t12 * t2 + t3 * t6 + t4 * t5 + t63 * t83) + (t1 * t101 - t102 * t2 - t5 * t52 + t53 * t6) * mrSges(7,3) + (-t94 * mrSges(4,3) + t77 * mrSges(5,1) - t21 * mrSges(6,2) + t20 * mrSges(6,1) - t95 * mrSges(5,3) + t159 / 0.2e1 + (t217 * qJD(2) + ((0.3e1 / 0.2e1 * t259 - 0.3e1 / 0.2e1 * t257 - t218) * t190 - 0.2e1 * t279 + t137 * mrSges(4,1) + t149 * mrSges(5,3) + (0.3e1 / 0.2e1 * Ifges(6,3) + 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(4,3) - t185 ^ 2 * Ifges(6,1) / 0.2e1 + Ifges(7,3) / 0.2e1 + t229 * pkin(7) + (t260 - t258 / 0.2e1) * t184) * t192) * qJD(1) + t286) * qJD(2) + t288) * t190 + (t79 * t270 + t184 * t78 / 0.2e1 - t94 * mrSges(4,1) - t77 * mrSges(5,2) + (mrSges(5,3) + t211) * t84 + t229 * t135 + t207 * mrSges(6,3) + (t216 * qJD(2) + ((t203 + t218) * t192 + Ifges(7,5) * t276 + Ifges(7,6) * t277 - 0.2e1 * t278 - t137 * mrSges(4,3) - t148 * mrSges(5,3)) * qJD(1) + t287) * qJD(2)) * t192 + m(4) * (t107 * t121 + t137 * t94) + m(6) * (-t106 * t99 + t20 * t54 + t21 * t55 + t28 * t31 + t29 * t32 - t254) + m(5) * (t100 * t108 + t106 * t114 + t122 * t77 + t148 * t95 + t88 * t89 - t254) + t250 * t106 + t122 * t243 + (Ifges(7,5) * t52 + Ifges(7,6) * t53) * t272 + t8 * t276 + t7 * t277 + (Ifges(7,1) * t52 + Ifges(7,4) * t53) * t280 + (Ifges(7,4) * t52 + Ifges(7,2) * t53) * t282 + (Ifges(7,4) * t102 + Ifges(7,2) * t101) * t284 + (Ifges(7,1) * t102 + Ifges(7,4) * t101) * t285; t154 * t11 + t135 * mrSges(4,3) - t133 * t141 - t105 * t134 + t64 * (-mrSges(7,1) * t125 - mrSges(7,2) * t126) + t95 * mrSges(5,2) - t97 * t24 - t84 * mrSges(5,1) - t45 * t90 - t44 * t91 + t72 * t26 + t73 * t27 + (-t112 / 0.2e1 + t92 / 0.2e1) * t18 + t290 * t48 + t299 * mrSges(7,3) + (-t28 * t44 - t29 * t45 - t188 * t84 + t206 * t181 + (t184 * t28 - t185 * t29) * qJD(5) + t289 * t99) * m(6) + (-t84 * qJ(3) - t100 * t133 - t105 * t89 - t114 * t289 + t193 * t95) * m(5) + (t1 * t73 + t154 * t64 + t2 * t72 + (qJD(3) - t97) * t63 + t290 * t6 + t297 * t5) * m(7) + t297 * t49 + (-m(4) * t121 - t131) * (pkin(2) * t241 - t166) + (((pkin(7) * mrSges(3,2) + (-Ifges(6,2) * t185 - t261) * t271 + (-Ifges(6,1) * t184 - t260) * t269 + (-mrSges(4,2) + mrSges(5,3)) * qJ(3) + t217) * qJD(2) + (t279 + (Ifges(3,4) / 0.2e1 + Ifges(5,4) / 0.2e1 + t203) * t190) * qJD(1) - t169 / 0.2e1 - t286) * t190 + ((t278 + (-Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t192) * qJD(1) + (-Ifges(3,1) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t241 - t170 / 0.2e1 + (Ifges(7,5) * t274 + Ifges(7,6) * t275 + Ifges(6,5) * t271 + Ifges(6,6) * t270 - t193 * mrSges(5,3) - pkin(2) * mrSges(4,2) + (-m(4) * pkin(2) + t300) * pkin(7) + t216) * qJD(2) - t287) * t192) * qJD(1) + (t111 / 0.2e1 - t93 / 0.2e1) * t19 + (t181 * t103 - qJD(5) * t90 - t21 * mrSges(6,3) - t78 / 0.2e1 - t84 * mrSges(6,1)) * t185 + (-t181 * t104 + qJD(5) * t91 + t20 * mrSges(6,3) - t79 / 0.2e1 + t84 * mrSges(6,2)) * t184 + (mrSges(7,1) * t251 + mrSges(7,2) * t252) * t63 - t250 * t130 + t188 * t96 + (Ifges(7,5) * t111 - Ifges(7,6) * t112) * t272 + (Ifges(7,5) * t93 - Ifges(7,6) * t92) * t273 + t8 * t274 + t7 * t275 + (Ifges(7,1) * t111 - Ifges(7,4) * t112) * t280 + (Ifges(7,1) * t93 - Ifges(7,4) * t92) * t281 + (Ifges(7,4) * t111 - Ifges(7,2) * t112) * t282 + (Ifges(7,4) * t93 - Ifges(7,2) * t92) * t283 + (-Ifges(7,4) * t126 + Ifges(7,2) * t125) * t284 + (-Ifges(7,1) * t126 + Ifges(7,4) * t125) * t285 + m(4) * (qJ(3) * t135 + qJD(3) * t144) - t221 * qJD(3); t185 * t103 - t184 * t104 - t125 * t27 - t126 * t26 + t252 * t49 + t251 * t48 + t221 * qJD(2) + ((-mrSges(5,3) + t229) * t238 + (-t184 * t90 - t185 * t91 + t131 - t134) * t190) * qJD(1) - m(4) * (qJD(2) * t144 - t121 * t241) + (-qJD(2) * t63 - t299) * m(7) + (-qJD(2) * t99 - t205 * t241 + t206) * m(6) + (qJD(2) * t114 - t241 * t89 + t95) * m(5); t184 * t103 + t185 * t104 - t125 * t26 + t126 * t27 + t251 * t49 - t252 * t48 + m(5) * t77 + m(6) * t207 + m(7) * (t1 * t126 - t111 * t6 + t112 * t5 - t125 * t2) - m(7) * (t5 * t92 - t6 * t93) + (-m(6) * (-t247 * t29 + t248 * t28) + (m(5) * t100 - t184 * t91 + t185 * t90 + t141) * t190 + (-m(5) * t114 + m(6) * t99 + m(7) * t63 + t230) * t192) * qJD(1) + t243; -t123 * t90 + t124 * t91 - t222 * t48 + t68 * t49 + t11 + t96 + (-t222 * t6 + t5 * t68 + t64) * m(7) + (-t123 * t29 + t124 * t28 - t84) * m(6); t159 - t63 * (mrSges(7,1) * t68 + mrSges(7,2) * t222) + (Ifges(7,1) * t222 - t268) * t281 + t18 * t280 + (Ifges(7,5) * t222 - Ifges(7,6) * t68) * t273 - t5 * t48 + t6 * t49 + (t222 * t5 + t6 * t68) * mrSges(7,3) + (-Ifges(7,2) * t68 + t19 + t62) * t283 + t288;];
tauc  = t9(:);
