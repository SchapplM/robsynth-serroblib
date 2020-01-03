% Calculate vector of inverse dynamics joint torques for
% S5RPRRR10
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR10_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR10_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:44
% EndTime: 2019-12-31 19:10:11
% DurationCPUTime: 14.51s
% Computational Cost: add. (7719->582), mult. (18522->772), div. (0->0), fcn. (13796->14), ass. (0->263)
t378 = mrSges(5,3) + mrSges(6,3);
t254 = qJD(1) * qJD(2);
t181 = qJ(2) * qJDD(1) + t254;
t338 = m(6) * pkin(4);
t200 = sin(qJ(4));
t207 = -pkin(8) - pkin(7);
t244 = qJD(4) * t207;
t196 = sin(pkin(9));
t197 = cos(pkin(9));
t201 = sin(qJ(3));
t205 = cos(qJ(3));
t364 = -t196 * t201 + t205 * t197;
t158 = t364 * qJD(1);
t281 = t158 * t200;
t164 = t196 * t205 + t197 * t201;
t159 = t164 * qJD(1);
t114 = pkin(3) * t159 - pkin(7) * t158;
t299 = pkin(6) + qJ(2);
t173 = t299 * t196;
t165 = qJD(1) * t173;
t174 = t299 * t197;
t166 = qJD(1) * t174;
t119 = -t165 * t205 - t201 * t166;
t204 = cos(qJ(4));
t64 = t200 * t114 + t204 * t119;
t377 = pkin(8) * t281 + t200 * t244 - t64;
t280 = t158 * t204;
t63 = t204 * t114 - t119 * t200;
t376 = -pkin(4) * t159 + pkin(8) * t280 + t204 * t244 - t63;
t260 = t196 ^ 2 + t197 ^ 2;
t186 = pkin(4) * t204 + pkin(3);
t195 = qJ(4) + qJ(5);
t190 = sin(t195);
t191 = cos(t195);
t227 = -mrSges(5,1) * t204 + mrSges(5,2) * t200;
t375 = -m(5) * pkin(3) - m(6) * t186 - mrSges(6,1) * t191 + mrSges(6,2) * t190 + t227;
t374 = -m(5) * pkin(7) + m(6) * t207 - t378;
t257 = qJD(4) * t200;
t373 = t257 - t281;
t185 = pkin(2) * t197 + pkin(1);
t172 = -qJD(1) * t185 + qJD(2);
t348 = Ifges(4,5) * qJD(3);
t372 = -t172 * mrSges(4,2) + t119 * mrSges(4,3) - t348 / 0.2e1;
t161 = t164 * qJD(3);
t120 = -t165 * t201 + t166 * t205;
t203 = cos(qJ(5));
t199 = sin(qJ(5));
t130 = qJD(3) * t204 - t159 * t200;
t113 = qJD(3) * pkin(7) + t120;
t89 = -pkin(3) * t158 - pkin(7) * t159 + t172;
t50 = t113 * t204 + t200 * t89;
t43 = pkin(8) * t130 + t50;
t284 = t199 * t43;
t152 = qJD(4) - t158;
t131 = qJD(3) * t200 + t159 * t204;
t49 = -t113 * t200 + t204 * t89;
t42 = -pkin(8) * t131 + t49;
t38 = pkin(4) * t152 + t42;
t13 = t203 * t38 - t284;
t283 = t203 * t43;
t14 = t199 * t38 + t283;
t347 = Ifges(4,6) * qJD(3);
t371 = -t172 * mrSges(4,1) - t49 * mrSges(5,1) - t13 * mrSges(6,1) + t50 * mrSges(5,2) + t14 * mrSges(6,2) + t120 * mrSges(4,3) + t347 / 0.2e1;
t369 = t200 * t338;
t145 = qJD(5) + t152;
t234 = t203 * t130 - t131 * t199;
t366 = t152 * Ifges(5,3);
t367 = t130 * Ifges(5,6);
t81 = t130 * t199 + t131 * t203;
t368 = t131 * Ifges(5,5) + t81 * Ifges(6,5) + Ifges(6,6) * t234 + t145 * Ifges(6,3) + t366 + t367;
t148 = Ifges(4,4) * t158;
t365 = t158 * Ifges(4,2);
t350 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t130 - mrSges(5,2) * t131 - mrSges(4,3) * t159;
t160 = t364 * qJD(3);
t256 = qJD(4) * t204;
t216 = t200 * t160 + t164 * t256;
t363 = -m(5) - m(6) - m(4);
t240 = m(3) * qJ(2) + mrSges(3,3);
t362 = mrSges(2,2) - mrSges(4,3) - t240;
t117 = qJD(1) * t160 + qJDD(1) * t164;
t251 = qJDD(1) * t197;
t252 = qJDD(1) * t196;
t118 = -qJD(1) * t161 - t201 * t252 + t205 * t251;
t171 = -qJDD(1) * t185 + qJDD(2);
t56 = -pkin(3) * t118 - pkin(7) * t117 + t171;
t235 = pkin(6) * qJDD(1) + t181;
t146 = t235 * t196;
t147 = t235 * t197;
t258 = qJD(3) * t205;
t259 = qJD(3) * t201;
t61 = -t201 * t146 + t205 * t147 - t165 * t258 - t166 * t259;
t57 = qJDD(3) * pkin(7) + t61;
t11 = -t113 * t257 + t200 * t56 + t204 * t57 + t89 * t256;
t73 = -qJD(4) * t131 + qJDD(3) * t204 - t117 * t200;
t10 = pkin(8) * t73 + t11;
t111 = qJDD(4) - t118;
t12 = -qJD(4) * t50 - t200 * t57 + t204 * t56;
t72 = qJD(4) * t130 + qJDD(3) * t200 + t117 * t204;
t9 = pkin(4) * t111 - pkin(8) * t72 + t12;
t2 = qJD(5) * t13 + t10 * t203 + t199 * t9;
t3 = -qJD(5) * t14 - t10 * t199 + t203 * t9;
t359 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t358 = t12 * mrSges(5,1) - t11 * mrSges(5,2);
t194 = pkin(9) + qJ(3);
t188 = sin(t194);
t189 = cos(t194);
t229 = mrSges(4,1) * t189 - mrSges(4,2) * t188;
t230 = -mrSges(3,1) * t197 + mrSges(3,2) * t196;
t357 = m(3) * pkin(1) + t378 * t188 + mrSges(2,1) + t229 - t230;
t24 = qJD(5) * t234 + t199 * t73 + t203 * t72;
t337 = t24 / 0.2e1;
t25 = -qJD(5) * t81 - t199 * t72 + t203 * t73;
t336 = t25 / 0.2e1;
t329 = t72 / 0.2e1;
t328 = t73 / 0.2e1;
t108 = qJDD(5) + t111;
t323 = t108 / 0.2e1;
t322 = t111 / 0.2e1;
t176 = t207 * t200;
t177 = t207 * t204;
t128 = t176 * t203 + t177 * t199;
t355 = qJD(5) * t128 + t199 * t376 + t203 * t377;
t129 = t176 * t199 - t177 * t203;
t354 = -qJD(5) * t129 - t199 * t377 + t203 * t376;
t353 = mrSges(5,1) + t338;
t218 = t199 * t200 - t203 * t204;
t344 = qJD(4) + qJD(5);
t122 = t344 * t218;
t94 = t218 * t158;
t352 = -t122 + t94;
t168 = t199 * t204 + t200 * t203;
t123 = t344 * t168;
t93 = t168 * t158;
t351 = -t123 + t93;
t349 = pkin(4) * t373 - t120;
t104 = t218 * t164;
t116 = -pkin(3) * t364 - pkin(7) * t164 - t185;
t126 = -t173 * t201 + t174 * t205;
t121 = t204 * t126;
t71 = t200 * t116 + t121;
t346 = -t205 * t173 - t174 * t201;
t345 = t11 * t204 - t12 * t200;
t112 = -qJD(3) * pkin(3) - t119;
t74 = -pkin(4) * t130 + t112;
t342 = -t74 * mrSges(6,1) + t14 * mrSges(6,3);
t341 = t74 * mrSges(6,2) - t13 * mrSges(6,3);
t340 = Ifges(6,4) * t337 + Ifges(6,2) * t336 + Ifges(6,6) * t323;
t339 = Ifges(6,1) * t337 + Ifges(6,4) * t336 + Ifges(6,5) * t323;
t335 = Ifges(5,1) * t329 + Ifges(5,4) * t328 + Ifges(5,5) * t322;
t308 = Ifges(6,4) * t81;
t36 = Ifges(6,2) * t234 + Ifges(6,6) * t145 + t308;
t334 = -t36 / 0.2e1;
t333 = t36 / 0.2e1;
t77 = Ifges(6,4) * t234;
t37 = Ifges(6,1) * t81 + Ifges(6,5) * t145 + t77;
t332 = -t37 / 0.2e1;
t331 = t37 / 0.2e1;
t327 = -t234 / 0.2e1;
t326 = t234 / 0.2e1;
t325 = -t81 / 0.2e1;
t324 = t81 / 0.2e1;
t321 = -t130 / 0.2e1;
t320 = -t131 / 0.2e1;
t319 = t131 / 0.2e1;
t318 = -t145 / 0.2e1;
t317 = t145 / 0.2e1;
t316 = -t152 / 0.2e1;
t315 = -t158 / 0.2e1;
t313 = t159 / 0.2e1;
t310 = t204 / 0.2e1;
t307 = pkin(4) * t131;
t303 = g(3) * t188;
t298 = mrSges(5,3) * t130;
t297 = mrSges(5,3) * t131;
t296 = Ifges(4,4) * t159;
t295 = Ifges(5,4) * t200;
t294 = Ifges(5,4) * t204;
t288 = t131 * Ifges(5,4);
t279 = t164 * t200;
t278 = t164 * t204;
t202 = sin(qJ(1));
t276 = t190 * t202;
t206 = cos(qJ(1));
t275 = t190 * t206;
t274 = t191 * t202;
t273 = t191 * t206;
t268 = t200 * t202;
t267 = t200 * t206;
t266 = t202 * t204;
t265 = t204 * t160;
t264 = t204 * t206;
t141 = t189 * t276 + t273;
t142 = -t189 * t274 + t275;
t262 = -t141 * mrSges(6,1) + t142 * mrSges(6,2);
t143 = -t189 * t275 + t274;
t144 = t189 * t273 + t276;
t261 = t143 * mrSges(6,1) - t144 * mrSges(6,2);
t250 = Ifges(6,5) * t24 + Ifges(6,6) * t25 + Ifges(6,3) * t108;
t249 = Ifges(5,5) * t72 + Ifges(5,6) * t73 + Ifges(5,3) * t111;
t127 = Ifges(5,4) * t130;
t67 = Ifges(5,1) * t131 + Ifges(5,5) * t152 + t127;
t245 = t67 * t310;
t238 = -t257 / 0.2e1;
t237 = -t118 * mrSges(4,1) + t117 * mrSges(4,2);
t115 = pkin(3) * t161 - pkin(7) * t160;
t90 = qJD(2) * t364 + qJD(3) * t346;
t236 = t204 * t115 - t200 * t90;
t70 = t204 * t116 - t126 * t200;
t232 = pkin(3) * t189 + pkin(7) * t188;
t231 = -mrSges(3,1) * t251 + mrSges(3,2) * t252;
t226 = mrSges(5,1) * t200 + mrSges(5,2) * t204;
t225 = -mrSges(6,1) * t190 - mrSges(6,2) * t191;
t224 = Ifges(5,1) * t204 - t295;
t223 = -Ifges(5,2) * t200 + t294;
t222 = Ifges(5,5) * t204 - Ifges(5,6) * t200;
t47 = -pkin(4) * t364 - pkin(8) * t278 + t70;
t51 = -pkin(8) * t279 + t71;
t27 = -t199 * t51 + t203 * t47;
t28 = t199 * t47 + t203 * t51;
t87 = -mrSges(5,2) * t152 + t298;
t88 = mrSges(5,1) * t152 - t297;
t220 = -t200 * t88 + t204 * t87;
t219 = t186 * t189 - t188 * t207;
t217 = t250 + t359;
t62 = -t146 * t205 - t201 * t147 + t165 * t259 - t166 * t258;
t155 = -t189 * t267 + t266;
t153 = t189 * t268 + t264;
t215 = t164 * t257 - t265;
t214 = t112 * t226;
t31 = t200 * t115 + t116 * t256 - t126 * t257 + t204 * t90;
t58 = -qJDD(3) * pkin(3) - t62;
t91 = qJD(2) * t164 + qJD(3) * t126;
t187 = -qJDD(1) * pkin(1) + qJDD(2);
t156 = t189 * t264 + t268;
t154 = -t189 * t266 + t267;
t136 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t158;
t103 = t168 * t164;
t97 = t159 * Ifges(4,1) + t148 + t348;
t96 = t296 + t347 + t365;
t92 = pkin(4) * t279 - t346;
t66 = t130 * Ifges(5,2) + t152 * Ifges(5,6) + t288;
t60 = mrSges(6,1) * t145 - mrSges(6,3) * t81;
t59 = -mrSges(6,2) * t145 + mrSges(6,3) * t234;
t55 = pkin(4) * t216 + t91;
t46 = -mrSges(5,2) * t111 + mrSges(5,3) * t73;
t45 = mrSges(5,1) * t111 - mrSges(5,3) * t72;
t41 = -mrSges(6,1) * t234 + mrSges(6,2) * t81;
t40 = t104 * t344 - t168 * t160;
t39 = -t123 * t164 - t160 * t218;
t34 = -mrSges(5,1) * t73 + mrSges(5,2) * t72;
t33 = -pkin(4) * t73 + t58;
t32 = -qJD(4) * t71 + t236;
t29 = t72 * Ifges(5,4) + t73 * Ifges(5,2) + t111 * Ifges(5,6);
t26 = -pkin(8) * t216 + t31;
t19 = -pkin(8) * t265 + pkin(4) * t161 + (-t121 + (pkin(8) * t164 - t116) * t200) * qJD(4) + t236;
t18 = t203 * t42 - t284;
t17 = -t199 * t42 - t283;
t16 = -mrSges(6,2) * t108 + mrSges(6,3) * t25;
t15 = mrSges(6,1) * t108 - mrSges(6,3) * t24;
t8 = -mrSges(6,1) * t25 + mrSges(6,2) * t24;
t5 = -qJD(5) * t28 + t19 * t203 - t199 * t26;
t4 = qJD(5) * t27 + t19 * t199 + t203 * t26;
t1 = [m(4) * (t120 * t90 + t126 * t61 - t171 * t185) + (-t11 * t279 - t12 * t278 + t215 * t49 - t216 * t50) * mrSges(5,3) + 0.2e1 * t260 * t181 * mrSges(3,3) + (t171 * mrSges(4,2) - t62 * mrSges(4,3) + Ifges(4,1) * t117 + Ifges(4,4) * t118 + Ifges(4,5) * qJDD(3) + t222 * t322 + t223 * t328 + t224 * t329 + t226 * t58 + t238 * t67) * t164 - (-m(4) * t62 + m(5) * t58 - qJDD(3) * mrSges(4,1) + mrSges(4,3) * t117 + t34) * t346 + t33 * (mrSges(6,1) * t103 - mrSges(6,2) * t104) - (-Ifges(4,1) * t313 - t148 / 0.2e1 - t97 / 0.2e1 - t245 + t372) * t160 + m(6) * (t13 * t5 + t14 * t4 + t2 * t28 + t27 * t3 + t33 * t92 + t55 * t74) + (Ifges(6,4) * t39 + Ifges(6,2) * t40) * t326 - (t250 + t249) * t364 / 0.2e1 - (mrSges(4,1) * t171 - mrSges(4,3) * t61 - Ifges(4,4) * t117 + Ifges(5,5) * t329 + Ifges(6,5) * t337 - Ifges(4,2) * t118 - Ifges(4,6) * qJDD(3) + Ifges(5,6) * t328 + Ifges(6,6) * t336 + Ifges(5,3) * t322 + Ifges(6,3) * t323 + t358 + t359) * t364 + t27 * t15 + t28 * t16 + (-t103 * t2 + t104 * t3 - t13 * t39 + t14 * t40) * mrSges(6,3) + m(5) * (t11 * t71 + t12 * t70 + t31 * t50 + t32 * t49) + t90 * t136 + (-m(4) * t119 + m(5) * t112 - t350) * t91 + t126 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t118) - t216 * t66 / 0.2e1 + t31 * t87 + t32 * t88 + t92 * t8 + t74 * (-mrSges(6,1) * t40 + mrSges(6,2) * t39) + t70 * t45 + t71 * t46 + t4 * t59 + t5 * t60 + t55 * t41 + t130 * (-Ifges(5,4) * t215 - Ifges(5,2) * t216) / 0.2e1 + (Ifges(6,5) * t39 + Ifges(6,6) * t40) * t317 + (Ifges(3,4) * t196 + Ifges(3,2) * t197) * t251 + (Ifges(3,1) * t196 + Ifges(3,4) * t197) * t252 + (-Ifges(6,5) * t104 - Ifges(6,6) * t103) * t323 + t152 * (-Ifges(5,5) * t215 - Ifges(5,6) * t216) / 0.2e1 + (-Ifges(5,1) * t215 - Ifges(5,4) * t216) * t319 + (Ifges(6,1) * t39 + Ifges(6,4) * t40) * t324 + (-t268 * t338 - t156 * mrSges(5,1) - t144 * mrSges(6,1) - t155 * mrSges(5,2) - t143 * mrSges(6,2) + t363 * (t206 * t185 + t202 * t299) + t362 * t202 + (-m(5) * t232 - m(6) * t219 - t357) * t206) * g(2) + (-t154 * mrSges(5,1) - t142 * mrSges(6,1) - t153 * mrSges(5,2) - t141 * mrSges(6,2) + (t299 * t363 + t362 - t369) * t206 + (-m(6) * (-t185 - t219) - m(5) * (-t185 - t232) + m(4) * t185 + t357) * t202) * g(1) + (t368 / 0.2e1 - Ifges(4,4) * t313 - t96 / 0.2e1 - t365 / 0.2e1 + t367 / 0.2e1 + t366 / 0.2e1 + Ifges(6,3) * t317 + Ifges(5,5) * t319 + Ifges(6,5) * t324 + Ifges(6,6) * t326 - t371) * t161 + m(3) * (-pkin(1) * t187 + (t181 + t254) * qJ(2) * t260) + Ifges(2,3) * qJDD(1) + (-Ifges(6,1) * t104 - Ifges(6,4) * t103) * t337 + (-Ifges(6,4) * t104 - Ifges(6,2) * t103) * t336 + t112 * (mrSges(5,1) * t216 - mrSges(5,2) * t215) + t39 * t331 + t40 * t333 + t278 * t335 - t104 * t339 - t103 * t340 + t187 * t230 - pkin(1) * t231 - t185 * t237 - t29 * t279 / 0.2e1; t237 + t231 + (-t41 + t350) * t159 + m(3) * t187 + t351 * t60 + t352 * t59 - t218 * t15 + t168 * t16 + t200 * t46 + t204 * t45 + (-t136 - t220) * t158 + t220 * qJD(4) + (-g(1) * t202 + g(2) * t206) * (m(3) - t363) - t240 * t260 * qJD(1) ^ 2 + (t13 * t351 + t14 * t352 - t159 * t74 + t168 * t2 - t218 * t3) * m(6) + (t11 * t200 - t112 * t159 + t12 * t204 + t152 * (-t49 * t200 + t50 * t204)) * m(5) + (t119 * t159 - t120 * t158 + t171) * m(4); (t214 + t245) * qJD(4) + (-pkin(3) * t58 - t112 * t120 - t49 * t63 - t50 * t64) * m(5) + t349 * t41 + t350 * t120 + (-t13 * t94 + t14 * t93 - t168 * t3 - t2 * t218) * mrSges(6,3) + t33 * (mrSges(6,1) * t218 + mrSges(6,2) * t168) + (Ifges(6,5) * t168 - Ifges(6,6) * t218) * t323 + (Ifges(6,4) * t168 - Ifges(6,2) * t218) * t336 + (Ifges(6,1) * t168 - Ifges(6,4) * t218) * t337 - t218 * t340 + (t148 + t97) * t315 + (t130 * t223 + t131 * t224 + t152 * t222) * qJD(4) / 0.2e1 - t74 * (mrSges(6,1) * t93 - mrSges(6,2) * t94) + (-Ifges(6,5) * t94 - Ifges(6,6) * t93) * t318 + (-Ifges(6,1) * t94 - Ifges(6,4) * t93) * t325 + (-Ifges(6,4) * t94 - Ifges(6,2) * t93) * t327 - (Ifges(6,1) * t324 + Ifges(6,4) * t326 + Ifges(6,5) * t317 + t331 + t341) * t122 - (Ifges(6,4) * t324 + Ifges(6,2) * t326 + Ifges(6,6) * t317 + t333 + t342) * t123 + (t222 * t316 + t223 * t321 + t224 * t320 - t214 + t372) * t158 + (-t373 * t50 + (-t256 + t280) * t49 + t345) * mrSges(5,3) + (g(1) * t206 + g(2) * t202) * ((mrSges(4,2) + t374) * t189 + (mrSges(4,1) - t375) * t188) + (t188 * t374 + t189 * t375 - t229) * g(3) + t29 * t310 - pkin(3) * t34 - (Ifges(4,1) * t158 - t296 + t368) * t159 / 0.2e1 - t119 * t136 + t128 * t15 + t129 * t16 + Ifges(4,5) * t117 + Ifges(4,6) * t118 - t64 * t87 - t63 * t88 - t61 * mrSges(4,2) + t62 * mrSges(4,1) + (t238 + t281 / 0.2e1) * t66 + (m(5) * ((-t200 * t50 - t204 * t49) * qJD(4) + t345) - t88 * t256 - t87 * t257 - t200 * t45 + t204 * t46) * pkin(7) + (Ifges(5,5) * t320 + Ifges(6,5) * t325 - Ifges(4,2) * t315 + Ifges(5,6) * t321 + Ifges(6,6) * t327 + Ifges(5,3) * t316 + Ifges(6,3) * t318 + t371) * t159 + t96 * t313 + t354 * t60 + t355 * t59 + (t128 * t3 + t129 * t2 + t13 * t354 + t14 * t355 - t186 * t33 + t349 * t74) * m(6) - t186 * t8 + (Ifges(5,5) * t200 + Ifges(5,6) * t204) * t322 + (Ifges(5,2) * t204 + t295) * t328 + (Ifges(5,1) * t200 + t294) * t329 - t94 * t332 - t93 * t334 + t200 * t335 + t168 * t339 + t58 * t227 + Ifges(4,3) * qJDD(3) - t67 * t280 / 0.2e1; (-mrSges(5,2) * t154 + t153 * t353 - t262) * g(2) + (mrSges(5,2) * t156 - t155 * t353 - t261) * g(1) + t358 + (Ifges(6,1) * t325 + Ifges(6,4) * t327 + Ifges(6,5) * t318 + t332 - t341) * t234 - (Ifges(6,4) * t325 + Ifges(6,2) * t327 + Ifges(6,6) * t318 + t334 - t342) * t81 + (t298 - t87) * t49 + (-t225 + t226 + t369) * t303 + t217 - t112 * (mrSges(5,1) * t131 + mrSges(5,2) * t130) - t18 * t59 - t17 * t60 + (t297 + t88) * t50 + t249 + (-Ifges(5,2) * t131 + t127 + t67) * t321 + ((-t199 * t60 + t203 * t59) * qJD(5) + t203 * t15 + t199 * t16) * pkin(4) - t41 * t307 - m(6) * (t13 * t17 + t14 * t18 + t307 * t74) + (Ifges(5,5) * t130 - Ifges(5,6) * t131) * t316 + t66 * t319 + (Ifges(5,1) * t130 - t288) * t320 + (t199 * t2 + t203 * t3 + (-t13 * t199 + t14 * t203) * qJD(5)) * t338; -t74 * (mrSges(6,1) * t81 + mrSges(6,2) * t234) + (Ifges(6,1) * t234 - t308) * t325 + t36 * t324 + (Ifges(6,5) * t234 - Ifges(6,6) * t81) * t318 - t13 * t59 + t14 * t60 - g(1) * t261 - g(2) * t262 - t225 * t303 + (t13 * t234 + t14 * t81) * mrSges(6,3) + t217 + (-Ifges(6,2) * t81 + t37 + t77) * t327;];
tau = t1;
