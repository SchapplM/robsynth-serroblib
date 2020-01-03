% Calculate vector of inverse dynamics joint torques for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:16
% EndTime: 2020-01-03 11:49:40
% DurationCPUTime: 10.82s
% Computational Cost: add. (3827->428), mult. (9236->568), div. (0->0), fcn. (6284->10), ass. (0->228)
t338 = Ifges(5,4) + Ifges(6,4);
t339 = Ifges(5,1) + Ifges(6,1);
t337 = Ifges(5,5) + Ifges(6,5);
t336 = Ifges(5,2) + Ifges(6,2);
t335 = Ifges(5,6) + Ifges(6,6);
t196 = sin(pkin(8));
t199 = sin(qJ(3));
t202 = cos(qJ(3));
t254 = qJD(1) * qJD(3);
t136 = (qJDD(1) * t202 - t199 * t254) * t196;
t137 = (-qJDD(1) * t199 - t202 * t254) * t196;
t198 = sin(qJ(4));
t201 = cos(qJ(4));
t219 = t198 * t199 - t201 * t202;
t215 = t219 * qJD(4);
t261 = qJD(1) * t196;
t52 = -t136 * t198 + t137 * t201 + t215 * t261;
t360 = t335 * t52;
t152 = t198 * t202 + t199 * t201;
t216 = qJD(1) * t152;
t109 = t196 * t216;
t51 = -qJD(4) * t109 + t136 * t201 + t137 * t198;
t359 = t337 * t51;
t197 = cos(pkin(8));
t251 = qJDD(1) * t197;
t173 = qJDD(3) - t251;
t163 = qJDD(4) + t173;
t358 = (Ifges(5,3) + Ifges(6,3)) * t163;
t253 = qJDD(1) * qJ(2);
t255 = qJD(1) * qJD(2);
t170 = t253 + t255;
t357 = t170 + t255;
t341 = -mrSges(5,1) - mrSges(6,1);
t340 = mrSges(5,2) + mrSges(6,2);
t355 = t338 * t109;
t238 = t202 * t261;
t239 = t199 * t261;
t111 = -t198 * t239 + t201 * t238;
t354 = t338 * t111;
t260 = qJD(1) * t197;
t175 = qJD(3) - t260;
t164 = qJD(4) + t175;
t353 = -t336 * t109 + t335 * t164 + t354;
t352 = t339 * t111 + t337 * t164 - t355;
t351 = Ifges(4,5) * t136;
t350 = Ifges(4,6) * t137;
t349 = Ifges(4,3) * t173;
t193 = t196 ^ 2;
t348 = t357 * t193;
t226 = pkin(2) * t197 + pkin(6) * t196;
t160 = -pkin(1) - t226;
t135 = qJDD(1) * t160 + qJDD(2);
t278 = t197 * t199;
t138 = qJD(1) * t160 + qJD(2);
t241 = qJ(2) * t260;
t97 = t138 * t199 + t202 * t241;
t50 = -qJD(3) * t97 + t202 * t135 - t170 * t278;
t21 = pkin(3) * t173 - pkin(7) * t136 + t50;
t123 = t202 * t138;
t243 = qJ(2) * t278;
t279 = t196 * t202;
t250 = pkin(7) * t279;
t214 = -t243 - t250;
t79 = qJD(1) * t214 + t123;
t62 = pkin(3) * t175 + t79;
t80 = -pkin(7) * t239 + t97;
t74 = t201 * t80;
t24 = t198 * t62 + t74;
t227 = qJD(3) * t243;
t258 = qJD(3) * t202;
t276 = t197 * t202;
t49 = -qJD(1) * t227 + t199 * t135 + t138 * t258 + t170 * t276;
t27 = pkin(7) * t137 + t49;
t6 = -qJD(4) * t24 - t198 * t27 + t201 * t21;
t2 = pkin(4) * t163 - qJ(5) * t51 - qJD(5) * t111 + t6;
t256 = qJD(4) * t201;
t257 = qJD(4) * t198;
t5 = t198 * t21 + t201 * t27 + t62 * t256 - t257 * t80;
t3 = qJ(5) * t52 - qJD(5) * t109 + t5;
t347 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2);
t346 = -t50 * mrSges(4,1) + t49 * mrSges(4,2);
t316 = m(5) * pkin(3);
t36 = -mrSges(6,2) * t163 + mrSges(6,3) * t52;
t37 = -mrSges(5,2) * t163 + mrSges(5,3) * t52;
t334 = t36 + t37;
t291 = mrSges(6,3) * t109;
t83 = -mrSges(6,2) * t164 - t291;
t292 = mrSges(5,3) * t109;
t84 = -mrSges(5,2) * t164 - t292;
t333 = t83 + t84;
t290 = mrSges(6,3) * t111;
t85 = mrSges(6,1) * t164 - t290;
t287 = t111 * mrSges(5,3);
t86 = mrSges(5,1) * t164 - t287;
t332 = t85 + t86;
t331 = -mrSges(4,1) - t316;
t103 = t111 * qJ(5);
t72 = t198 * t80;
t23 = t201 * t62 - t72;
t13 = -t103 + t23;
t323 = -qJD(3) - qJD(4);
t330 = (t260 + t323) * t219;
t95 = t323 * t152;
t329 = t197 * t216 + t95;
t174 = qJ(2) * t276;
t107 = t199 * t160 + t174;
t194 = t197 ^ 2;
t328 = t194 + t193;
t327 = t358 + t359 + t360;
t195 = qJ(3) + qJ(4);
t186 = sin(t195);
t203 = cos(qJ(1));
t269 = t203 * t186;
t187 = cos(t195);
t200 = sin(qJ(1));
t272 = t200 * t187;
t126 = t197 * t269 - t272;
t275 = t197 * t203;
t127 = t186 * t200 + t187 * t275;
t326 = t126 * t341 - t340 * t127;
t277 = t197 * t200;
t124 = -t186 * t277 - t187 * t203;
t125 = t197 * t272 - t269;
t325 = t124 * t341 + t340 * t125;
t324 = mrSges(5,1) * t186 + t187 * t340;
t322 = -m(5) - m(6) - m(3) - m(4);
t321 = t194 * t357;
t299 = pkin(3) * t199;
t156 = pkin(4) * t186 + t299;
t320 = m(6) * t156 - mrSges(2,2) + mrSges(3,3);
t288 = Ifges(4,4) * t202;
t289 = Ifges(4,4) * t199;
t319 = (qJ(2) * (mrSges(4,1) * t202 - mrSges(4,2) * t199) - t199 * (-Ifges(4,2) * t202 - t289) / 0.2e1 + t202 * (-Ifges(4,1) * t199 - t288) / 0.2e1) * t193;
t190 = t202 * pkin(3);
t157 = pkin(4) * t187 + t190;
t204 = -pkin(7) - pkin(6);
t224 = mrSges(3,1) * t197 - mrSges(3,2) * t196;
t318 = -m(4) * t226 - mrSges(2,1) - t224 + (-m(5) * (t190 + pkin(2)) - m(6) * (pkin(2) + t157)) * t197 + (m(5) * t204 - mrSges(5,3) - m(6) * (qJ(5) - t204) - mrSges(6,3) - mrSges(4,3)) * t196;
t217 = (-t199 * Ifges(4,2) + t288) * t196;
t218 = (t202 * Ifges(4,1) - t289) * t196;
t317 = t202 * (t175 * Ifges(4,6) + qJD(1) * t217) + t199 * (t175 * Ifges(4,5) + qJD(1) * t218);
t205 = qJD(1) ^ 2;
t315 = m(6) * pkin(4);
t307 = t111 / 0.2e1;
t298 = pkin(3) * t201;
t297 = pkin(4) * t111;
t296 = g(1) * t196;
t29 = t201 * t79 - t72;
t149 = t202 * t160;
t87 = -t250 + t149 + (-qJ(2) * t199 - pkin(3)) * t197;
t280 = t196 * t199;
t98 = -pkin(7) * t280 + t107;
t39 = t198 * t87 + t201 * t98;
t283 = qJ(2) * t205;
t282 = qJ(5) * t109;
t185 = t196 * qJ(2);
t155 = t196 * t170;
t274 = t199 * t200;
t273 = t199 * t203;
t271 = t200 * t202;
t270 = t202 * t203;
t259 = qJD(2) * t197;
t264 = t160 * t258 + t202 * t259;
t263 = t348 * qJ(2);
t139 = pkin(3) * t239 + qJ(2) * t261;
t184 = t196 * qJD(2);
t236 = t196 * t258;
t146 = pkin(3) * t236 + t184;
t150 = pkin(3) * t280 + t185;
t252 = qJDD(1) * t196;
t247 = mrSges(4,3) * t280;
t242 = t349 + t350 + t351;
t237 = t199 * t259;
t235 = -t52 * mrSges(6,1) + t51 * mrSges(6,2);
t93 = -pkin(3) * t137 + t155;
t28 = -t198 * t79 - t74;
t38 = -t198 * t98 + t201 * t87;
t230 = pkin(3) * t238;
t229 = mrSges(4,3) * t239;
t228 = mrSges(4,3) * t238;
t225 = -mrSges(3,1) * t251 + mrSges(3,2) * t252;
t223 = mrSges(4,1) * t199 + mrSges(4,2) * t202;
t132 = -mrSges(4,2) * t175 - t229;
t133 = mrSges(4,1) * t175 - t228;
t220 = t132 * t202 - t133 * t199;
t142 = t197 * t273 - t271;
t140 = -t197 * t274 - t270;
t75 = qJD(3) * t214 + t264;
t76 = -t237 + (-t174 + (pkin(7) * t196 - t160) * t199) * qJD(3);
t9 = t198 * t76 + t201 * t75 + t87 * t256 - t257 * t98;
t212 = t175 * t196 * (-Ifges(4,5) * t199 - Ifges(4,6) * t202);
t10 = -qJD(4) * t39 - t198 * t75 + t201 * t76;
t12 = pkin(4) * t164 + t13;
t14 = t24 - t282;
t77 = pkin(4) * t109 + qJD(5) + t139;
t206 = -t12 * t291 - t23 * t292 + t14 * t290 + t24 * t287 - t77 * (mrSges(6,1) * t111 - mrSges(6,2) * t109) - t139 * (mrSges(5,1) * t111 - mrSges(5,2) * t109) + t327 - (-t339 * t109 - t354) * t111 / 0.2e1 + t353 * t307 - (-t109 * t337 - t111 * t335) * t164 / 0.2e1 + (-t336 * t111 + t352 - t355) * t109 / 0.2e1 + t347;
t183 = -qJDD(1) * pkin(1) + qJDD(2);
t181 = pkin(4) + t298;
t177 = t193 * t283;
t144 = t223 * t196;
t143 = t197 * t270 + t274;
t141 = t197 * t271 - t273;
t134 = t223 * t261;
t131 = t219 * t196;
t130 = t152 * t196;
t106 = t149 - t243;
t100 = -mrSges(4,2) * t173 + mrSges(4,3) * t137;
t99 = mrSges(4,1) * t173 - mrSges(4,3) * t136;
t96 = -t199 * t241 + t123;
t92 = t230 + t297;
t90 = -qJD(3) * t107 - t237;
t89 = -t227 + t264;
t88 = pkin(4) * t130 + t150;
t71 = (qJD(3) * t219 + t215) * t196;
t70 = t95 * t196;
t66 = mrSges(5,1) * t109 + mrSges(5,2) * t111;
t65 = mrSges(6,1) * t109 + mrSges(6,2) * t111;
t53 = -pkin(4) * t71 + t146;
t35 = mrSges(5,1) * t163 - mrSges(5,3) * t51;
t34 = mrSges(6,1) * t163 - mrSges(6,3) * t51;
t25 = -qJ(5) * t130 + t39;
t22 = -pkin(4) * t197 + qJ(5) * t131 + t38;
t18 = -pkin(4) * t52 + qJDD(5) + t93;
t17 = -t103 + t29;
t16 = t28 + t282;
t8 = -qJ(5) * t70 + qJD(5) * t131 + t10;
t7 = qJ(5) * t71 - qJD(5) * t130 + t9;
t1 = [t352 * t70 / 0.2e1 + t353 * t71 / 0.2e1 + (t253 * t328 + t321 + t348) * mrSges(3,3) + (-t274 * t316 - t143 * mrSges(4,1) + t142 * mrSges(4,2) + t322 * (t203 * pkin(1) + t200 * qJ(2)) - t320 * t200 + t341 * t127 + t340 * t126 + t318 * t203) * g(2) + (t273 * t316 - t141 * mrSges(4,1) - t140 * mrSges(4,2) + t322 * (t200 * pkin(1) - qJ(2) * t203) + t320 * t203 + t341 * t125 - t340 * t124 + t318 * t200) * g(3) + (-mrSges(4,1) * t137 + mrSges(4,2) * t136) * t185 + t134 * t184 + (t335 * t71 + t337 * t70) * t164 / 0.2e1 - (t336 * t71 + t338 * t70) * t109 / 0.2e1 + (t338 * t71 + t339 * t70) * t307 + (-t12 * t70 + t14 * t71) * mrSges(6,3) + (-t23 * t70 + t24 * t71) * mrSges(5,3) + (-t236 * t97 - t279 * t50) * mrSges(4,3) + t319 * t254 + m(3) * (-pkin(1) * t183 + qJ(2) * t321 + t263) + t137 * t217 / 0.2e1 + (-mrSges(5,2) * t93 - mrSges(6,2) * t18 + mrSges(5,3) * t6 + mrSges(6,3) * t2 - t163 * t337 - t338 * t52 - t339 * t51) * t131 + t146 * t66 + t150 * (-mrSges(5,1) * t52 + mrSges(5,2) * t51) + t139 * (-mrSges(5,1) * t71 + mrSges(5,2) * t70) + t89 * t132 + t90 * t133 + t106 * t99 + t107 * t100 + t7 * t83 + t9 * t84 + t8 * t85 + t10 * t86 + t77 * (-mrSges(6,1) * t71 + mrSges(6,2) * t70) + t53 * t65 + t22 * t34 + t25 * t36 + t38 * t35 + t39 * t37 - (Ifges(4,4) * t136 + Ifges(4,2) * t137 + Ifges(4,6) * t173) * t280 / 0.2e1 + t144 * t155 + (Ifges(4,1) * t136 + Ifges(4,4) * t137 + Ifges(4,5) * t173) * t279 / 0.2e1 + m(4) * (t106 * t50 + t107 * t49 + t89 * t97 + t90 * t96 + t263) - t183 * t224 - pkin(1) * t225 + m(5) * (t10 * t23 + t139 * t146 + t150 * t93 + t24 * t9 + t38 * t6 + t39 * t5) + m(6) * (t12 * t8 + t14 * t7 + t18 * t88 + t2 * t22 + t25 * t3 + t53 * t77) + t88 * t235 + (mrSges(5,1) * t93 + mrSges(6,1) * t18 - mrSges(5,3) * t5 - mrSges(6,3) * t3 - t163 * t335 - t336 * t52 - t338 * t51) * t130 + (Ifges(3,2) * t251 + Ifges(3,4) * t252 - t358 / 0.2e1 - t360 / 0.2e1 - t359 / 0.2e1 - t242 / 0.2e1 - t327 / 0.2e1 - t350 / 0.2e1 - t351 / 0.2e1 - t349 / 0.2e1 + t346 - t347) * t197 + t136 * t218 / 0.2e1 + (Ifges(3,4) * t251 + Ifges(3,1) * t252 - t317 * qJD(3) / 0.2e1 + t173 * (Ifges(4,5) * t202 - Ifges(4,6) * t199) / 0.2e1) * t196 + (t212 / 0.2e1 + t96 * t247) * qJD(3) - t49 * t247 + Ifges(2,3) * qJDD(1); t220 * qJD(3) + (-t220 * t197 + (-t134 - t65 - t66) * t196) * qJD(1) + t199 * t100 - t328 * t205 * mrSges(3,3) + t225 + t334 * t152 - (t34 + t35) * t219 + t202 * t99 + t330 * t333 + t329 * t332 - (g(2) * t203 + g(3) * t200) * t322 + (t12 * t329 + t14 * t330 + t152 * t3 - t2 * t219 - t77 * t261) * m(6) + (-t139 * t261 + t152 * t5 - t219 * t6 + t23 * t329 + t24 * t330) * m(5) + (t199 * t49 + t202 * t50 - t177 + t175 * (-t199 * t96 + t202 * t97)) * m(4) + (-t194 * t283 - t177 + t183) * m(3); -qJD(1) * t212 / 0.2e1 + (t133 + t228) * t97 + (t198 * t5 + t201 * t6 + (-t198 * t23 + t201 * t24) * qJD(4)) * t316 + t35 * t298 + (-t229 - t132) * t96 - t346 + (mrSges(4,2) * t141 - m(6) * (-t156 * t277 - t157 * t203) + t331 * t140 + t325) * g(2) + (-mrSges(4,2) * t143 - m(6) * (t156 * t275 - t157 * t200) + t331 * t142 + t326) * g(3) + (m(5) * t299 + mrSges(6,1) * t186 + t324) * t296 - t319 * t205 + t317 * t261 / 0.2e1 + t242 + t181 * t34 + g(1) * t144 - t17 * t83 - t29 * t84 - t16 * t85 - t28 * t86 - t92 * t65 + t206 - t66 * t230 - m(5) * (t139 * t230 + t23 * t28 + t24 * t29) + (-t12 * t16 - t14 * t17 + t156 * t296 + t181 * t2 - t77 * t92) * m(6) + ((t198 * t3 + (-t12 * t198 + t14 * t201) * qJD(4)) * m(6) - t332 * t257 + t333 * t256 + t334 * t198) * pkin(3); -m(6) * (t77 * t297 + (-t12 + t13) * t14) - t13 * t83 - t23 * t84 + t14 * t85 + t24 * t86 - t65 * t297 + t2 * t315 + pkin(4) * t34 + t206 + (-(-mrSges(6,1) - t315) * t186 + t324) * t296 + (-t126 * t315 + t326) * g(3) + (-t124 * t315 + t325) * g(2); t109 * t83 + t111 * t85 + (g(1) * t197 + t109 * t14 + t111 * t12 + t18 + (-g(2) * t200 + g(3) * t203) * t196) * m(6) + t235;];
tau = t1;
