% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRPPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:09:59
% EndTime: 2018-11-23 15:10:06
% DurationCPUTime: 6.94s
% Computational Cost: add. (3946->547), mult. (10368->761), div. (0->0), fcn. (7071->10), ass. (0->254)
t331 = Ifges(5,1) + Ifges(6,1);
t185 = cos(pkin(11));
t192 = cos(qJ(2));
t184 = sin(pkin(6));
t257 = qJD(1) * t184;
t236 = t192 * t257;
t189 = sin(qJ(2));
t237 = t189 * t257;
t183 = sin(pkin(11));
t191 = cos(qJ(3));
t270 = t183 * t191;
t112 = -t185 * t237 + t236 * t270;
t188 = sin(qJ(3));
t238 = -pkin(8) * t183 - pkin(4);
t267 = t185 * t191;
t243 = pkin(9) * t267;
t215 = pkin(3) * t188 - qJ(4) * t191;
t135 = qJD(3) * t215 - qJD(4) * t188;
t271 = t135 * t185;
t333 = -t112 - t271 + (-t243 + (-pkin(5) + t238) * t188) * qJD(3);
t266 = t185 * t192;
t113 = (t183 * t189 + t191 * t266) * t257;
t127 = t183 * t135;
t244 = pkin(9) * t270;
t268 = t185 * t188;
t248 = qJD(5) * t191;
t250 = qJD(3) * t188;
t313 = qJ(5) * t250 - t248;
t332 = t113 - t127 - (-pkin(8) * t268 + t244) * qJD(3) - t313;
t325 = Ifges(6,4) + Ifges(5,5);
t187 = sin(qJ(6));
t190 = cos(qJ(6));
t160 = -pkin(3) * t191 - qJ(4) * t188 - pkin(2);
t173 = pkin(8) * t270;
t182 = t191 * pkin(4);
t292 = pkin(9) * t188;
t71 = pkin(5) * t191 + t173 + t182 + (-t160 - t292) * t185;
t118 = pkin(8) * t267 + t183 * t160;
t109 = -qJ(5) * t191 + t118;
t85 = t183 * t292 + t109;
t26 = -t187 * t85 + t190 * t71;
t330 = qJD(6) * t26 + t187 * t333 - t190 * t332;
t27 = t187 * t71 + t190 * t85;
t329 = -qJD(6) * t27 + t187 * t332 + t190 * t333;
t328 = -qJD(3) / 0.2e1;
t327 = qJD(3) / 0.2e1;
t326 = -Ifges(5,4) + Ifges(6,5);
t287 = -pkin(9) + qJ(4);
t162 = t287 * t183;
t163 = t287 * t185;
t103 = t162 * t187 + t163 * t190;
t151 = t183 * t190 - t185 * t187;
t300 = pkin(4) + pkin(5);
t198 = (-t188 * t300 - t243) * qJD(2);
t155 = qJD(2) * pkin(8) + t237;
t146 = t188 * t155;
t186 = cos(pkin(6));
t265 = t186 * t191;
t171 = qJD(1) * t265;
t114 = -t146 + t171;
t153 = t215 * qJD(2);
t61 = -t183 * t114 + t153 * t185;
t37 = t198 - t61;
t253 = qJD(2) * t188;
t62 = t185 * t114 + t183 * t153;
t49 = qJ(5) * t253 + t62;
t43 = qJD(2) * t244 + t49;
t324 = qJD(4) * t151 - qJD(6) * t103 + t187 * t43 - t190 * t37;
t102 = t162 * t190 - t163 * t187;
t209 = t183 * t187 + t185 * t190;
t323 = qJD(4) * t209 + qJD(6) * t102 - t187 * t37 - t190 * t43;
t281 = Ifges(6,5) * t183;
t285 = Ifges(5,4) * t183;
t322 = t185 * t331 + t281 - t285;
t204 = t209 * t191;
t201 = qJD(3) * t204;
t148 = -t185 * qJD(3) + t183 * t253;
t247 = t183 * qJD(3);
t149 = t185 * t253 + t247;
t81 = t148 * t190 - t149 * t187;
t40 = qJD(2) * t201 + qJD(6) * t81;
t306 = t40 / 0.2e1;
t205 = t151 * t191;
t202 = qJD(3) * t205;
t210 = t148 * t187 + t149 * t190;
t41 = qJD(2) * t202 - qJD(6) * t210;
t305 = t41 / 0.2e1;
t320 = t149 / 0.2e1;
t246 = qJD(2) * qJD(3);
t229 = t188 * t246;
t319 = -t229 / 0.2e1;
t318 = t246 / 0.2e1;
t316 = t253 / 0.2e1;
t231 = Ifges(4,6) * t328;
t315 = Ifges(4,5) * t327;
t314 = -qJD(2) / 0.2e1;
t251 = qJD(2) * t191;
t175 = qJD(6) + t251;
t254 = qJD(2) * t184;
t230 = qJD(1) * t254;
t224 = t192 * t230;
t258 = qJD(3) * t171 + t191 * t224;
t65 = (qJD(4) - t146) * qJD(3) + t258;
t89 = (t135 + t237) * qJD(2);
t23 = -t183 * t65 + t185 * t89;
t14 = qJD(3) * t198 - t23;
t24 = t183 * t89 + t185 * t65;
t242 = qJ(5) * t229 + t24;
t15 = (pkin(9) * t247 - qJD(5)) * t251 + t242;
t256 = qJD(1) * t188;
t235 = t186 * t256;
t115 = t155 * t191 + t235;
t110 = qJD(3) * qJ(4) + t115;
t116 = qJD(2) * t160 - t236;
t44 = -t183 * t110 + t116 * t185;
t35 = pkin(4) * t251 + qJD(5) - t44;
t18 = pkin(5) * t251 - pkin(9) * t149 + t35;
t45 = t185 * t110 + t183 * t116;
t36 = -qJ(5) * t251 + t45;
t25 = pkin(9) * t148 + t36;
t5 = t18 * t190 - t187 * t25;
t1 = qJD(6) * t5 + t14 * t187 + t15 * t190;
t6 = t18 * t187 + t190 * t25;
t2 = -qJD(6) * t6 + t14 * t190 - t15 * t187;
t312 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t40 + Ifges(7,6) * t41;
t228 = t191 * t246;
t222 = t185 * t228;
t223 = t183 * t228;
t125 = mrSges(6,1) * t223 - mrSges(6,3) * t222;
t126 = mrSges(5,1) * t223 + mrSges(5,2) * t222;
t13 = -t41 * mrSges(7,1) + t40 * mrSges(7,2);
t311 = t125 + t126 - t13;
t156 = -qJD(2) * pkin(2) - t236;
t179 = Ifges(4,4) * t251;
t207 = qJD(3) * pkin(3) - qJD(4) + t114;
t280 = Ifges(6,5) * t185;
t216 = Ifges(6,3) * t183 + t280;
t284 = Ifges(5,4) * t185;
t217 = -Ifges(5,2) * t183 + t284;
t294 = t185 / 0.2e1;
t295 = t183 / 0.2e1;
t296 = -t183 / 0.2e1;
t199 = qJ(5) * t149 + t207;
t42 = pkin(4) * t148 - t199;
t310 = -(t183 * t36 - t185 * t35) * mrSges(6,2) - (t183 * t45 + t185 * t44) * mrSges(5,3) - t207 * (mrSges(5,1) * t183 + mrSges(5,2) * t185) + t156 * mrSges(4,2) + t42 * (mrSges(6,1) * t183 - mrSges(6,3) * t185) + Ifges(4,1) * t316 + t179 / 0.2e1 + t315 - t114 * mrSges(4,3) + (t149 * Ifges(6,5) - Ifges(6,6) * t251) * t295 + (t149 * Ifges(5,4) - Ifges(5,6) * t251) * t296 + t322 * t320 + (t149 * t331 - t325 * t251) * t294 + (Ifges(6,3) * t295 - Ifges(5,2) * t296 + t326 * t294 - t217 / 0.2e1 + t216 / 0.2e1) * t148;
t309 = 0.2e1 * (-Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t148 - (Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t149 - t156 * mrSges(4,1) - t36 * mrSges(6,3) - t44 * mrSges(5,1) - t6 * mrSges(7,2) - t231 - (Ifges(4,4) * t188 + Ifges(4,2) * t191) * t314 + t175 * Ifges(7,3) + t210 * Ifges(7,5) + t81 * Ifges(7,6) + t115 * mrSges(4,3) + t35 * mrSges(6,1) + t45 * mrSges(5,2) + t5 * mrSges(7,1) - t325 * t320 + (Ifges(5,3) + Ifges(6,2)) * t251 / 0.2e1;
t308 = Ifges(7,4) * t306 + Ifges(7,2) * t305 + Ifges(7,6) * t319;
t307 = Ifges(7,1) * t306 + Ifges(7,4) * t305 + Ifges(7,5) * t319;
t304 = -t81 / 0.2e1;
t303 = t81 / 0.2e1;
t302 = -t210 / 0.2e1;
t301 = t210 / 0.2e1;
t298 = -t175 / 0.2e1;
t297 = t175 / 0.2e1;
t293 = Ifges(7,4) * t210;
t30 = -mrSges(7,1) * t81 + mrSges(7,2) * t210;
t87 = mrSges(6,1) * t148 - mrSges(6,3) * t149;
t286 = t30 - t87;
t283 = Ifges(6,4) * t185;
t282 = Ifges(5,5) * t185;
t279 = Ifges(5,6) * t183;
t278 = Ifges(6,6) * t183;
t269 = t184 * t189;
t138 = t188 * t269 - t265;
t233 = t192 * t254;
t249 = qJD(3) * t191;
t70 = t155 * t249 + (qJD(3) * t186 + t233) * t256;
t277 = t138 * t70;
t275 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t148 + mrSges(5,2) * t149 + mrSges(4,3) * t253;
t272 = qJ(5) * t185;
t119 = -mrSges(6,2) * t148 - mrSges(6,3) * t251;
t120 = mrSges(5,2) * t251 - mrSges(5,3) * t148;
t264 = t119 + t120;
t121 = -mrSges(5,1) * t251 - mrSges(5,3) * t149;
t122 = mrSges(6,1) * t251 + mrSges(6,2) * t149;
t263 = t121 - t122;
t123 = qJD(2) * t205;
t137 = t151 * qJD(6);
t262 = t123 + t137;
t124 = qJD(2) * t204;
t136 = t209 * qJD(6);
t261 = t124 + t136;
t131 = (-mrSges(5,2) * t188 - mrSges(5,3) * t270) * t246;
t134 = (-mrSges(6,2) * t270 + mrSges(6,3) * t188) * t246;
t260 = t131 + t134;
t132 = (mrSges(5,1) * t188 - mrSges(5,3) * t267) * t246;
t166 = mrSges(6,2) * t222;
t133 = -mrSges(6,1) * t229 + t166;
t259 = -t132 + t133;
t255 = qJD(2) * t183;
t252 = qJD(2) * t189;
t245 = pkin(8) * t188 * t70;
t241 = pkin(8) * t250;
t234 = t184 * t252;
t232 = qJD(5) * t268;
t227 = qJ(5) * t183 + pkin(3);
t226 = t275 - t286;
t117 = t160 * t185 - t173;
t225 = t188 * t236;
t214 = pkin(4) * t183 - t272;
t55 = -mrSges(7,2) * t175 + mrSges(7,3) * t81;
t56 = mrSges(7,1) * t175 - mrSges(7,3) * t210;
t211 = -t187 * t56 + t190 * t55;
t139 = t186 * t188 + t191 * t269;
t94 = t139 * t183 + t184 * t266;
t95 = -t183 * t184 * t192 + t139 * t185;
t33 = -t187 * t95 + t190 * t94;
t34 = t187 * t94 + t190 * t95;
t93 = -t185 * t241 + t127;
t208 = pkin(8) + t214;
t206 = -t183 * t300 + t272;
t129 = t151 * t188;
t203 = -pkin(8) + t206;
t200 = -qJ(5) * t222 - t149 * qJD(5) + t188 * t224;
t16 = -qJD(2) * t248 + t242;
t28 = (t235 + (pkin(4) * t255 + t155) * t191) * qJD(3) + t200;
t197 = (Ifges(6,6) * t188 + t191 * t216) * t318 - (Ifges(5,6) * t188 + t191 * t217) * t246 / 0.2e1 + t28 * mrSges(6,1) + t70 * mrSges(5,1) - t16 * mrSges(6,2) - t24 * mrSges(5,3);
t17 = -pkin(4) * t229 - t23;
t196 = t70 * mrSges(5,2) + t17 * mrSges(6,2) - t23 * mrSges(5,3) - t28 * mrSges(6,3) + (t188 * t325 + t191 * t322) * t318;
t193 = qJD(2) ^ 2;
t169 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t251;
t157 = -pkin(4) * t185 - t227;
t152 = (-mrSges(4,1) * t191 + mrSges(4,2) * t188) * qJD(2);
t145 = (mrSges(4,1) * t188 + mrSges(4,2) * t191) * t246;
t141 = t185 * t300 + t227;
t130 = t209 * t188;
t128 = t208 * t188;
t111 = -t117 + t182;
t105 = t203 * t188;
t97 = qJD(3) * t139 + t188 * t233;
t96 = -qJD(3) * t138 + t191 * t233;
t92 = t183 * t241 + t271;
t86 = t208 * t249 - t232;
t79 = t238 * t250 - t271;
t78 = Ifges(7,4) * t81;
t69 = -t155 * t250 + t258;
t68 = t203 * t249 + t232;
t67 = t235 + (qJD(2) * t214 + t155) * t191;
t66 = t93 + t313;
t64 = -t136 * t188 + t202;
t63 = qJD(6) * t129 + t201;
t60 = t183 * t234 + t185 * t96;
t59 = t183 * t96 - t185 * t234;
t52 = -pkin(4) * t253 - t61;
t51 = -t235 + (qJD(2) * t206 - t155) * t191;
t32 = mrSges(7,2) * t229 + mrSges(7,3) * t41;
t31 = -mrSges(7,1) * t229 - mrSges(7,3) * t40;
t29 = -t148 * t300 + t199;
t22 = Ifges(7,1) * t210 + Ifges(7,5) * t175 + t78;
t21 = Ifges(7,2) * t81 + Ifges(7,6) * t175 + t293;
t19 = (t235 + (t255 * t300 + t155) * t191) * qJD(3) + t200;
t8 = -qJD(6) * t34 - t187 * t60 + t190 * t59;
t7 = qJD(6) * t33 + t187 * t59 + t190 * t60;
t3 = [-t139 * mrSges(4,3) * t229 + t96 * t169 + t33 * t31 + t34 * t32 + t7 * t55 + t8 * t56 + t260 * t95 + t259 * t94 + t264 * t60 - t263 * t59 + t226 * t97 + ((-mrSges(3,2) * t193 - t145) * t192 + (-mrSges(3,1) * t193 + qJD(2) * t152) * t189) * t184 + (mrSges(4,3) * t228 + t311) * t138 + m(6) * (t138 * t28 + t16 * t95 + t17 * t94 + t35 * t59 + t36 * t60 + t42 * t97) + m(5) * (-t207 * t97 - t23 * t94 + t24 * t95 - t44 * t59 + t45 * t60 + t277) + m(7) * (t1 * t34 + t138 * t19 + t2 * t33 - t29 * t97 + t5 * t8 + t6 * t7) + m(4) * (-t114 * t97 + t115 * t96 + t277 + t139 * t69 + (t156 - t236) * t234); t329 * t56 + (t1 * t27 - t105 * t19 + t2 * t26 + t330 * t6 + t329 * t5 + (t225 + t68) * t29) * m(7) + t330 * t55 - t152 * t237 + t263 * t112 - t264 * t113 + (t70 * mrSges(4,3) + pkin(8) * t126 + t196 * t185 + t197 * t183 + (mrSges(4,2) * t252 - t192 * t226) * t257) * t188 + (t109 * t16 + t111 * t17 + t128 * t28 + (-t225 + t86) * t42 + (-t113 + t66) * t36 + (-t112 + t79) * t35) * m(6) + m(4) * (-pkin(2) * t189 * t230 + pkin(8) * t191 * t69 + t245) + m(5) * (t117 * t23 + t118 * t24 + t44 * t92 + t45 * t93 + t245) - m(4) * (t156 * t189 + (-t114 * t188 + t115 * t191) * t192) * t257 + (((-0.3e1 / 0.2e1 * qJD(2) * Ifges(4,4) + (t278 + t283 - t279 + t282) * qJD(2) / 0.2e1) * t188 + (Ifges(7,5) * t130 + Ifges(7,6) * t129) * t314 + (-m(4) * t115 - t169) * pkin(8) + t231 - t309) * t188 + (((0.3e1 / 0.2e1 * Ifges(4,4) + (-0.3e1 / 0.2e1 * Ifges(5,5) - 0.3e1 / 0.2e1 * Ifges(6,4)) * t185 + (0.3e1 / 0.2e1 * Ifges(5,6) - 0.3e1 / 0.2e1 * Ifges(6,6)) * t183) * t191 + (0.3e1 / 0.2e1 * Ifges(4,1) - Ifges(7,3) - 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(6,2) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t185 ^ 2 + ((Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t183 + t326 * t185) * t183) * t188) * qJD(2) + (-m(4) * t114 - m(5) * t207 + t275) * pkin(8) + t315 + t310) * t191) * qJD(3) - m(5) * (-t112 * t44 + t113 * t45 - t207 * t225) + (-t16 * mrSges(6,3) + t24 * mrSges(5,2) + t17 * mrSges(6,1) - t23 * mrSges(5,1) + t69 * mrSges(4,3) + (-mrSges(4,1) * t252 - t169 * t192) * t257 + t312) * t191 + t26 * t31 + t27 * t32 + t63 * t22 / 0.2e1 + t64 * t21 / 0.2e1 + t29 * (-mrSges(7,1) * t64 + mrSges(7,2) * t63) + t68 * t30 + t86 * t87 + (Ifges(7,5) * t63 + Ifges(7,6) * t64) * t297 + (Ifges(7,1) * t63 + Ifges(7,4) * t64) * t301 + (Ifges(7,4) * t63 + Ifges(7,2) * t64) * t303 + (Ifges(7,4) * t130 + Ifges(7,2) * t129) * t305 + (Ifges(7,1) * t130 + Ifges(7,4) * t129) * t306 + t130 * t307 + t129 * t308 + t105 * t13 + t66 * t119 + t93 * t120 + t92 * t121 + t79 * t122 + t128 * t125 - t19 * (-mrSges(7,1) * t129 + mrSges(7,2) * t130) + t118 * t131 + t117 * t132 + t111 * t133 + t109 * t134 - pkin(2) * t145 + (t1 * t129 - t130 * t2 - t5 * t63 + t6 * t64) * mrSges(7,3); (-t124 / 0.2e1 - t136 / 0.2e1) * t22 + (-Ifges(7,5) * t136 - Ifges(7,6) * t137) * t297 + (-Ifges(7,1) * t136 - Ifges(7,4) * t137) * t301 + (-Ifges(7,4) * t136 - Ifges(7,2) * t137) * t303 + (-t123 / 0.2e1 - t137 / 0.2e1) * t21 - t275 * t115 + (mrSges(7,1) * t262 - mrSges(7,2) * t261) * t29 + (qJ(4) * t260 + qJD(4) * t264 - t197) * t185 + (t157 * t28 + (qJ(4) * t16 + qJD(4) * t36) * t185 + (qJ(4) * t17 + qJD(4) * t35 - qJD(5) * t42) * t183 - t35 * t52 - t36 * t49 - t42 * t67) * m(6) + (((Ifges(7,5) * t151 - Ifges(7,6) * t209) * t328 + Ifges(4,4) * t316 + t231 + ((Ifges(5,6) - Ifges(6,6)) * t185 + t325 * t183) * t327 + t309) * t188 + ((-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t253 - t179 / 0.2e1 + ((Ifges(5,2) * t185 + t285) * t296 + (-Ifges(6,3) * t185 + t281) * t295 + Ifges(4,5) / 0.2e1 + (t183 * t331 - t280 + t284) * t294) * qJD(3) + (t282 / 0.2e1 - t279 / 0.2e1 + t283 / 0.2e1 + t278 / 0.2e1) * t251 - t310) * t191) * qJD(2) + t324 * t56 + t323 * t55 + (t1 * t103 + t102 * t2 - t141 * t19 + t323 * t6 + t324 * t5 + (qJD(5) * t183 - t51) * t29) * m(7) + (t207 * t115 - t44 * t61 - t45 * t62 - pkin(3) * t70 + (-t183 * t44 + t185 * t45) * qJD(4) + (-t183 * t23 + t185 * t24) * qJ(4)) * m(5) + (-t1 * t209 - t151 * t2 + t261 * t5 - t262 * t6) * mrSges(7,3) + (Ifges(7,4) * t151 - Ifges(7,2) * t209) * t305 + (Ifges(7,1) * t151 - Ifges(7,4) * t209) * t306 - t19 * (mrSges(7,1) * t209 + mrSges(7,2) * t151) - t209 * t308 + (qJ(4) * t259 - qJD(4) * t263 + qJD(5) * t286 + t196) * t183 - t51 * t30 - t69 * mrSges(4,2) - t70 * mrSges(4,1) - t67 * t87 + (Ifges(7,5) * t124 + Ifges(7,6) * t123) * t298 + (Ifges(7,1) * t124 + Ifges(7,4) * t123) * t302 + (Ifges(7,4) * t124 + Ifges(7,2) * t123) * t304 + t151 * t307 + t102 * t31 + t103 * t32 - t49 * t119 - t62 * t120 - t61 * t121 - t52 * t122 - pkin(3) * t126 + t141 * t13 + t157 * t125 - t114 * t169; t264 * t148 + t263 * t149 + t81 * t55 - t210 * t56 + (-t210 * t5 + t6 * t81 + t19) * m(7) + (t148 * t36 - t149 * t35 + t28) * m(6) + (t148 * t45 + t149 * t44 + t70) * m(5) + t311; t187 * t32 + t190 * t31 + t166 - t286 * t149 + t211 * qJD(6) + (-mrSges(6,1) * t250 + (t119 + t211) * t191) * qJD(2) + (t1 * t187 - t149 * t29 + t190 * t2 - t175 * (t187 * t5 - t190 * t6)) * m(7) + (t149 * t42 + t251 * t36 + t17) * m(6); -Ifges(7,3) * t229 - t29 * (mrSges(7,1) * t210 + mrSges(7,2) * t81) + (Ifges(7,1) * t81 - t293) * t302 + t21 * t301 + (Ifges(7,5) * t81 - Ifges(7,6) * t210) * t298 - t5 * t55 + t6 * t56 + (t210 * t6 + t5 * t81) * mrSges(7,3) + (-Ifges(7,2) * t210 + t22 + t78) * t304 + t312;];
tauc  = t3(:);
