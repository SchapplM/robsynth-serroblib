% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:41:34
% EndTime: 2019-03-09 01:41:47
% DurationCPUTime: 10.44s
% Computational Cost: add. (4830->515), mult. (10467->629), div. (0->0), fcn. (7211->14), ass. (0->246)
t324 = mrSges(5,2) - mrSges(6,3);
t323 = -mrSges(6,2) + mrSges(5,1);
t153 = sin(pkin(10));
t155 = cos(pkin(10));
t322 = t153 ^ 2 + t155 ^ 2;
t154 = sin(pkin(9));
t133 = pkin(1) * t154 + qJ(3);
t113 = qJD(1) * qJD(3) + qJDD(1) * t133;
t151 = qJ(1) + pkin(9);
t144 = sin(t151);
t267 = g(2) * t144;
t281 = m(6) + m(7);
t321 = -m(5) - t281;
t315 = Ifges(6,5) - Ifges(5,6);
t146 = cos(t151);
t301 = g(1) * t146 + t267;
t159 = sin(qJ(4));
t272 = cos(qJ(4));
t116 = t153 * t272 + t159 * t155;
t108 = t116 * qJD(1);
t251 = pkin(7) * qJD(1);
t124 = t133 * qJD(1);
t99 = t153 * qJD(2) + t155 * t124;
t84 = t155 * t251 + t99;
t244 = t159 * t84;
t142 = t155 * qJD(2);
t83 = t142 + (-t124 - t251) * t153;
t41 = -t272 * t83 + t244;
t178 = pkin(5) * t108 + t41;
t320 = t178 + qJD(5);
t158 = sin(qJ(6));
t161 = cos(qJ(6));
t104 = qJD(6) + t108;
t211 = t272 * t155;
t200 = qJD(1) * t211;
t227 = t153 * t159;
t107 = qJD(1) * t227 - t200;
t85 = -qJD(4) * t158 + t107 * t161;
t52 = -mrSges(7,2) * t104 + mrSges(7,3) * t85;
t86 = qJD(4) * t161 + t107 * t158;
t53 = mrSges(7,1) * t104 - mrSges(7,3) * t86;
t183 = -t158 * t53 + t161 * t52;
t110 = t116 * qJD(4);
t204 = qJDD(1) * t272;
t217 = qJDD(1) * t159;
t72 = qJD(1) * t110 + t153 * t217 - t155 * t204;
t35 = qJD(6) * t85 + qJDD(4) * t161 + t158 * t72;
t224 = qJD(4) * t159;
t208 = t153 * t224;
t71 = qJD(1) * t208 - qJD(4) * t200 - t153 * t204 - t155 * t217;
t70 = qJDD(6) - t71;
t21 = mrSges(7,1) * t70 - mrSges(7,3) * t35;
t36 = -qJD(6) * t86 - qJDD(4) * t158 + t161 * t72;
t22 = -mrSges(7,2) * t70 + mrSges(7,3) * t36;
t319 = qJD(6) * t183 + t158 * t22 + t161 * t21;
t287 = t35 / 0.2e1;
t286 = t36 / 0.2e1;
t285 = t70 / 0.2e1;
t318 = m(3) + m(4);
t150 = pkin(10) + qJ(4);
t145 = cos(t150);
t266 = g(3) * t145;
t317 = t85 * Ifges(7,6);
t316 = Ifges(6,4) - Ifges(5,5);
t47 = -mrSges(7,1) * t85 + mrSges(7,2) * t86;
t258 = mrSges(6,1) * t107;
t95 = -qJD(4) * mrSges(6,3) + t258;
t314 = t47 - t95;
t11 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t57 = mrSges(6,1) * t72 - qJDD(4) * mrSges(6,3);
t313 = -t57 + t11;
t102 = Ifges(5,4) * t107;
t310 = t104 * Ifges(7,3);
t312 = t108 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t86 * Ifges(7,5) - t102 + t310 + t317;
t255 = mrSges(5,3) * t108;
t257 = mrSges(6,1) * t108;
t311 = -qJD(4) * t323 + t255 + t257;
t241 = qJDD(4) / 0.2e1;
t58 = -t71 * mrSges(6,1) + qJDD(4) * mrSges(6,2);
t309 = t319 + t58;
t143 = sin(t150);
t308 = t301 * t143;
t132 = t143 * qJ(5);
t230 = t145 * t146;
t307 = pkin(4) * t230 + t146 * t132;
t305 = t143 * t324 - t323 * t145;
t140 = t155 * qJDD(2);
t87 = -t113 * t153 + t140;
t88 = t153 * qJDD(2) + t155 * t113;
t302 = -t153 * t87 + t155 * t88;
t299 = 0.2e1 * t241;
t191 = mrSges(7,1) * t161 - mrSges(7,2) * t158;
t264 = t107 * pkin(5);
t42 = t159 * t83 + t272 * t84;
t40 = -qJD(4) * qJ(5) - t42;
t25 = -t40 - t264;
t271 = Ifges(7,4) * t86;
t31 = Ifges(7,2) * t85 + Ifges(7,6) * t104 + t271;
t298 = t25 * t191 - t161 * t31 / 0.2e1;
t260 = pkin(7) + t133;
t111 = t260 * t153;
t112 = t260 * t155;
t206 = qJD(4) * t272;
t297 = t159 * (qJD(3) * t153 + qJD(4) * t112) - qJD(3) * t211 + t111 * t206;
t79 = t140 + (-pkin(7) * qJDD(1) - t113) * t153;
t218 = qJDD(1) * t155;
t80 = pkin(7) * t218 + t88;
t215 = t159 * t79 + t83 * t206 + t272 * t80;
t15 = -qJDD(4) * qJ(5) + qJD(4) * (-qJD(5) + t244) - t215;
t194 = -mrSges(4,1) * t155 + mrSges(4,2) * t153;
t296 = -m(4) * pkin(2) - mrSges(3,1) + t194 + t305;
t39 = -qJD(4) * pkin(4) + qJD(5) + t41;
t295 = -m(6) * t39 - t311;
t256 = mrSges(5,3) * t107;
t93 = -qJD(4) * mrSges(5,2) - t256;
t294 = -m(6) * t40 + t93 - t95;
t137 = pkin(3) * t155 + pkin(2);
t156 = cos(pkin(9));
t270 = pkin(1) * t156;
t123 = -t137 - t270;
t103 = qJDD(1) * t123 + qJDD(3);
t166 = qJ(5) * t71 - qJD(5) * t108 + t103;
t280 = pkin(4) + pkin(8);
t12 = t280 * t72 + t166;
t18 = -t159 * t80 - t84 * t206 - t83 * t224 + t272 * t79;
t172 = qJDD(5) - t18;
t7 = -t71 * pkin(5) - qJDD(4) * t280 + t172;
t24 = -qJD(4) * t280 + t320;
t106 = qJD(1) * t123 + qJD(3);
t170 = -qJ(5) * t108 + t106;
t37 = t107 * t280 + t170;
t9 = -t158 * t37 + t161 * t24;
t1 = qJD(6) * t9 + t12 * t161 + t158 * t7;
t10 = t158 * t24 + t161 * t37;
t2 = -qJD(6) * t10 - t12 * t158 + t161 * t7;
t293 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t292 = (-g(1) * t230 - t145 * t267) * qJ(5);
t50 = pkin(4) * t107 + t170;
t291 = -t106 * mrSges(5,1) + t50 * mrSges(6,2);
t290 = -m(4) * qJ(3) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t289 = t9 * mrSges(7,1) + t106 * mrSges(5,2) - t10 * mrSges(7,2) - t50 * mrSges(6,3);
t288 = Ifges(7,1) * t287 + Ifges(7,4) * t286 + Ifges(7,5) * t285;
t284 = -t85 / 0.2e1;
t283 = -t86 / 0.2e1;
t282 = t86 / 0.2e1;
t279 = -t104 / 0.2e1;
t278 = -t107 / 0.2e1;
t277 = t107 / 0.2e1;
t276 = -t108 / 0.2e1;
t275 = t108 / 0.2e1;
t160 = sin(qJ(1));
t269 = pkin(1) * t160;
t265 = t1 * t158;
t134 = t145 * pkin(4);
t162 = cos(qJ(1));
t147 = t162 * pkin(1);
t263 = -qJD(4) / 0.2e1;
t262 = qJD(4) / 0.2e1;
t254 = mrSges(7,3) * t161;
t253 = Ifges(7,4) * t158;
t252 = Ifges(7,4) * t161;
t250 = t108 * Ifges(5,4);
t249 = t108 * Ifges(6,6);
t238 = qJ(5) * t107;
t237 = t108 * t158;
t236 = t110 * t158;
t235 = t110 * t161;
t115 = -t211 + t227;
t234 = t115 * t158;
t233 = t115 * t161;
t232 = t144 * t158;
t231 = t144 * t161;
t229 = t146 * t158;
t228 = t146 * t161;
t226 = t146 * t137 + t147;
t225 = t134 + t132;
t223 = qJD(6) * t158;
t222 = qJD(6) * t161;
t219 = qJDD(1) * t153;
t216 = Ifges(7,5) * t35 + Ifges(7,6) * t36 + Ifges(7,3) * t70;
t214 = t93 + t314;
t213 = m(4) - t321;
t138 = -pkin(2) - t270;
t205 = -t223 / 0.2e1;
t203 = -t137 - t132;
t63 = t272 * t111 + t112 * t159;
t201 = -m(7) * t280 - mrSges(7,3);
t157 = -pkin(7) - qJ(3);
t198 = -t144 * t157 + t226;
t197 = -mrSges(4,1) * t218 + mrSges(4,2) * t219;
t196 = t10 * t161 - t158 * t9;
t195 = t10 * t158 + t161 * t9;
t190 = mrSges(7,1) * t158 + mrSges(7,2) * t161;
t188 = Ifges(7,1) * t158 + t252;
t187 = Ifges(7,2) * t161 + t253;
t186 = Ifges(7,5) * t158 + Ifges(7,6) * t161;
t184 = -t153 * (-t124 * t153 + t142) + t155 * t99;
t171 = -qJ(5) * t116 + t123;
t45 = t115 * t280 + t171;
t48 = pkin(5) * t116 + t63;
t20 = t158 * t48 + t161 * t45;
t19 = -t158 * t45 + t161 * t48;
t182 = -t158 * t52 - t161 * t53;
t109 = -t155 * t206 + t208;
t181 = qJ(5) * t109 - qJD(5) * t116;
t64 = -t159 * t111 + t112 * t272;
t177 = t115 * t222 + t236;
t176 = t115 * t223 - t235;
t175 = -t182 + t311;
t169 = qJD(6) * t196 + t161 * t2 + t265;
t44 = qJD(3) * t116 + qJD(4) * t64;
t121 = qJDD(1) * t138 + qJDD(3);
t101 = Ifges(6,6) * t107;
t92 = -t143 * t232 + t228;
t91 = t143 * t231 + t229;
t90 = t143 * t229 + t231;
t89 = t143 * t228 - t232;
t82 = Ifges(7,4) * t85;
t69 = t71 * mrSges(5,2);
t68 = t71 * mrSges(6,3);
t67 = -mrSges(6,2) * t107 - mrSges(6,3) * t108;
t66 = pkin(4) * t108 + t238;
t61 = -t107 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t250;
t60 = Ifges(6,4) * qJD(4) - t108 * Ifges(6,2) + t101;
t59 = Ifges(6,5) * qJD(4) + t107 * Ifges(6,3) - t249;
t56 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t72;
t55 = qJDD(4) * mrSges(5,1) + mrSges(5,3) * t71;
t54 = pkin(4) * t115 + t171;
t51 = pkin(4) * t110 + t181;
t49 = -t115 * pkin(5) + t64;
t46 = t108 * t280 + t238;
t38 = t110 * t280 + t181;
t32 = Ifges(7,1) * t86 + Ifges(7,5) * t104 + t82;
t29 = -t109 * pkin(5) + t44;
t28 = -pkin(5) * t110 - t297;
t27 = t42 - t264;
t23 = pkin(4) * t72 + t166;
t17 = -t224 * t84 + t215;
t16 = -qJDD(4) * pkin(4) + t172;
t14 = t158 * t27 + t161 * t46;
t13 = -t158 * t46 + t161 * t27;
t8 = -pkin(5) * t72 - t15;
t5 = t35 * Ifges(7,4) + t36 * Ifges(7,2) + t70 * Ifges(7,6);
t4 = -qJD(6) * t20 - t158 * t38 + t161 * t29;
t3 = qJD(6) * t19 + t158 * t29 + t161 * t38;
t6 = [t3 * t52 + t4 * t53 + t28 * t47 + t49 * t11 + t19 * t21 + t20 * t22 + (m(5) * t103 + t72 * mrSges(5,1) - t69) * t123 + (t113 * t322 + t302) * mrSges(4,3) + t121 * t194 + t138 * t197 + (-m(5) * t18 + m(6) * t16 - t55 + t58) * t63 + (qJD(6) * t32 + t5) * t233 / 0.2e1 - t71 * Ifges(5,1) * t116 + m(7) * (t1 * t20 + t10 * t3 + t19 * t2 + t25 * t28 + t4 * t9 + t49 * t8) + (t1 * t233 - t10 * t176 - t177 * t9 - t2 * t234) * mrSges(7,3) + t85 * (Ifges(7,4) * t177 - Ifges(7,2) * t176) / 0.2e1 + (-t312 / 0.2e1 + Ifges(6,2) * t276 + Ifges(6,6) * t277 - Ifges(5,4) * t278 - Ifges(7,5) * t282 - Ifges(5,1) * t275 - t310 / 0.2e1 - t317 / 0.2e1 - t39 * mrSges(6,1) - t41 * mrSges(5,3) + t60 / 0.2e1 + t316 * t262 - t289) * t109 + (mrSges(2,1) * t160 - t92 * mrSges(7,1) + mrSges(2,2) * t162 + t91 * mrSges(7,2) + t318 * t269 + t321 * (-t146 * t157 - t269) + (-m(7) * pkin(5) + t290) * t146 + (m(5) * t137 - m(7) * t203 - t145 * t201 - m(6) * (t203 - t134) - t296) * t144) * g(1) + (t103 * mrSges(5,1) + t15 * mrSges(6,1) - t23 * mrSges(6,2) - t17 * mrSges(5,3) + t186 * t285 + t187 * t286 + t188 * t287 - t8 * t191 + t31 * t205 + (Ifges(6,3) + Ifges(5,2)) * t72 + (Ifges(6,6) + Ifges(5,4)) * t71 + t315 * t299) * t115 + m(4) * (t184 * qJD(3) + t121 * t138 + t133 * t302) + (m(5) * t41 - t295) * t44 + (-m(5) * t42 - t294) * t297 + (0.2e1 * (mrSges(3,1) * t156 - mrSges(3,2) * t154) * pkin(1) + m(3) * (t154 ^ 2 + t156 ^ 2) * pkin(1) ^ 2 + Ifges(3,3) + Ifges(2,3)) * qJDD(1) + t31 * t235 / 0.2e1 + t32 * t236 / 0.2e1 + (-Ifges(5,4) * t72 + Ifges(5,5) * qJDD(4) + t216) * t116 / 0.2e1 + (t16 * mrSges(6,1) + t103 * mrSges(5,2) - t18 * mrSges(5,3) - t23 * mrSges(6,3) + Ifges(5,5) * t241 + Ifges(7,5) * t287 + Ifges(7,6) * t286 + Ifges(7,3) * t285 + (-Ifges(5,4) / 0.2e1 - Ifges(6,6)) * t72 - Ifges(6,2) * t71 - t299 * Ifges(6,4) + t293) * t116 + (m(5) * t17 - m(6) * t15 + t56 - t57) * t64 + t104 * (Ifges(7,5) * t177 - Ifges(7,6) * t176) / 0.2e1 + (Ifges(4,4) * t153 + Ifges(4,2) * t155) * t218 + (Ifges(4,1) * t153 + Ifges(4,4) * t155) * t219 + t25 * (mrSges(7,1) * t176 + mrSges(7,2) * t177) + t234 * t288 + t51 * t67 + t54 * (-t72 * mrSges(6,2) + t68) + (Ifges(6,6) * t276 + Ifges(6,3) * t277 - Ifges(5,2) * t278 - Ifges(5,4) * t275 + t40 * mrSges(6,1) - t42 * mrSges(5,3) + t59 / 0.2e1 - t61 / 0.2e1 + t315 * t262 - t291) * t110 + (-m(6) * (t198 + t307) - m(5) * t198 - m(7) * (pkin(8) * t230 + t226 + t307) - t90 * mrSges(7,1) - t89 * mrSges(7,2) - mrSges(7,3) * t230 - mrSges(2,1) * t162 + mrSges(2,2) * t160 - t318 * t147 + t296 * t146 + (-m(7) * (pkin(5) - t157) + t290) * t144) * g(2) + m(6) * (t23 * t54 + t50 * t51) + (Ifges(7,1) * t177 - Ifges(7,4) * t176) * t282; m(3) * qJDD(2) + (t56 + t313) * t116 - t214 * t109 + t175 * t110 + (-m(3) - t213) * g(3) + (-t55 + t309) * t115 + m(5) * (-t109 * t42 + t110 * t41 - t115 * t18 + t116 * t17) + m(6) * (t109 * t40 + t110 * t39 + t115 * t16 - t116 * t15) + m(4) * (t153 * t88 + t155 * t87) + m(7) * (-t109 * t25 + t110 * t195 + t115 * t169 + t116 * t8); -t158 * t21 + t161 * t22 + t68 - t69 + t323 * t72 + t182 * qJD(6) - t322 * qJD(1) ^ 2 * mrSges(4,3) + t214 * t107 - t175 * t108 + t197 + (-g(1) * t144 + g(2) * t146) * t213 + (t1 * t161 - t104 * t195 + t107 * t25 - t2 * t158) * m(7) + (-t107 * t40 - t108 * t39 + t23) * m(6) + (t107 * t42 - t108 * t41 + t103) * m(5) + (-qJD(1) * t184 + t121) * m(4); -t14 * t52 - t13 * t53 - pkin(4) * t58 + t16 * mrSges(6,2) - t17 * mrSges(5,2) + t18 * mrSges(5,1) - t15 * mrSges(6,3) + (Ifges(7,5) * t161 - Ifges(7,6) * t158) * t285 + (-Ifges(7,2) * t158 + t252) * t286 + (t249 + t61) * t275 + t8 * t190 - (t104 * t186 + t187 * t85 + t188 * t86) * qJD(6) / 0.2e1 + (-m(7) * t169 - t319) * t280 + (qJ(5) * t8 - t10 * t14 - t13 * t9 + t320 * t25 + t292) * m(7) + (-t250 + t59) * t276 + (t101 + t60) * t278 + ((-t190 + t324) * t145 + (m(6) * pkin(4) - t201 + t323) * t143) * t301 + (Ifges(5,3) + Ifges(6,1)) * qJDD(4) + t298 * qJD(6) + (-m(6) * t225 - m(7) * (pkin(8) * t145 + t225) - t190 * t143 + t305) * g(3) + (t255 + t295) * t42 + t178 * t47 + t39 * t258 + (-t237 / 0.2e1 + t205) * t32 + (-pkin(4) * t16 - qJ(5) * t15 - qJD(5) * t40 - t50 * t66 + t292) * m(6) + (t256 + t294) * t41 + (Ifges(7,1) * t161 - t253) * t287 + t161 * t288 - t66 * t67 - t2 * t254 + (-Ifges(5,2) * t108 - t102 + t312) * t277 + t313 * qJ(5) + t314 * qJD(5) + t315 * t72 + (Ifges(6,3) * t278 - t10 * t254 + t186 * t279 + t187 * t284 + t188 * t283 + t263 * t315 + t291 + t298) * t108 + (-Ifges(5,1) * t276 - Ifges(7,5) * t283 + Ifges(6,2) * t275 - Ifges(7,6) * t284 - Ifges(7,3) * t279 + t263 * t316 + t289) * t107 + t316 * t71 - t40 * t257 + (-t266 - t10 * t222 - t265 + (t223 + t237) * t9) * mrSges(7,3) - t158 * t5 / 0.2e1; -t314 * qJD(4) + t281 * t266 + (t183 + t67) * t108 + (-qJD(4) * t25 + t108 * t196 + t169 - t308) * m(7) + (qJD(4) * t40 + t108 * t50 + t16 - t308) * m(6) + t309; -t25 * (mrSges(7,1) * t86 + mrSges(7,2) * t85) + (Ifges(7,1) * t85 - t271) * t283 + t31 * t282 + (Ifges(7,5) * t85 - Ifges(7,6) * t86) * t279 - t9 * t52 + t10 * t53 - g(1) * (mrSges(7,1) * t89 - mrSges(7,2) * t90) - g(2) * (mrSges(7,1) * t91 + mrSges(7,2) * t92) + t191 * t266 + (t10 * t86 + t85 * t9) * mrSges(7,3) + t216 + (-Ifges(7,2) * t86 + t32 + t82) * t284 + t293;];
tau  = t6;
