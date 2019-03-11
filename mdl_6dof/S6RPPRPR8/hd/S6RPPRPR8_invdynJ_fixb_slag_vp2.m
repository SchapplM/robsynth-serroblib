% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:20
% EndTime: 2019-03-09 01:55:35
% DurationCPUTime: 10.24s
% Computational Cost: add. (4773->517), mult. (9748->627), div. (0->0), fcn. (6436->10), ass. (0->244)
t331 = -mrSges(6,2) + mrSges(5,1);
t150 = sin(pkin(9));
t151 = cos(pkin(9));
t225 = t150 ^ 2 + t151 ^ 2;
t203 = t225 * mrSges(4,3);
t153 = -pkin(1) - qJ(3);
t330 = -qJD(1) * qJD(3) + qJDD(1) * t153;
t154 = sin(qJ(6));
t156 = cos(qJ(6));
t157 = cos(qJ(4));
t280 = sin(qJ(4));
t106 = t157 * t150 + t151 * t280;
t211 = t280 * t150;
t123 = qJD(1) * t211;
t224 = qJD(4) * t157;
t208 = t151 * t224;
t71 = qJD(1) * t208 - qJD(4) * t123 + qJDD(1) * t106;
t100 = t106 * qJD(1);
t78 = -qJD(4) * t154 + t100 * t156;
t30 = qJD(6) * t78 + qJDD(4) * t156 + t154 * t71;
t102 = t106 * qJD(4);
t218 = qJDD(1) * t151;
t70 = qJD(1) * t102 + qJDD(1) * t211 - t157 * t218;
t69 = qJDD(6) - t70;
t15 = mrSges(7,1) * t69 - mrSges(7,3) * t30;
t79 = qJD(4) * t156 + t100 * t154;
t31 = -qJD(6) * t79 - qJDD(4) * t154 + t156 * t71;
t16 = -mrSges(7,2) * t69 + mrSges(7,3) * t31;
t234 = t151 * t157;
t101 = qJD(1) * t234 - t123;
t94 = qJD(6) + t101;
t48 = -mrSges(7,2) * t94 + mrSges(7,3) * t78;
t49 = mrSges(7,1) * t94 - mrSges(7,3) * t79;
t316 = -t154 * t49 + t156 * t48;
t329 = qJD(6) * t316 + t156 * t15 + t154 * t16;
t324 = Ifges(6,5) - Ifges(5,6);
t145 = pkin(9) + qJ(4);
t137 = sin(t145);
t138 = cos(t145);
t328 = mrSges(5,2) * t138 + t331 * t137;
t107 = -t211 + t234;
t295 = t30 / 0.2e1;
t294 = t31 / 0.2e1;
t293 = t69 / 0.2e1;
t288 = -m(6) - m(7);
t327 = t78 * Ifges(7,6);
t326 = t94 * Ifges(7,3);
t325 = Ifges(6,4) - Ifges(5,5);
t40 = -mrSges(7,1) * t78 + mrSges(7,2) * t79;
t265 = mrSges(6,1) * t100;
t88 = -qJD(4) * mrSges(6,3) + t265;
t323 = t40 - t88;
t56 = mrSges(6,1) * t71 - qJDD(4) * mrSges(6,3);
t322 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t71 - t56;
t93 = Ifges(5,4) * t100;
t321 = t101 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t79 * Ifges(7,5) + t326 + t327 - t93;
t259 = t101 * mrSges(6,1);
t264 = mrSges(5,3) * t101;
t320 = -t331 * qJD(4) + t259 + t264;
t242 = qJDD(4) / 0.2e1;
t158 = cos(qJ(1));
t276 = g(2) * t158;
t155 = sin(qJ(1));
t277 = g(1) * t155;
t315 = t276 - t277;
t319 = t138 * t315;
t119 = qJD(1) * t153 + qJD(2);
t202 = -pkin(7) * qJD(1) + t119;
t90 = t202 * t150;
t213 = t280 * t90;
t91 = t202 * t151;
t52 = -t157 * t91 + t213;
t37 = -t101 * pkin(5) - t52;
t318 = qJD(5) - t37;
t147 = qJD(1) * qJD(2);
t121 = -qJDD(1) * qJ(2) - t147;
t312 = 0.2e1 * t242;
t189 = mrSges(7,1) * t156 - mrSges(7,2) * t154;
t271 = t79 * Ifges(7,4);
t25 = t78 * Ifges(7,2) + t94 * Ifges(7,6) + t271;
t273 = t100 * pkin(5);
t53 = t157 * t90 + t280 * t91;
t46 = -qJD(4) * qJ(5) - t53;
t35 = -t46 - t273;
t310 = t35 * t189 - t156 * t25 / 0.2e1;
t267 = -pkin(7) + t153;
t110 = t267 * t150;
t111 = t267 * t151;
t206 = qJD(4) * t280;
t309 = qJD(3) * t106 + t110 * t206 - t111 * t224;
t109 = qJDD(2) + t330;
t198 = -pkin(7) * qJDD(1) + t109;
t83 = t198 * t150;
t84 = t198 * t151;
t216 = t157 * t83 + t91 * t224 + t280 * t84;
t17 = -qJDD(4) * qJ(5) + qJD(4) * (t213 - qJD(5)) - t216;
t287 = pkin(4) + pkin(8);
t32 = -qJD(4) * t287 + t318;
t135 = qJD(1) * qJ(2) + qJD(3);
t139 = t150 * pkin(3);
t112 = qJD(1) * t139 + t135;
t173 = -qJ(5) * t101 + t112;
t36 = t100 * t287 + t173;
t11 = -t154 * t36 + t156 * t32;
t12 = t154 * t32 + t156 * t36;
t182 = t11 * t154 - t12 * t156;
t114 = qJDD(3) - t121;
t219 = qJDD(1) * t150;
t104 = pkin(3) * t219 + t114;
t162 = qJ(5) * t70 - qJD(5) * t101 + t104;
t10 = t287 * t71 + t162;
t23 = t157 * t84 - t91 * t206 - t90 * t224 - t280 * t83;
t168 = qJDD(5) - t23;
t8 = -pkin(5) * t70 - qJDD(4) * t287 + t168;
t2 = -qJD(6) * t12 - t10 * t154 + t156 * t8;
t1 = qJD(6) * t11 + t10 * t156 + t154 * t8;
t274 = t1 * t154;
t308 = qJD(6) * t182 - t156 * t2 - t274;
t43 = -qJD(4) * pkin(4) + qJD(5) + t52;
t307 = -m(6) * t43 - t320;
t260 = t100 * mrSges(5,3);
t86 = -qJD(4) * mrSges(5,2) - t260;
t306 = -m(6) * t46 + t86 - t88;
t103 = -t150 * t206 + t208;
t22 = -t206 * t90 + t216;
t305 = -t102 * t52 - t103 * t53 - t106 * t22 - t107 * t23;
t19 = -qJDD(4) * pkin(4) + t168;
t304 = -t102 * t43 + t103 * t46 + t106 * t17 + t107 * t19;
t303 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t57 = -t70 * mrSges(6,1) + qJDD(4) * mrSges(6,2);
t302 = t329 + t57;
t47 = pkin(4) * t100 + t173;
t301 = t112 * mrSges(5,1) - t47 * mrSges(6,2);
t127 = t138 * qJ(5);
t192 = mrSges(4,1) * t150 + mrSges(4,2) * t151;
t300 = -t192 - mrSges(3,3) - m(7) * (pkin(8) * t137 - t127) - t137 * mrSges(7,3) - (-m(6) * qJ(5) - mrSges(6,3)) * t138 + mrSges(2,2) - t328;
t299 = -m(7) * pkin(5) - mrSges(2,1) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t298 = m(4) * t135 + m(5) * t112 + mrSges(5,1) * t100 + mrSges(5,2) * t101 + t192 * qJD(1);
t297 = t11 * mrSges(7,1) + t112 * mrSges(5,2) - t12 * mrSges(7,2) - t47 * mrSges(6,3);
t296 = Ifges(7,1) * t295 + Ifges(7,4) * t294 + Ifges(7,5) * t293;
t292 = -t78 / 0.2e1;
t291 = -t79 / 0.2e1;
t290 = t79 / 0.2e1;
t289 = -t94 / 0.2e1;
t286 = -t100 / 0.2e1;
t285 = t100 / 0.2e1;
t284 = -t101 / 0.2e1;
t283 = t101 / 0.2e1;
t279 = pkin(4) * t137;
t278 = pkin(4) * t138;
t275 = g(3) * t137;
t9 = -t71 * pkin(5) - t17;
t272 = t106 * t9;
t270 = -qJD(4) / 0.2e1;
t269 = qJD(4) / 0.2e1;
t263 = mrSges(7,3) * t156;
t262 = Ifges(7,4) * t154;
t261 = Ifges(7,4) * t156;
t258 = t101 * Ifges(5,4);
t257 = t101 * Ifges(6,6);
t241 = qJ(5) * t100;
t240 = qJDD(1) * pkin(1);
t239 = t101 * t154;
t238 = t103 * t154;
t237 = t106 * t154;
t236 = t106 * t156;
t235 = t137 * t155;
t233 = t154 * t155;
t232 = t154 * t158;
t231 = t155 * t156;
t230 = t156 * t103;
t229 = t156 * t158;
t128 = qJ(2) + t139;
t228 = qJ(5) * t235 + t155 * t278;
t227 = mrSges(4,1) * t219 + mrSges(4,2) * t218;
t226 = t158 * pkin(1) + t155 * qJ(2);
t223 = qJD(6) * t154;
t222 = qJD(6) * t156;
t217 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t69;
t215 = t86 + t323;
t214 = -m(4) - m(5) + t288;
t205 = -t223 / 0.2e1;
t141 = t158 * qJ(2);
t204 = -pkin(1) * t155 + t141;
t201 = t225 * t119;
t200 = t225 * t109;
t197 = -m(7) * t287 - mrSges(7,3);
t194 = -qJ(5) * t107 + t128;
t72 = t110 * t280 - t157 * t111;
t191 = mrSges(5,1) * t138 - mrSges(5,2) * t137;
t188 = mrSges(7,1) * t154 + mrSges(7,2) * t156;
t187 = -t138 * mrSges(6,2) + t137 * mrSges(6,3);
t186 = Ifges(7,1) * t154 + t261;
t185 = Ifges(7,2) * t156 + t262;
t184 = Ifges(7,5) * t154 + Ifges(7,6) * t156;
t183 = t11 * t156 + t12 * t154;
t42 = t106 * t287 + t194;
t50 = t107 * pkin(5) + t72;
t21 = t154 * t50 + t156 * t42;
t20 = -t154 * t42 + t156 * t50;
t179 = -t154 * t48 - t156 * t49;
t152 = -pkin(7) - qJ(3);
t178 = t158 * t139 + t155 * t152 + t204;
t176 = t155 * t139 - t152 * t158 + t226;
t66 = -mrSges(6,2) * t100 - mrSges(6,3) * t101;
t174 = -t316 - t66;
t73 = t157 * t110 + t111 * t280;
t171 = t106 * t222 + t238;
t170 = t106 * t223 - t230;
t167 = qJ(5) * t102 - qJD(5) * t107 + qJD(2);
t166 = -t179 + t320;
t45 = qJD(3) * t107 + qJD(4) * t73;
t160 = qJD(1) ^ 2;
t136 = qJDD(2) - t240;
t98 = -t138 * t233 + t229;
t97 = -t138 * t231 - t232;
t96 = -t138 * t232 - t231;
t95 = -t138 * t229 + t233;
t92 = Ifges(6,6) * t100;
t74 = Ifges(7,4) * t78;
t68 = t70 * mrSges(5,2);
t67 = t70 * mrSges(6,3);
t64 = pkin(4) * t101 + t241;
t63 = pkin(4) * t106 + t194;
t60 = -t100 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t258;
t59 = Ifges(6,4) * qJD(4) - t101 * Ifges(6,2) + t92;
t58 = Ifges(6,5) * qJD(4) + t100 * Ifges(6,3) - t257;
t54 = qJDD(4) * mrSges(5,1) + mrSges(5,3) * t70;
t51 = -t106 * pkin(5) + t73;
t41 = pkin(4) * t103 + t167;
t39 = t101 * t287 + t241;
t38 = t53 - t273;
t34 = -t102 * pkin(5) + t45;
t33 = -t103 * pkin(5) - t309;
t29 = t103 * t287 + t167;
t26 = Ifges(7,1) * t79 + Ifges(7,5) * t94 + t74;
t18 = pkin(4) * t71 + t162;
t14 = t154 * t38 + t156 * t39;
t13 = -t154 * t39 + t156 * t38;
t7 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t5 = t30 * Ifges(7,4) + t31 * Ifges(7,2) + t69 * Ifges(7,6);
t4 = -qJD(6) * t21 - t154 * t29 + t156 * t34;
t3 = qJD(6) * t20 + t154 * t34 + t156 * t29;
t6 = [(-m(3) * t204 - m(4) * t141 - m(5) * t178 - t96 * mrSges(7,1) - t95 * mrSges(7,2) + t288 * (t158 * t279 + t178) + (-m(4) * t153 - t299) * t155 + t300 * t158) * g(1) + (-m(5) * t176 - t98 * mrSges(7,1) - t97 * mrSges(7,2) + (-m(4) - m(3)) * t226 + t288 * (pkin(4) * t235 + t176) + (-m(4) * qJ(3) + t299) * t158 + t300 * t155) * g(2) + (m(5) * t22 - m(6) * t17 + t322) * t73 + (-Ifges(5,4) * t283 + Ifges(6,6) * t284 + Ifges(6,3) * t285 - Ifges(5,2) * t286 + t58 / 0.2e1 - t60 / 0.2e1 + t324 * t269 + t301) * t103 + (qJD(6) * t26 + t5) * t236 / 0.2e1 + m(4) * (qJ(2) * t114 - qJD(3) * t201 + t153 * t200) + (m(5) * t52 - t307) * t45 + (-m(5) * t53 - t306) * t309 + t304 * mrSges(6,1) + t305 * mrSges(5,3) + (t1 * t236 - t11 * t171 - t12 * t170 - t2 * t237) * mrSges(7,3) + (-t240 + t136) * mrSges(3,2) + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (m(5) * t104 + t71 * mrSges(5,1) - t68) * t128 + m(7) * (t1 * t21 + t11 * t4 + t12 * t3 + t2 * t20 + t33 * t35 + t51 * t9) + t298 * qJD(2) + (-m(5) * t23 + m(6) * t19 - t54 + t57) * t72 + (Ifges(4,1) * t151 - Ifges(4,4) * t150) * t218 + (-t109 - t330) * t203 + (-t327 / 0.2e1 - t326 / 0.2e1 - Ifges(5,1) * t283 + Ifges(6,2) * t284 + Ifges(6,6) * t285 - Ifges(5,4) * t286 - Ifges(7,5) * t290 + t59 / 0.2e1 + t325 * t269 - t297 - t321 / 0.2e1) * t102 + (t104 * mrSges(5,1) - t18 * mrSges(6,2) + t184 * t293 + t185 * t294 + t186 * t295 + t25 * t205 + (Ifges(6,3) + Ifges(5,2)) * t71 + (Ifges(6,6) + Ifges(5,4)) * t70 + t324 * t312) * t106 + t94 * (Ifges(7,5) * t171 - Ifges(7,6) * t170) / 0.2e1 + t78 * (Ifges(7,4) * t171 - Ifges(7,2) * t170) / 0.2e1 + (Ifges(7,1) * t171 - Ifges(7,4) * t170) * t290 + m(6) * (t18 * t63 + t41 * t47) - t189 * t272 + t26 * t238 / 0.2e1 + t25 * t230 / 0.2e1 + qJ(2) * t227 + m(3) * (-pkin(1) * t136 + (-t121 + t147) * qJ(2)) - (Ifges(4,4) * t151 - Ifges(4,2) * t150) * t219 + t20 * t15 + t21 * t16 + (-Ifges(5,4) * t71 + Ifges(5,5) * qJDD(4) + t217) * t107 / 0.2e1 + (t104 * mrSges(5,2) - t18 * mrSges(6,3) + Ifges(5,5) * t242 + Ifges(7,5) * t295 + Ifges(7,6) * t294 + Ifges(7,3) * t293 + (-Ifges(5,4) / 0.2e1 - Ifges(6,6)) * t71 - Ifges(6,2) * t70 - t312 * Ifges(6,4) + t303) * t107 + t35 * (mrSges(7,1) * t170 + mrSges(7,2) * t171) + t33 * t40 + t3 * t48 + t4 * t49 + t51 * t7 + t41 * t66 + t63 * (-t71 * mrSges(6,2) + t67) + t237 * t296 + t114 * t192 - 0.2e1 * t121 * mrSges(3,3) - t70 * Ifges(5,1) * t107; (-m(3) * qJ(2) - mrSges(3,3)) * t160 + (t7 + t322) * t106 + t215 * t103 + (mrSges(3,2) - t203) * qJDD(1) + t166 * t102 + (-t302 + t54) * t107 - m(5) * t305 - m(6) * t304 + m(3) * t136 + m(7) * (t183 * t102 + t103 * t35 + t107 * t308 + t272) + m(4) * t200 + (-m(6) * t47 + m(7) * t182 + t174 - t298) * qJD(1) + t315 * (m(3) - t214); -t154 * t15 + t156 * t16 + t67 - t68 + t331 * t71 + t179 * qJD(6) - t160 * t203 + t215 * t100 - t166 * t101 + t227 + (g(1) * t158 + g(2) * t155) * t214 + (t1 * t156 + t100 * t35 - t2 * t154 - t183 * t94) * m(7) + (-t46 * t100 - t43 * t101 + t18) * m(6) + (t100 * t53 - t101 * t52 + t104) * m(5) + (qJD(1) * t201 + t114) * m(4); t323 * qJD(5) + (Ifges(6,3) * t286 - t12 * t263 + t184 * t289 + t185 * t292 + t186 * t291 + t270 * t324 - t301 + t310) * t101 + t324 * t71 + (-Ifges(5,1) * t284 - Ifges(7,5) * t291 + Ifges(6,2) * t283 - Ifges(7,6) * t292 - Ifges(7,3) * t289 + t270 * t325 + t297) * t100 + t325 * t70 + (-Ifges(5,2) * t101 + t321 - t93) * t285 + (t92 + t59) * t286 - (t184 * t94 + t185 * t78 + t186 * t79) * qJD(6) / 0.2e1 + t310 * qJD(6) + (t260 + t306) * t52 + (t264 + t307) * t53 + (t7 - t56) * qJ(5) + (-pkin(4) * t19 - g(1) * t228 - qJ(5) * t17 - qJD(5) * t46 - t47 * t64) * m(6) + (-m(6) * (t127 - t279) - m(7) * t127 - t197 * t137 + (-mrSges(6,3) - t188) * t138 + t328) * g(3) + (-m(7) * t228 + (-t187 - (m(7) * pkin(8) + mrSges(7,3)) * t138 - t188 * t137) * t155) * g(1) + (-t12 * t222 - t274 + (t223 + t239) * t11) * mrSges(7,3) + (-t258 + t58) * t284 + (qJ(5) * t9 - t11 * t13 - t12 * t14 + t318 * t35) * m(7) + (t257 + t60) * t283 + (m(7) * t308 - t329) * t287 + (Ifges(5,3) + Ifges(6,1)) * qJDD(4) + (-t197 * t138 - (-m(7) * qJ(5) - t188) * t137 - m(6) * (-qJ(5) * t137 - t278) + t187 + t191) * t276 + (-t239 / 0.2e1 + t205) * t26 - t191 * t277 - t2 * t263 - t46 * t259 - t17 * mrSges(6,3) + t19 * mrSges(6,2) - t22 * mrSges(5,2) + t23 * mrSges(5,1) - t37 * t40 - t14 * t48 - t13 * t49 - pkin(4) * t57 - t64 * t66 + t43 * t265 + (Ifges(7,5) * t156 - Ifges(7,6) * t154) * t293 + (-Ifges(7,2) * t154 + t261) * t294 + (Ifges(7,1) * t156 - t262) * t295 + t156 * t296 + t9 * t188 - t154 * t5 / 0.2e1; -t323 * qJD(4) - t174 * t101 + t288 * t275 + (-qJD(4) * t35 - t101 * t182 - t308 - t319) * m(7) + (qJD(4) * t46 + t101 * t47 + t19 - t319) * m(6) + t302; -t35 * (mrSges(7,1) * t79 + mrSges(7,2) * t78) + (Ifges(7,1) * t78 - t271) * t291 + t25 * t290 + (Ifges(7,5) * t78 - Ifges(7,6) * t79) * t289 - t11 * t48 + t12 * t49 - g(1) * (mrSges(7,1) * t97 - mrSges(7,2) * t98) - g(2) * (-mrSges(7,1) * t95 + mrSges(7,2) * t96) - t189 * t275 + (t11 * t78 + t12 * t79) * mrSges(7,3) + t217 + (-Ifges(7,2) * t79 + t26 + t74) * t292 + t303;];
tau  = t6;
