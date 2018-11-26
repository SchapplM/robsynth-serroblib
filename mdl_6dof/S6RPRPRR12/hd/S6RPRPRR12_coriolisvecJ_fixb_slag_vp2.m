% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2018-11-23 16:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:10:09
% EndTime: 2018-11-23 16:10:14
% DurationCPUTime: 5.46s
% Computational Cost: add. (5404->507), mult. (11134->680), div. (0->0), fcn. (6250->6), ass. (0->241)
t181 = sin(qJ(5));
t182 = sin(qJ(3));
t185 = cos(qJ(3));
t198 = -pkin(9) * t181 * t185 - pkin(5) * t182;
t234 = qJD(5) * t181;
t186 = -pkin(3) - pkin(8);
t269 = pkin(9) - t186;
t174 = t185 * qJD(1);
t239 = qJD(1) * t182;
t140 = pkin(3) * t174 + qJ(4) * t239;
t105 = pkin(8) * t174 + t140;
t187 = -pkin(1) - pkin(7);
t156 = qJD(1) * t187 + qJD(2);
t143 = t182 * t156;
t110 = -pkin(4) * t239 + t143;
t184 = cos(qJ(5));
t55 = -t105 * t181 + t184 * t110;
t334 = -qJD(1) * t198 + t269 * t234 - t55;
t148 = t269 * t184;
t222 = t184 * t174;
t56 = t184 * t105 + t181 * t110;
t333 = pkin(9) * t222 + qJD(5) * t148 + t56;
t144 = t185 * t156;
t111 = -pkin(4) * t174 + t144;
t304 = qJD(4) - t111;
t85 = qJD(3) * t186 + t304;
t247 = qJ(4) * t185;
t206 = pkin(8) * t182 - t247;
t240 = pkin(3) * t239 + qJD(1) * qJ(2);
t95 = qJD(1) * t206 + t240;
t47 = -t181 * t95 + t184 * t85;
t48 = t181 * t85 + t184 * t95;
t204 = t181 * t47 - t184 * t48;
t266 = Ifges(6,4) * t181;
t208 = Ifges(6,2) * t184 + t266;
t265 = Ifges(6,4) * t184;
t210 = Ifges(6,1) * t181 + t265;
t213 = mrSges(6,1) * t184 - mrSges(6,2) * t181;
t261 = Ifges(6,6) * t184;
t264 = Ifges(6,5) * t181;
t277 = -t181 / 0.2e1;
t163 = t174 + qJD(5);
t278 = -t163 / 0.2e1;
t236 = qJD(3) * t184;
t136 = t181 * t239 + t236;
t282 = -t136 / 0.2e1;
t238 = qJD(3) * t181;
t135 = t184 * t239 - t238;
t283 = -t135 / 0.2e1;
t267 = Ifges(6,4) * t136;
t59 = t135 * Ifges(6,2) + t163 * Ifges(6,6) + t267;
t128 = Ifges(6,4) * t135;
t60 = t136 * Ifges(6,1) + t163 * Ifges(6,5) + t128;
t178 = qJD(3) * qJ(4);
t97 = t110 + t178;
t332 = t204 * mrSges(6,3) + (t261 + t264) * t278 + t208 * t283 + t210 * t282 + t97 * t213 + t60 * t277 - t184 * t59 / 0.2e1;
t229 = qJD(5) + qJD(6);
t331 = t185 * t229 + qJD(1);
t180 = sin(qJ(6));
t183 = cos(qJ(6));
t37 = pkin(9) * t135 + t48;
t251 = t183 * t37;
t36 = -pkin(9) * t136 + t47;
t34 = pkin(5) * t163 + t36;
t10 = t180 * t34 + t251;
t215 = t183 * t135 - t136 * t180;
t70 = t135 * t180 + t136 * t183;
t275 = Ifges(7,4) * t70;
t155 = t174 + t229;
t280 = -t155 / 0.2e1;
t290 = -t70 / 0.2e1;
t292 = -t215 / 0.2e1;
t63 = Ifges(7,4) * t215;
t31 = Ifges(7,1) * t70 + Ifges(7,5) * t155 + t63;
t61 = -pkin(5) * t135 + t97;
t252 = t180 * t37;
t9 = t183 * t34 - t252;
t330 = (Ifges(7,5) * t215 - Ifges(7,6) * t70) * t280 + (t10 * t70 + t215 * t9) * mrSges(7,3) + (-Ifges(7,2) * t70 + t31 + t63) * t292 - t61 * (mrSges(7,1) * t70 + mrSges(7,2) * t215) + (Ifges(7,1) * t215 - t275) * t290;
t270 = pkin(4) - t187;
t118 = -qJ(4) * t174 + t240;
t125 = -t143 - t178;
t316 = qJD(3) / 0.2e1;
t317 = -qJD(3) / 0.2e1;
t318 = -qJD(1) / 0.2e1;
t329 = Ifges(4,6) * t316 + (Ifges(4,4) * t185 - t182 * Ifges(4,2)) * qJD(1) / 0.2e1 + Ifges(5,5) * t317 + (-Ifges(5,6) * t185 + t182 * Ifges(5,3)) * t318 + t118 * mrSges(5,2) - t125 * mrSges(5,1) + t332;
t237 = qJD(3) * t182;
t102 = (-pkin(4) * qJD(1) + t156) * t237;
t201 = (qJD(3) * pkin(8) - qJD(4)) * t185;
t177 = qJD(1) * qJD(2);
t231 = qJD(1) * qJD(3);
t219 = t185 * t231;
t220 = t182 * t231;
t224 = pkin(3) * t219 + qJ(4) * t220 + t177;
t72 = qJD(1) * t201 + t224;
t16 = -qJD(5) * t48 + t184 * t102 - t181 * t72;
t230 = qJD(3) * qJD(5);
t233 = qJD(5) * t184;
t235 = qJD(3) * t185;
t87 = -t181 * t230 + (t181 * t235 + t182 * t233) * qJD(1);
t11 = -pkin(5) * t220 - pkin(9) * t87 + t16;
t15 = t181 * t102 + t184 * t72 + t85 * t233 - t234 * t95;
t192 = -t182 * t234 + t184 * t235;
t88 = qJD(1) * t192 - t184 * t230;
t14 = pkin(9) * t88 + t15;
t2 = qJD(6) * t9 + t11 * t180 + t14 * t183;
t25 = qJD(6) * t215 + t180 * t88 + t183 * t87;
t26 = -qJD(6) * t70 - t180 * t87 + t183 * t88;
t306 = qJD(6) * t10;
t3 = t11 * t183 - t14 * t180 - t306;
t328 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t25 + Ifges(7,6) * t26;
t327 = qJD(1) * (qJ(2) * (m(3) + m(4)) + mrSges(3,3));
t326 = qJD(4) - t144;
t243 = t183 * t184;
t244 = t180 * t181;
t199 = -t243 + t244;
t200 = t180 * t184 + t183 * t181;
t101 = t200 * t174;
t74 = t229 * t200;
t249 = -t74 - t101;
t100 = -t174 * t244 + t183 * t222;
t232 = qJD(6) * t180;
t73 = -t180 * t234 - t181 * t232 + t229 * t243;
t250 = t73 + t100;
t325 = -t10 * t250 + t199 * t3 - t2 * t200 - t249 * t9;
t324 = (-mrSges(4,1) + mrSges(5,2)) * qJD(3);
t294 = t25 / 0.2e1;
t293 = t26 / 0.2e1;
t30 = Ifges(7,2) * t215 + Ifges(7,6) * t155 + t275;
t322 = t30 / 0.2e1;
t321 = -t220 / 0.2e1;
t146 = t269 * t181;
t79 = -t146 * t183 - t148 * t180;
t315 = -qJD(6) * t79 + t333 * t180 + t334 * t183;
t78 = t146 * t180 - t148 * t183;
t314 = qJD(6) * t78 + t334 * t180 - t333 * t183;
t194 = qJD(3) * t200;
t309 = t182 * t194 + t331 * t199;
t193 = qJD(3) * t199;
t308 = -t182 * t193 + t331 * t200;
t223 = -pkin(5) * t184 - pkin(4);
t307 = pkin(5) * t233 - t174 * t223 + t326;
t106 = t199 * t182;
t305 = t182 * pkin(3) + qJ(2);
t127 = t206 + t305;
t149 = t270 * t185;
t129 = t181 * t149;
t71 = t184 * t127 + t129;
t205 = t15 * t181 + t16 * t184;
t300 = t16 * mrSges(6,1) - t15 * mrSges(6,2) + Ifges(6,5) * t87 + Ifges(6,6) * t88 + t328;
t298 = Ifges(7,4) * t294 + Ifges(7,2) * t293 + Ifges(7,6) * t321;
t297 = Ifges(7,1) * t294 + Ifges(7,4) * t293 + Ifges(7,5) * t321;
t296 = -Ifges(7,5) / 0.2e1;
t295 = -Ifges(7,6) / 0.2e1;
t291 = t215 / 0.2e1;
t289 = t70 / 0.2e1;
t288 = t87 / 0.2e1;
t287 = t88 / 0.2e1;
t285 = m(7) * t61;
t281 = t136 / 0.2e1;
t279 = t155 / 0.2e1;
t276 = t184 / 0.2e1;
t274 = pkin(9) * t182;
t271 = -Ifges(6,3) - Ifges(7,3);
t268 = Ifges(4,4) * t182;
t263 = Ifges(6,5) * t184;
t262 = Ifges(6,6) * t181;
t260 = qJ(2) * mrSges(4,1);
t259 = qJ(2) * mrSges(4,2);
t153 = mrSges(5,1) * t239 - qJD(3) * mrSges(5,3);
t76 = -mrSges(6,1) * t135 + mrSges(6,2) * t136;
t248 = t76 - t153;
t246 = qJD(3) * mrSges(4,2);
t113 = -qJD(3) * qJD(4) - t156 * t235;
t151 = -mrSges(4,3) * t239 - t246;
t242 = t151 - t153;
t241 = (mrSges(5,1) + mrSges(4,3)) * t174 + t324;
t35 = -mrSges(7,1) * t215 + mrSges(7,2) * t70;
t228 = -t35 - t248;
t227 = Ifges(5,4) / 0.2e1 - Ifges(4,5) / 0.2e1;
t226 = Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1;
t225 = -0.3e1 / 0.2e1 * Ifges(5,6) - 0.3e1 / 0.2e1 * Ifges(4,4);
t221 = pkin(3) * t235 + qJ(4) * t237 + qJD(2);
t218 = -t127 - t274;
t133 = t270 * t237;
t89 = t201 + t221;
t217 = -t184 * t133 - t181 * t89;
t116 = -qJD(3) * pkin(3) + t326;
t216 = -t116 + t144;
t134 = t270 * t235;
t212 = mrSges(6,1) * t181 + mrSges(6,2) * t184;
t211 = Ifges(6,1) * t184 - t266;
t209 = -Ifges(6,2) * t181 + t265;
t130 = t184 * t149;
t54 = pkin(5) * t185 + t181 * t218 + t130;
t57 = t184 * t274 + t71;
t27 = -t180 * t57 + t183 * t54;
t28 = t180 * t54 + t183 * t57;
t64 = -mrSges(6,1) * t220 - mrSges(6,3) * t87;
t65 = mrSges(6,2) * t220 + mrSges(6,3) * t88;
t203 = t181 * t65 + t184 * t64;
t91 = -mrSges(6,2) * t163 + mrSges(6,3) * t135;
t92 = mrSges(6,1) * t163 - mrSges(6,3) * t136;
t202 = -t181 * t92 + t184 * t91;
t90 = -pkin(4) * t219 - t113;
t32 = -t127 * t234 - t181 * t133 + t149 * t233 + t184 * t89;
t191 = -m(6) * t204 + t202;
t169 = Ifges(5,6) * t239;
t190 = t10 * mrSges(7,2) + t118 * mrSges(5,3) + t48 * mrSges(6,2) + Ifges(4,5) * t317 + (Ifges(4,1) * t185 - t268) * t318 + Ifges(5,4) * t316 - Ifges(5,2) * t174 / 0.2e1 + t169 / 0.2e1 - t155 * Ifges(7,3) - t70 * Ifges(7,5) - t215 * Ifges(7,6) - t163 * Ifges(6,3) - t136 * Ifges(6,5) - t135 * Ifges(6,6) - t47 * mrSges(6,1) - t9 * mrSges(7,1);
t170 = t182 * t187;
t164 = pkin(5) * t181 + qJ(4);
t147 = -pkin(4) * t182 + t170;
t145 = -t247 + t305;
t142 = qJD(1) * (mrSges(4,1) * t182 + mrSges(4,2) * t185);
t141 = (-mrSges(5,2) * t182 - mrSges(5,3) * t185) * qJD(1);
t112 = t182 * t223 + t170;
t109 = t200 * t185;
t108 = t200 * t182;
t107 = t199 * t185;
t103 = -qJD(4) * t185 + t221;
t86 = -qJD(4) * t174 + t224;
t75 = -pkin(5) * t192 - t134;
t68 = -t127 * t181 + t130;
t53 = mrSges(7,1) * t155 - mrSges(7,3) * t70;
t52 = -mrSges(7,2) * t155 + mrSges(7,3) * t215;
t51 = -pkin(5) * t88 + t90;
t49 = -mrSges(6,1) * t88 + mrSges(6,2) * t87;
t45 = t87 * Ifges(6,1) + t88 * Ifges(6,4) - Ifges(6,5) * t220;
t44 = t87 * Ifges(6,4) + t88 * Ifges(6,2) - Ifges(6,6) * t220;
t42 = -t182 * t74 - t185 * t193;
t40 = -t106 * t229 + t185 * t194;
t33 = -qJD(5) * t71 + t217;
t22 = pkin(9) * t192 + t32;
t21 = mrSges(7,2) * t220 + mrSges(7,3) * t26;
t20 = -mrSges(7,1) * t220 - mrSges(7,3) * t25;
t19 = t198 * qJD(3) + (t184 * t218 - t129) * qJD(5) + t217;
t13 = t183 * t36 - t252;
t12 = -t180 * t36 - t251;
t8 = -mrSges(7,1) * t26 + mrSges(7,2) * t25;
t5 = -qJD(6) * t28 - t180 * t22 + t183 * t19;
t4 = qJD(6) * t27 + t180 * t19 + t183 * t22;
t1 = [(t10 * t42 - t106 * t2 - t108 * t3 - t40 * t9) * mrSges(7,3) + (mrSges(4,1) * t177 + t181 * t45 / 0.2e1 + t44 * t276 - t90 * t213 - t86 * mrSges(5,2) + t210 * t288 + t208 * t287 + (-m(5) * t187 + mrSges(5,1)) * t113 + (t15 * t184 - t16 * t181) * mrSges(6,3) + (t59 * t277 + t60 * t276 + t97 * t212 + t135 * t209 / 0.2e1 + t211 * t281 + t163 * (-t262 + t263) / 0.2e1 + (-t181 * t48 - t184 * t47) * mrSges(6,3)) * qJD(5) + ((-m(5) * t216 + t241) * t187 + t190 + (t145 * mrSges(5,3) - 0.2e1 * t259 + t108 * t296 - t106 * t295 + (-t264 / 0.2e1 - t261 / 0.2e1 - t225) * t182 + (-0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,2) + t271) * t185) * qJD(1) + t227 * qJD(3) + t216 * mrSges(5,1)) * qJD(3)) * t182 + (Ifges(7,4) * t108 - Ifges(7,2) * t106) * t293 + (Ifges(7,1) * t108 - Ifges(7,4) * t106) * t294 + t51 * (mrSges(7,1) * t106 + mrSges(7,2) * t108) + m(6) * (-t134 * t97 + t147 * t90 + t15 * t71 + t16 * t68 + t32 * t48 + t33 * t47) + m(5) * (t103 * t118 + t145 * t86) + m(7) * (t10 * t4 + t112 * t51 + t2 * t28 + t27 * t3 + t5 * t9 + t61 * t75) + t27 * t20 + (t142 + 0.2e1 * t327) * qJD(2) + (mrSges(4,2) * t177 - t86 * mrSges(5,3) + ((-m(5) * t125 + t242) * t187 + t226 * qJD(3) + (-t145 * mrSges(5,2) + t185 * t225 + 0.2e1 * t260) * qJD(1) - t329) * qJD(3) + t300) * t185 + t42 * t322 + t108 * t297 - t106 * t298 + (Ifges(7,5) * t40 + Ifges(7,6) * t42) * t279 + (Ifges(7,1) * t40 + Ifges(7,4) * t42) * t289 + (Ifges(7,4) * t40 + Ifges(7,2) * t42) * t291 + t28 * t21 + t40 * t31 / 0.2e1 + t4 * t52 + t5 * t53 + t61 * (-mrSges(7,1) * t42 + mrSges(7,2) * t40) + t68 * t64 + t71 * t65 + t75 * t35 + t32 * t91 + t33 * t92 + t112 * t8 - t134 * t76 + t103 * t141 + t147 * t49; t107 * t20 - t109 * t21 + t308 * t53 + t309 * t52 + (t49 + t8 + (t181 * t91 + t184 * t92 + t241) * qJD(3) + m(5) * (qJD(3) * t116 - t113) + m(6) * (t236 * t47 + t238 * t48 + t90)) * t182 + (t92 * t234 - t91 * t233 + m(6) * (-t233 * t48 + t234 * t47 - t205) + (t151 + t285 + m(5) * (-t125 - t143) + m(6) * t97 - t228) * qJD(3) - t203) * t185 + (-m(5) * t118 - t141 - t142 - t191 - t327) * qJD(1) + (t10 * t309 + t107 * t3 - t109 * t2 + t51 * t182 + t308 * t9) * m(7); (-t73 / 0.2e1 - t100 / 0.2e1) * t30 + t307 * t35 + (-qJ(4) * t113 - qJD(4) * t125 - t118 * t140 + (-pkin(3) * t237 - t116 * t182 + t125 * t185) * t156) * m(5) + t248 * qJD(4) + (mrSges(7,1) * t250 + mrSges(7,2) * t249) * t61 + ((-t242 - t246) * t185 + (t324 - t241) * t182) * t156 + (qJ(4) * t90 + t186 * t205 + t304 * t97 - t47 * t55 - t48 * t56) * m(6) + t314 * t52 + (t10 * t314 + t164 * t51 + t2 * t79 + t3 * t78 + t307 * t61 + t315 * t9) * m(7) + t315 * t53 + (t186 * t191 + t332) * qJD(5) - t200 * t298 + (-Ifges(7,4) * t199 - Ifges(7,2) * t200) * t293 + (-Ifges(7,1) * t199 - Ifges(7,4) * t200) * t294 + t51 * (mrSges(7,1) * t200 - mrSges(7,2) * t199) - t199 * t297 + (-Ifges(7,5) * t74 - Ifges(7,6) * t73) * t279 + t325 * mrSges(7,3) + (((-t260 + (Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1) * t185) * qJD(1) + (-qJ(4) * mrSges(5,1) + t226) * qJD(3) + t329) * t185 + ((-Ifges(4,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1 - Ifges(5,3) / 0.2e1) * t174 + (t259 - t268 / 0.2e1) * qJD(1) + (-t263 / 0.2e1 + t262 / 0.2e1 - t199 * t296 - t200 * t295 + pkin(3) * mrSges(5,1) + t227) * qJD(3) - t169 / 0.2e1 - t190 + t116 * mrSges(5,1)) * t182) * qJD(1) + t90 * t212 + t203 * t186 - t205 * mrSges(6,3) + t45 * t276 + t44 * t277 + (Ifges(7,5) * t101 + Ifges(7,6) * t100) * t280 + t209 * t287 + t211 * t288 + (Ifges(7,1) * t101 + Ifges(7,4) * t100) * t290 + (Ifges(7,4) * t101 + Ifges(7,2) * t100) * t292 + (-Ifges(7,1) * t74 - Ifges(7,4) * t73) * t289 + (-Ifges(7,4) * t74 - Ifges(7,2) * t73) * t291 + (-t74 / 0.2e1 - t101 / 0.2e1) * t31 + qJ(4) * t49 + t78 * t20 + t79 * t21 - t56 * t91 - t55 * t92 - t111 * t76 - t113 * mrSges(5,3) - t140 * t141 + t164 * t8; t200 * t21 - t199 * t20 + t249 * t53 + t250 * t52 + t202 * qJD(5) + (t141 + t202) * t174 + ((m(5) * t156 - mrSges(5,1) * qJD(1)) * t182 + t228) * qJD(3) - m(5) * (-qJD(3) * t125 - t118 * t174) + t203 + (-qJD(3) * t61 - t325) * m(7) + (-qJD(3) * t97 - t163 * t204 + t205) * m(6); t271 * t220 + (-m(7) * t9 * t232 + (m(7) * t2 - t53 * qJD(6) + t21) * t180 + (t20 + t52 * qJD(6) + m(7) * (t3 + t306)) * t183 + (-t35 - t285) * t136) * pkin(5) + (t135 * t47 + t136 * t48) * mrSges(6,3) - m(7) * (t10 * t13 + t12 * t9) + t70 * t322 + (Ifges(6,5) * t135 - Ifges(6,6) * t136) * t278 + t59 * t281 + (Ifges(6,1) * t135 - t267) * t282 + t300 - t13 * t52 - t12 * t53 - t47 * t91 + t48 * t92 - t97 * (mrSges(6,1) * t136 + mrSges(6,2) * t135) + (-Ifges(6,2) * t136 + t128 + t60) * t283 + t330; -Ifges(7,3) * t220 + t10 * t53 + t30 * t289 - t9 * t52 + t328 + t330;];
tauc  = t1(:);
