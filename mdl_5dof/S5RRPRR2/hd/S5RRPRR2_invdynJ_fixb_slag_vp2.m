% Calculate vector of inverse dynamics joint torques for
% S5RRPRR2
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:05
% EndTime: 2019-12-05 18:27:41
% DurationCPUTime: 13.81s
% Computational Cost: add. (8523->496), mult. (20473->672), div. (0->0), fcn. (15331->16), ass. (0->228)
t230 = cos(qJ(5));
t226 = sin(qJ(5));
t223 = sin(pkin(9));
t224 = cos(pkin(9));
t228 = sin(qJ(2));
t232 = cos(qJ(2));
t175 = -t223 * t228 + t224 * t232;
t162 = t175 * qJD(1);
t274 = qJD(1) * t232;
t275 = qJD(1) * t228;
t163 = -t223 * t274 - t224 * t275;
t227 = sin(qJ(4));
t231 = cos(qJ(4));
t249 = t231 * t162 + t163 * t227;
t334 = pkin(8) * t249;
t225 = -qJ(3) - pkin(6);
t191 = t225 * t232;
t182 = qJD(1) * t191;
t166 = t223 * t182;
t189 = t225 * t228;
t181 = qJD(1) * t189;
t172 = qJD(2) * pkin(2) + t181;
t127 = t224 * t172 + t166;
t297 = pkin(7) * t163;
t91 = qJD(2) * pkin(3) + t127 + t297;
t277 = t224 * t182;
t128 = t223 * t172 - t277;
t298 = pkin(7) * t162;
t96 = t128 + t298;
t48 = t227 * t91 + t231 * t96;
t38 = t48 + t334;
t286 = t226 * t38;
t221 = qJD(2) + qJD(4);
t121 = t162 * t227 - t163 * t231;
t345 = pkin(8) * t121;
t47 = -t227 * t96 + t231 * t91;
t37 = t47 - t345;
t36 = pkin(4) * t221 + t37;
t13 = t230 * t36 - t286;
t285 = t230 * t38;
t14 = t226 * t36 + t285;
t266 = qJD(1) * qJD(2);
t255 = t228 * t266;
t265 = qJDD(1) * t232;
t183 = -t255 + t265;
t184 = qJDD(1) * t228 + t232 * t266;
t136 = t183 * t223 + t184 * t224;
t174 = t184 * pkin(6);
t271 = qJD(3) * t228;
t125 = qJDD(2) * pkin(2) - qJ(3) * t184 - qJD(1) * t271 - t174;
t207 = pkin(6) * t265;
t273 = qJD(2) * t228;
t263 = pkin(6) * t273;
t270 = qJD(3) * t232;
t131 = qJ(3) * t183 + t207 + (-t263 + t270) * qJD(1);
t76 = t224 * t125 - t131 * t223;
t50 = qJDD(2) * pkin(3) - pkin(7) * t136 + t76;
t135 = t183 * t224 - t184 * t223;
t77 = t223 * t125 + t224 * t131;
t56 = pkin(7) * t135 + t77;
t10 = -qJD(4) * t48 - t227 * t56 + t231 * t50;
t219 = qJDD(2) + qJDD(4);
t54 = qJD(4) * t249 + t135 * t227 + t136 * t231;
t6 = pkin(4) * t219 - pkin(8) * t54 + t10;
t55 = -qJD(4) * t121 + t135 * t231 - t136 * t227;
t268 = qJD(4) * t231;
t269 = qJD(4) * t227;
t9 = t227 * t50 + t231 * t56 + t91 * t268 - t269 * t96;
t7 = pkin(8) * t55 + t9;
t2 = qJD(5) * t13 + t226 * t6 + t230 * t7;
t338 = -t121 * t226 + t230 * t249;
t20 = qJD(5) * t338 + t226 * t55 + t230 * t54;
t74 = t121 * t230 + t226 * t249;
t21 = -qJD(5) * t74 - t226 * t54 + t230 * t55;
t212 = qJDD(5) + t219;
t215 = qJD(5) + t221;
t3 = -qJD(5) * t14 - t226 * t7 + t230 * t6;
t301 = Ifges(6,4) * t74;
t311 = t74 / 0.2e1;
t32 = Ifges(6,2) * t338 + Ifges(6,6) * t215 + t301;
t68 = Ifges(6,4) * t338;
t33 = Ifges(6,1) * t74 + Ifges(6,5) * t215 + t68;
t218 = t232 * pkin(2);
t206 = t218 + pkin(1);
t186 = -qJD(1) * t206 + qJD(3);
t137 = -pkin(3) * t162 + t186;
t82 = -pkin(4) * t249 + t137;
t362 = -(Ifges(6,1) * t338 - t301) * t74 / 0.2e1 - (Ifges(6,5) * t338 - Ifges(6,6) * t74) * t215 / 0.2e1 + (t13 * t338 + t14 * t74) * mrSges(6,3) - t82 * (mrSges(6,1) * t74 + mrSges(6,2) * t338) + t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t20 + Ifges(6,6) * t21 + Ifges(6,3) * t212 + t32 * t311 - (-Ifges(6,2) * t74 + t33 + t68) * t338 / 0.2e1;
t113 = Ifges(5,4) * t249;
t289 = Ifges(5,4) * t121;
t67 = Ifges(5,1) * t121 + Ifges(5,5) * t221 + t113;
t361 = t10 * mrSges(5,1) - t9 * mrSges(5,2) + Ifges(5,5) * t54 + Ifges(5,6) * t55 + Ifges(5,3) * t219 - (Ifges(5,5) * t249 - Ifges(5,6) * t121) * t221 / 0.2e1 + (t121 * t48 + t249 * t47) * mrSges(5,3) - (-Ifges(5,2) * t121 + t113 + t67) * t249 / 0.2e1 - t137 * (mrSges(5,1) * t121 + mrSges(5,2) * t249) - (Ifges(5,1) * t249 - t289) * t121 / 0.2e1 + t362;
t204 = pkin(2) * t224 + pkin(3);
t300 = pkin(2) * t223;
t159 = t204 * t227 + t231 * t300;
t133 = -t181 * t223 + t277;
t98 = t133 - t298;
t134 = t224 * t181 + t166;
t99 = t134 + t297;
t342 = -t159 * qJD(4) + t227 * t99 - t231 * t98;
t158 = t231 * t204 - t227 * t300;
t341 = t158 * qJD(4) - t227 * t98 - t231 * t99;
t358 = t345 + t341;
t357 = t334 + t342;
t293 = mrSges(3,2) * t232;
t355 = mrSges(3,1) * t228 + t293;
t303 = t228 / 0.2e1;
t190 = -t232 * mrSges(3,1) + t228 * mrSges(3,2);
t222 = qJ(2) + pkin(9);
t213 = sin(t222);
t214 = cos(t222);
t216 = qJ(4) + t222;
t202 = sin(t216);
t203 = cos(t216);
t205 = qJ(5) + t216;
t198 = sin(t205);
t199 = cos(t205);
t251 = t199 * mrSges(6,1) - t198 * mrSges(6,2);
t240 = -t203 * mrSges(5,1) + t202 * mrSges(5,2) - t251;
t348 = mrSges(4,1) * t214 - mrSges(4,2) * t213 - t190 - t240;
t156 = pkin(4) + t158;
t107 = t156 * t226 + t159 * t230;
t344 = -qJD(5) * t107 - t226 * t358 + t230 * t357;
t106 = t156 * t230 - t159 * t226;
t343 = qJD(5) * t106 + t226 * t357 + t230 * t358;
t229 = sin(qJ(1));
t233 = cos(qJ(1));
t319 = g(1) * t233 + g(2) * t229;
t66 = Ifges(5,2) * t249 + Ifges(5,6) * t221 + t289;
t335 = t66 / 0.2e1;
t282 = qJDD(1) * pkin(1);
t138 = t224 * t189 + t191 * t223;
t176 = t223 * t232 + t224 * t228;
t111 = -pkin(7) * t176 + t138;
t139 = t223 * t189 - t224 * t191;
t112 = pkin(7) * t175 + t139;
t65 = t227 * t111 + t231 * t112;
t173 = -pkin(6) * t255 + t207;
t321 = t173 * t232 + t174 * t228;
t276 = pkin(3) * t214 + t218;
t261 = pkin(4) * t203 + t276;
t318 = mrSges(2,1) + m(6) * (pkin(1) + t261) + m(5) * (pkin(1) + t276) + m(4) * t206 + m(3) * pkin(1) + t348;
t220 = -pkin(7) + t225;
t317 = mrSges(2,2) + m(6) * (-pkin(8) + t220) - mrSges(6,3) + m(5) * t220 - mrSges(5,3) + m(4) * t225 - mrSges(4,3) - m(3) * pkin(6) - mrSges(3,3);
t316 = m(4) * pkin(2);
t308 = t121 / 0.2e1;
t306 = -t163 / 0.2e1;
t299 = pkin(2) * t228;
t294 = -qJD(2) / 0.2e1;
t292 = Ifges(3,4) * t228;
t291 = Ifges(3,4) * t232;
t290 = Ifges(4,4) * t163;
t284 = t232 * Ifges(3,2);
t254 = qJD(2) * t225;
t160 = t228 * t254 + t270;
t161 = t232 * t254 - t271;
t110 = t224 * t160 + t223 * t161;
t272 = qJD(2) * t232;
t210 = pkin(2) * t273;
t257 = -t55 * mrSges(5,1) + t54 * mrSges(5,2);
t256 = -t21 * mrSges(6,1) + t20 * mrSges(6,2);
t140 = pkin(2) * t275 - pkin(3) * t163;
t164 = t176 * qJD(2);
t141 = pkin(3) * t164 + t210;
t252 = -t135 * mrSges(4,1) + t136 * mrSges(4,2);
t64 = t231 * t111 - t112 * t227;
t109 = -t160 * t223 + t224 * t161;
t144 = -pkin(3) * t175 - t206;
t185 = -pkin(3) * t213 - t299;
t247 = mrSges(6,1) * t198 + mrSges(6,2) * t199;
t246 = t284 + t292;
t245 = Ifges(3,5) * t232 - Ifges(3,6) * t228;
t130 = t175 * t227 + t176 * t231;
t41 = -pkin(8) * t130 + t64;
t129 = t175 * t231 - t176 * t227;
t42 = pkin(8) * t129 + t65;
t26 = -t226 * t42 + t230 * t41;
t27 = t226 * t41 + t230 * t42;
t80 = t129 * t230 - t130 * t226;
t81 = t129 * t226 + t130 * t230;
t243 = pkin(1) * t355;
t155 = -pkin(2) * t183 + qJDD(3) - t282;
t242 = t228 * (Ifges(3,1) * t232 - t292);
t165 = t175 * qJD(2);
t88 = -pkin(7) * t165 + t109;
t89 = -pkin(7) * t164 + t110;
t30 = t111 * t268 - t112 * t269 + t227 * t88 + t231 * t89;
t97 = -pkin(3) * t135 + t155;
t236 = mrSges(5,1) * t202 + mrSges(5,2) * t203 + t247;
t31 = -qJD(4) * t65 - t227 * t89 + t231 * t88;
t208 = Ifges(3,4) * t274;
t188 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t274;
t187 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t275;
t171 = Ifges(3,1) * t275 + Ifges(3,5) * qJD(2) + t208;
t170 = Ifges(3,6) * qJD(2) + qJD(1) * t246;
t157 = Ifges(4,4) * t162;
t143 = qJD(2) * mrSges(4,1) + mrSges(4,3) * t163;
t142 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t162;
t126 = -mrSges(4,1) * t162 - mrSges(4,2) * t163;
t123 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t136;
t122 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t135;
t117 = -t163 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t157;
t116 = t162 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t290;
t101 = mrSges(5,1) * t221 - mrSges(5,3) * t121;
t100 = -mrSges(5,2) * t221 + mrSges(5,3) * t249;
t92 = -pkin(4) * t129 + t144;
t83 = pkin(4) * t121 + t140;
t79 = -qJD(4) * t130 - t164 * t231 - t165 * t227;
t78 = qJD(4) * t129 - t164 * t227 + t165 * t231;
t75 = -mrSges(5,1) * t249 + mrSges(5,2) * t121;
t61 = mrSges(6,1) * t215 - mrSges(6,3) * t74;
t60 = -mrSges(6,2) * t215 + mrSges(6,3) * t338;
t59 = -pkin(4) * t79 + t141;
t44 = -mrSges(5,2) * t219 + mrSges(5,3) * t55;
t43 = mrSges(5,1) * t219 - mrSges(5,3) * t54;
t35 = -pkin(4) * t55 + t97;
t34 = -mrSges(6,1) * t338 + mrSges(6,2) * t74;
t29 = -qJD(5) * t81 - t226 * t78 + t230 * t79;
t28 = qJD(5) * t80 + t226 * t79 + t230 * t78;
t25 = -pkin(8) * t78 + t31;
t24 = pkin(8) * t79 + t30;
t16 = t230 * t37 - t286;
t15 = -t226 * t37 - t285;
t12 = -mrSges(6,2) * t212 + mrSges(6,3) * t21;
t11 = mrSges(6,1) * t212 - mrSges(6,3) * t20;
t5 = -qJD(5) * t27 - t226 * t24 + t230 * t25;
t4 = qJD(5) * t26 + t226 * t25 + t230 * t24;
t1 = [(mrSges(6,2) * t35 - mrSges(6,3) * t3 + Ifges(6,1) * t20 + Ifges(6,4) * t21 + Ifges(6,5) * t212) * t81 + qJD(2) ^ 2 * t245 / 0.2e1 + t183 * t246 / 0.2e1 + ((t232 * t183 + t184 * t228) * pkin(6) + t321) * mrSges(3,3) + t126 * t210 + t29 * t32 / 0.2e1 + t28 * t33 / 0.2e1 + t26 * t11 + t27 * t12 + (Ifges(5,1) * t78 + Ifges(5,4) * t79) * t308 + (Ifges(6,1) * t28 + Ifges(6,4) * t29) * t311 + (0.2e1 * Ifges(3,5) * t303 + Ifges(4,5) * t176 + Ifges(3,6) * t232 + Ifges(4,6) * t175 - pkin(6) * t355) * qJDD(2) + (Ifges(3,1) * t184 + Ifges(3,4) * t183) * t303 + (m(3) * t321 - t187 * t272) * pkin(6) + t232 * (Ifges(3,4) * t184 + Ifges(3,2) * t183) / 0.2e1 + (mrSges(5,2) * t97 - mrSges(5,3) * t10 + Ifges(5,1) * t54 + Ifges(5,4) * t55 + Ifges(5,5) * t219) * t130 + (-t13 * t28 + t14 * t29) * mrSges(6,3) + t338 * (Ifges(6,4) * t28 + Ifges(6,2) * t29) / 0.2e1 + (-mrSges(6,1) * t35 + mrSges(6,3) * t2 + Ifges(6,4) * t20 + Ifges(6,2) * t21 + Ifges(6,6) * t212) * t80 + (t232 * (-Ifges(3,2) * t228 + t291) + t242) * t266 / 0.2e1 + t249 * (Ifges(5,4) * t78 + Ifges(5,2) * t79) / 0.2e1 + (Ifges(4,1) * t165 - Ifges(4,4) * t164) * t306 + (-t127 * t165 - t128 * t164 + t175 * t77 - t176 * t76) * mrSges(4,3) + t162 * (Ifges(4,4) * t165 - Ifges(4,2) * t164) / 0.2e1 + qJD(2) * (Ifges(4,5) * t165 - Ifges(4,6) * t164) / 0.2e1 + t186 * (mrSges(4,1) * t164 + mrSges(4,2) * t165) + t59 * t34 + t4 * t60 + t5 * t61 + t64 * t43 + t65 * t44 + t78 * t67 / 0.2e1 - t206 * t252 + t82 * (-mrSges(6,1) * t29 + mrSges(6,2) * t28) + t30 * t100 + t92 * t256 + t144 * t257 + t31 * t101 + t79 * t335 + m(4) * (t109 * t127 + t110 * t128 + t138 * t76 + t139 * t77 - t155 * t206 + t186 * t210) - t188 * t263 - t243 * t266 + t137 * (-mrSges(5,1) * t79 + mrSges(5,2) * t78) + t138 * t123 + t139 * t122 + t141 * t75 + t171 * t272 / 0.2e1 + (-t47 * t78 + t48 * t79) * mrSges(5,3) - t170 * t273 / 0.2e1 + t110 * t142 + t109 * t143 - t164 * t116 / 0.2e1 + t165 * t117 / 0.2e1 - t190 * t282 + t135 * (Ifges(4,4) * t176 + Ifges(4,2) * t175) + t136 * (Ifges(4,1) * t176 + Ifges(4,4) * t175) + t155 * (-mrSges(4,1) * t175 + mrSges(4,2) * t176) + t184 * (t228 * Ifges(3,1) + t291) / 0.2e1 + t215 * (Ifges(6,5) * t28 + Ifges(6,6) * t29) / 0.2e1 + t221 * (Ifges(5,5) * t78 + Ifges(5,6) * t79) / 0.2e1 + m(5) * (t10 * t64 + t137 * t141 + t144 * t97 + t30 * t48 + t31 * t47 + t65 * t9) + m(6) * (t13 * t5 + t14 * t4 + t2 * t27 + t26 * t3 + t35 * t92 + t59 * t82) + Ifges(2,3) * qJDD(1) + (-mrSges(5,1) * t97 + mrSges(5,3) * t9 + Ifges(5,4) * t54 + Ifges(5,2) * t55 + Ifges(5,6) * t219) * t129 + (t229 * t318 + t233 * t317) * g(1) + (t229 * t317 - t233 * t318) * g(2) + (m(3) * t282 + mrSges(3,1) * t183 - mrSges(3,2) * t184) * pkin(1); -m(4) * (t127 * t133 + t128 * t134) + (t122 * t223 + t123 * t224) * pkin(2) + (t223 * t77 + t224 * t76) * t316 + (Ifges(4,5) * t162 + Ifges(4,6) * t163) * t294 + t116 * t306 + (t127 * t162 - t128 * t163) * mrSges(4,3) + (-m(4) * t218 - t348) * g(3) + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + t361 + (t170 * t303 + t245 * t294 + (t284 * t303 - t242 / 0.2e1 + t243) * qJD(1) + (t232 * t187 + t228 * t188) * pkin(6) + (-m(4) * t186 - t126) * t299 - (t171 + t208) * t232 / 0.2e1) * qJD(1) - (Ifges(4,2) * t163 + t117 + t157) * t162 / 0.2e1 + t121 * t335 + t76 * mrSges(4,1) - t77 * mrSges(4,2) + t319 * (-m(5) * t185 - m(6) * (-pkin(4) * t202 + t185) + mrSges(4,1) * t213 + t293 + mrSges(4,2) * t214 + (mrSges(3,1) + t316) * t228 + t236) - t83 * t34 + t341 * t100 + t342 * t101 + (-g(3) * t276 + t10 * t158 - t137 * t140 + t159 * t9 + t341 * t48 + t342 * t47) * m(5) + t106 * t11 + t107 * t12 + t343 * t60 + t344 * t61 + (-g(3) * t261 + t106 * t3 + t107 * t2 + t344 * t13 + t343 * t14 - t82 * t83) * m(6) + Ifges(4,6) * t135 + Ifges(4,5) * t136 - t140 * t75 - t134 * t142 - t133 * t143 + t158 * t43 + t159 * t44 - t173 * mrSges(3,2) - t174 * mrSges(3,1) + Ifges(3,6) * t183 + Ifges(3,5) * t184 - t186 * (-mrSges(4,1) * t163 + mrSges(4,2) * t162) + t163 * (Ifges(4,1) * t162 + t290) / 0.2e1; -t249 * t100 + t121 * t101 - t162 * t142 - t163 * t143 - t338 * t60 + t74 * t61 + t252 + t256 + t257 + (-g(1) * t229 + g(2) * t233) * (m(4) + m(5) + m(6)) + (t13 * t74 - t14 * t338 + t35) * m(6) + (t121 * t47 - t249 * t48 + t97) * m(5) + (-t127 * t163 - t128 * t162 + t155) * m(4); t319 * t236 - m(6) * (t13 * t15 + t14 * t16) + t66 * t308 + t240 * g(3) - t16 * t60 - t15 * t61 - t47 * t100 + t48 * t101 + (t230 * t11 + t226 * t12 - t121 * t34 + (-g(3) * t203 - t121 * t82 + t2 * t226 + t319 * t202 + t230 * t3) * m(6) + (-t226 * t61 + t230 * t60 + (-t13 * t226 + t14 * t230) * m(6)) * qJD(5)) * pkin(4) + t361; -g(3) * t251 - t13 * t60 + t14 * t61 + t319 * t247 + t362;];
tau = t1;
