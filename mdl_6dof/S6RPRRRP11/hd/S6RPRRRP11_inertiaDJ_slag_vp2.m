% Calculate time derivative of joint inertia matrix for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP11_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:34:38
% EndTime: 2019-03-09 06:34:54
% DurationCPUTime: 8.11s
% Computational Cost: add. (13344->667), mult. (39072->928), div. (0->0), fcn. (41025->12), ass. (0->266)
t223 = cos(qJ(5));
t338 = Ifges(7,6) + Ifges(6,6);
t341 = t338 * t223;
t220 = sin(qJ(5));
t340 = (Ifges(6,5) + Ifges(7,5)) * t220;
t215 = sin(pkin(7));
t218 = cos(pkin(7));
t219 = cos(pkin(6));
t214 = sin(pkin(12));
t216 = sin(pkin(6));
t217 = cos(pkin(12));
t289 = t216 * t217;
t315 = pkin(1) * t219;
t276 = qJ(2) * t289 + t214 * t315;
t123 = (t215 * t219 + t218 * t289) * pkin(9) + t276;
t222 = sin(qJ(3));
t225 = cos(qJ(3));
t203 = t217 * t315;
t293 = t214 * t216;
t127 = pkin(2) * t219 + t203 + (-pkin(9) * t218 - qJ(2)) * t293;
t139 = (-pkin(9) * t214 * t215 - pkin(2) * t217 - pkin(1)) * t216;
t232 = t127 * t218 + t139 * t215;
t81 = -t222 * t123 + t232 * t225;
t224 = cos(qJ(4));
t271 = qJD(4) * t224;
t250 = t223 * t271;
t221 = sin(qJ(4));
t269 = qJD(5) * t221;
t253 = t220 * t269;
t228 = t250 - t253;
t251 = t220 * t271;
t268 = qJD(5) * t223;
t227 = t221 * t268 + t251;
t287 = t218 * t225;
t290 = t215 * t225;
t337 = t216 * (-t214 * t222 + t217 * t287) + t219 * t290;
t336 = -2 * Ifges(4,4);
t288 = t218 * t222;
t291 = t215 * t222;
t126 = t219 * t291 + (t214 * t225 + t217 * t288) * t216;
t120 = t126 * qJD(3);
t152 = -t215 * t289 + t218 * t219;
t101 = t126 * t221 - t152 * t224;
t119 = t337 * qJD(3);
t88 = -qJD(4) * t101 + t119 * t224;
t102 = t126 * t224 + t152 * t221;
t90 = t102 * t223 - t220 * t337;
t52 = -qJD(5) * t90 + t120 * t223 - t220 * t88;
t89 = -t102 * t220 - t223 * t337;
t53 = qJD(5) * t89 + t120 * t220 + t223 * t88;
t87 = qJD(4) * t102 + t119 * t221;
t8 = Ifges(7,5) * t53 + Ifges(7,6) * t52 + Ifges(7,3) * t87;
t9 = Ifges(6,5) * t53 + Ifges(6,6) * t52 + Ifges(6,3) * t87;
t335 = t8 + t9;
t331 = -m(6) * pkin(11) - mrSges(6,3);
t180 = -pkin(4) * t224 - pkin(11) * t221 - pkin(3);
t284 = t223 * t224;
t204 = pkin(10) * t284;
t142 = t220 * t180 + t204;
t183 = -mrSges(6,1) * t223 + mrSges(6,2) * t220;
t329 = -m(6) * pkin(4) - mrSges(5,1) + t183;
t182 = -mrSges(7,1) * t223 + mrSges(7,2) * t220;
t206 = -pkin(5) * t223 - pkin(4);
t328 = m(7) * t206 + t182;
t327 = 2 * m(5);
t326 = 0.2e1 * m(6);
t325 = 0.2e1 * m(7);
t324 = 0.2e1 * pkin(10);
t323 = -2 * mrSges(4,3);
t322 = -2 * mrSges(7,3);
t321 = m(6) / 0.2e1;
t320 = m(5) * pkin(3);
t317 = m(7) * pkin(5);
t272 = qJD(4) * t221;
t275 = qJD(2) * t216;
t69 = (-t214 * t288 + t217 * t225) * t275 + t81 * qJD(3);
t100 = -t127 * t215 + t218 * t139;
t74 = -pkin(3) * t337 - pkin(10) * t126 + t100;
t115 = t225 * t123;
t82 = t127 * t288 + t139 * t291 + t115;
t79 = pkin(10) * t152 + t82;
t256 = t214 * t275;
t238 = t215 * t256;
t94 = pkin(3) * t120 - pkin(10) * t119 + t238;
t20 = -t221 * t69 + t224 * t94 - t79 * t271 - t74 * t272;
t18 = -pkin(4) * t120 - t20;
t316 = m(6) * t18;
t314 = pkin(10) * t220;
t312 = -qJ(6) - pkin(11);
t22 = -mrSges(6,1) * t52 + mrSges(6,2) * t53;
t66 = mrSges(5,1) * t120 - mrSges(5,3) * t88;
t311 = t22 - t66;
t35 = t221 * t74 + t224 * t79;
t33 = -pkin(11) * t337 + t35;
t78 = -t152 * pkin(3) - t81;
t51 = t101 * pkin(4) - t102 * pkin(11) + t78;
t15 = t220 * t51 + t223 * t33;
t58 = -mrSges(6,1) * t89 + mrSges(6,2) * t90;
t92 = -mrSges(5,1) * t337 - mrSges(5,3) * t102;
t310 = t58 - t92;
t309 = Ifges(5,4) * t221;
t308 = Ifges(5,4) * t224;
t307 = Ifges(6,4) * t220;
t306 = Ifges(6,4) * t223;
t305 = Ifges(7,4) * t220;
t304 = Ifges(7,4) * t223;
t303 = t120 * Ifges(5,5);
t302 = t120 * Ifges(5,6);
t301 = t337 * Ifges(5,6);
t300 = t221 * mrSges(5,2);
t70 = (t214 * t287 + t217 * t222) * t275 + (t222 * t232 + t115) * qJD(3);
t299 = t225 * t70;
t298 = -t224 * mrSges(5,1) - mrSges(4,1) + t300;
t153 = -t224 * t218 + t221 * t291;
t273 = qJD(3) * t225;
t254 = t215 * t273;
t129 = -qJD(4) * t153 + t224 * t254;
t296 = t129 * t224;
t154 = t218 * t221 + t224 * t291;
t130 = qJD(4) * t154 + t221 * t254;
t295 = t130 * t221;
t294 = t153 * t130;
t286 = t220 * t221;
t285 = t221 * t223;
t283 = Ifges(4,5) * t119 - Ifges(4,6) * t120;
t233 = -Ifges(7,2) * t220 + t304;
t146 = -Ifges(7,6) * t224 + t221 * t233;
t234 = -Ifges(6,2) * t220 + t306;
t147 = -Ifges(6,6) * t224 + t221 * t234;
t282 = -t146 - t147;
t235 = Ifges(7,1) * t223 - t305;
t148 = -Ifges(7,5) * t224 + t221 * t235;
t236 = Ifges(6,1) * t223 - t307;
t149 = -Ifges(6,5) * t224 + t221 * t236;
t281 = t148 + t149;
t178 = (pkin(4) * t221 - pkin(11) * t224) * qJD(4);
t280 = t220 * t178 + t180 * t268;
t279 = t223 * t178 + t272 * t314;
t278 = Ifges(7,5) * t250 + Ifges(7,3) * t272;
t277 = Ifges(6,5) * t250 + Ifges(6,3) * t272;
t270 = qJD(5) * t220;
t162 = mrSges(7,1) * t270 + mrSges(7,2) * t268;
t274 = qJD(3) * t222;
t267 = qJD(6) * t223;
t265 = Ifges(5,5) * t88 - Ifges(5,6) * t87 + Ifges(5,3) * t120;
t264 = pkin(5) * t270;
t10 = Ifges(7,4) * t53 + Ifges(7,2) * t52 + Ifges(7,6) * t87;
t11 = Ifges(6,4) * t53 + Ifges(6,2) * t52 + Ifges(6,6) * t87;
t263 = -t10 / 0.2e1 - t11 / 0.2e1;
t12 = Ifges(7,1) * t53 + Ifges(7,4) * t52 + Ifges(7,5) * t87;
t13 = Ifges(6,1) * t53 + Ifges(6,4) * t52 + Ifges(6,5) * t87;
t262 = t12 / 0.2e1 + t13 / 0.2e1;
t39 = Ifges(7,4) * t90 + Ifges(7,2) * t89 + Ifges(7,6) * t101;
t40 = Ifges(6,4) * t90 + Ifges(6,2) * t89 + Ifges(6,6) * t101;
t261 = -t39 / 0.2e1 - t40 / 0.2e1;
t41 = Ifges(7,1) * t90 + Ifges(7,4) * t89 + Ifges(7,5) * t101;
t42 = Ifges(6,1) * t90 + Ifges(6,4) * t89 + Ifges(6,5) * t101;
t260 = t41 / 0.2e1 + t42 / 0.2e1;
t259 = mrSges(7,1) + t317;
t255 = t215 * t274;
t188 = Ifges(7,2) * t223 + t305;
t109 = -t188 * t269 + (Ifges(7,6) * t221 + t224 * t233) * qJD(4);
t189 = Ifges(6,2) * t223 + t307;
t110 = -t189 * t269 + (Ifges(6,6) * t221 + t224 * t234) * qJD(4);
t249 = t109 / 0.2e1 + t110 / 0.2e1;
t191 = Ifges(7,1) * t220 + t304;
t111 = -t191 * t269 + (Ifges(7,5) * t221 + t224 * t235) * qJD(4);
t192 = Ifges(6,1) * t220 + t306;
t112 = -t192 * t269 + (Ifges(6,5) * t221 + t224 * t236) * qJD(4);
t248 = t111 / 0.2e1 + t112 / 0.2e1;
t247 = t146 / 0.2e1 + t147 / 0.2e1;
t246 = t148 / 0.2e1 + t149 / 0.2e1;
t211 = Ifges(7,5) * t268;
t212 = Ifges(6,5) * t268;
t245 = t211 / 0.2e1 + t212 / 0.2e1 - t338 * t270 / 0.2e1;
t167 = t233 * qJD(5);
t168 = t234 * qJD(5);
t244 = t167 / 0.2e1 + t168 / 0.2e1;
t170 = t235 * qJD(5);
t171 = t236 * qJD(5);
t243 = t170 / 0.2e1 + t171 / 0.2e1;
t242 = t340 / 0.2e1 + t341 / 0.2e1;
t241 = t188 / 0.2e1 + t189 / 0.2e1;
t240 = -t192 / 0.2e1 - t191 / 0.2e1;
t21 = -t52 * mrSges(7,1) + t53 * mrSges(7,2);
t14 = -t220 * t33 + t223 * t51;
t34 = -t221 * t79 + t224 * t74;
t239 = qJD(5) * t312;
t237 = mrSges(6,1) * t220 + mrSges(6,2) * t223;
t32 = pkin(4) * t337 - t34;
t131 = -t220 * t154 - t223 * t290;
t231 = -t223 * t154 + t220 * t290;
t19 = t221 * t94 + t224 * t69 + t74 * t271 - t272 * t79;
t17 = pkin(11) * t120 + t19;
t30 = t87 * pkin(4) - t88 * pkin(11) + t70;
t3 = t223 * t17 + t220 * t30 + t51 * t268 - t270 * t33;
t230 = t153 * t271 + t295;
t121 = mrSges(7,1) * t227 + mrSges(7,2) * t228;
t4 = -qJD(5) * t15 - t17 * t220 + t223 * t30;
t213 = Ifges(5,5) * t271;
t193 = Ifges(5,1) * t221 + t308;
t190 = Ifges(5,2) * t224 + t309;
t185 = t312 * t223;
t181 = t312 * t220;
t179 = (pkin(5) * t220 + pkin(10)) * t221;
t177 = -mrSges(6,1) * t224 - mrSges(6,3) * t285;
t176 = -mrSges(7,1) * t224 - mrSges(7,3) * t285;
t175 = mrSges(6,2) * t224 - mrSges(6,3) * t286;
t174 = mrSges(7,2) * t224 - mrSges(7,3) * t286;
t172 = (Ifges(5,1) * t224 - t309) * qJD(4);
t169 = (-Ifges(5,2) * t221 + t308) * qJD(4);
t164 = (mrSges(5,1) * t221 + mrSges(5,2) * t224) * qJD(4);
t163 = t237 * qJD(5);
t161 = t223 * t180;
t158 = t237 * t221;
t157 = (mrSges(7,1) * t220 + mrSges(7,2) * t223) * t221;
t151 = -qJD(6) * t220 + t223 * t239;
t150 = t220 * t239 + t267;
t145 = -Ifges(6,3) * t224 + (Ifges(6,5) * t223 - Ifges(6,6) * t220) * t221;
t144 = -Ifges(7,3) * t224 + (Ifges(7,5) * t223 - Ifges(7,6) * t220) * t221;
t141 = -t224 * t314 + t161;
t140 = pkin(5) * t227 + pkin(10) * t271;
t138 = -mrSges(6,2) * t272 - mrSges(6,3) * t227;
t137 = -mrSges(7,2) * t272 - mrSges(7,3) * t227;
t136 = mrSges(6,1) * t272 - mrSges(6,3) * t228;
t135 = mrSges(7,1) * t272 - mrSges(7,3) * t228;
t133 = -qJ(6) * t286 + t142;
t124 = -qJ(6) * t285 + t161 + (-pkin(5) - t314) * t224;
t122 = mrSges(6,1) * t227 + mrSges(6,2) * t228;
t108 = -Ifges(6,5) * t253 - Ifges(6,6) * t227 + t277;
t107 = -Ifges(7,5) * t253 - Ifges(7,6) * t227 + t278;
t106 = -qJD(5) * t142 + t279;
t105 = (-t223 * t272 - t224 * t270) * pkin(10) + t280;
t104 = mrSges(4,1) * t152 - mrSges(4,3) * t126;
t103 = -mrSges(4,2) * t152 + mrSges(4,3) * t337;
t99 = (-pkin(10) * qJD(4) - qJ(6) * qJD(5)) * t285 + (-qJD(6) * t221 + (-pkin(10) * qJD(5) - qJ(6) * qJD(4)) * t224) * t220 + t280;
t98 = qJD(5) * t231 - t220 * t129 + t223 * t255;
t97 = qJD(5) * t131 + t223 * t129 + t220 * t255;
t96 = mrSges(4,1) * t120 + mrSges(4,2) * t119;
t95 = -t221 * t267 + (pkin(5) * t221 - qJ(6) * t284) * qJD(4) + (-t204 + (qJ(6) * t221 - t180) * t220) * qJD(5) + t279;
t91 = mrSges(5,2) * t337 - mrSges(5,3) * t101;
t80 = mrSges(5,1) * t101 + mrSges(5,2) * t102;
t65 = -mrSges(5,2) * t120 - mrSges(5,3) * t87;
t64 = Ifges(5,1) * t102 - Ifges(5,4) * t101 - Ifges(5,5) * t337;
t63 = Ifges(5,4) * t102 - Ifges(5,2) * t101 - t301;
t62 = mrSges(6,1) * t101 - mrSges(6,3) * t90;
t61 = mrSges(7,1) * t101 - mrSges(7,3) * t90;
t60 = -mrSges(6,2) * t101 + mrSges(6,3) * t89;
t59 = -mrSges(7,2) * t101 + mrSges(7,3) * t89;
t57 = -mrSges(7,1) * t89 + mrSges(7,2) * t90;
t56 = mrSges(5,1) * t87 + mrSges(5,2) * t88;
t55 = Ifges(5,1) * t88 - Ifges(5,4) * t87 + t303;
t54 = Ifges(5,4) * t88 - Ifges(5,2) * t87 + t302;
t38 = Ifges(6,5) * t90 + Ifges(6,6) * t89 + Ifges(6,3) * t101;
t37 = Ifges(7,5) * t90 + Ifges(7,6) * t89 + Ifges(7,3) * t101;
t27 = mrSges(6,1) * t87 - mrSges(6,3) * t53;
t26 = mrSges(7,1) * t87 - mrSges(7,3) * t53;
t25 = -mrSges(6,2) * t87 + mrSges(6,3) * t52;
t24 = -mrSges(7,2) * t87 + mrSges(7,3) * t52;
t23 = -pkin(5) * t89 + t32;
t7 = qJ(6) * t89 + t15;
t6 = pkin(5) * t101 - qJ(6) * t90 + t14;
t5 = -pkin(5) * t52 + t18;
t2 = qJ(6) * t52 + qJD(6) * t89 + t3;
t1 = pkin(5) * t87 - qJ(6) * t53 - qJD(6) * t90 + t4;
t16 = [0.2e1 * (t217 * (-mrSges(3,2) * t219 + mrSges(3,3) * t289) + m(3) * (t276 * t217 + (qJ(2) * t293 - t203) * t214)) * t275 + (-t54 + t335) * t101 + (t12 + t13) * t90 + (t10 + t11) * t89 + (t37 + t38 - t63) * t87 + (t42 + t41) * t53 + (t39 + t40) * t52 + (t126 * t336 + Ifges(5,5) * t102 - Ifges(4,6) * t152 - Ifges(5,6) * t101 + t323 * t82 - ((2 * Ifges(4,2)) + Ifges(5,3)) * t337) * t120 + (0.2e1 * Ifges(4,1) * t126 + Ifges(4,5) * t152 + t323 * t81 - t336 * t337) * t119 + 0.2e1 * (-mrSges(4,1) * t337 + mrSges(4,2) * t126) * t238 - t337 * t265 + (t1 * t6 + t2 * t7 + t23 * t5) * t325 + (t14 * t4 + t15 * t3 + t18 * t32) * t326 + (t19 * t35 + t20 * t34 + t70 * t78) * t327 + 0.2e1 * t69 * t103 - 0.2e1 * t70 * t104 + t102 * t55 + 0.2e1 * t100 * t96 + 0.2e1 * t19 * t91 + 0.2e1 * t20 * t92 + t88 * t64 + 0.2e1 * t78 * t56 + 0.2e1 * t70 * t80 + 0.2e1 * t2 * t59 + 0.2e1 * t3 * t60 + 0.2e1 * t1 * t61 + 0.2e1 * t4 * t62 + 0.2e1 * t35 * t65 + 0.2e1 * t34 * t66 + 0.2e1 * t5 * t57 + 0.2e1 * t18 * t58 + 0.2e1 * t32 * t22 + 0.2e1 * t6 * t26 + 0.2e1 * t14 * t27 + 0.2e1 * t23 * t21 + 0.2e1 * t7 * t24 + 0.2e1 * t15 * t25 + 0.2e1 * m(4) * (t100 * t238 + t69 * t82 - t70 * t81) + t152 * t283 - 0.2e1 * (mrSges(3,1) * t219 - mrSges(3,3) * t293) * t256; t129 * t91 + t154 * t65 + t218 * t96 + (t61 + t62) * t98 + (t59 + t60) * t97 - (t24 + t25) * t231 + (t26 + t27) * t131 + (t21 + t311) * t153 + (t57 + t310) * t130 + m(5) * (t129 * t35 - t130 * t34 - t153 * t20 + t154 * t19) + m(7) * (t1 * t131 + t130 * t23 + t153 * t5 - t2 * t231 + t6 * t98 + t7 * t97) + m(6) * (t130 * t32 + t131 * t4 + t14 * t98 + t15 * t97 + t153 * t18 - t231 * t3) + (-t225 * t56 + (-t119 * t225 - t120 * t222) * mrSges(4,3) + (t225 * t103 + (-t104 + t80) * t222) * qJD(3) + m(4) * (t218 * t256 + t222 * t69 + t273 * t82 - t274 * t81 - t299) + m(5) * (t274 * t78 - t299)) * t215; (-t215 ^ 2 * t222 * t273 + t154 * t129 + t294) * t327 + 0.4e1 * (m(7) / 0.2e1 + t321) * (t131 * t98 - t231 * t97 + t294); (t303 / 0.2e1 + t55 / 0.2e1 - t20 * mrSges(5,3) + t262 * t223 + t263 * t220 + (-t220 * t260 + t223 * t261) * qJD(5) + (t301 / 0.2e1 - t63 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1 - t35 * mrSges(5,3)) * qJD(4) + (-qJD(4) * t91 + t316 + m(5) * (-qJD(4) * t35 - t20) + t311) * pkin(10)) * t221 + m(6) * (t105 * t15 + t106 * t14 + t141 * t4 + t142 * t3) - t337 * t213 / 0.2e1 + m(7) * (t1 * t124 + t133 * t2 + t140 * t23 + t179 * t5 + t6 * t95 + t7 * t99) + (-t190 / 0.2e1 + t144 / 0.2e1 + t145 / 0.2e1) * t87 + t88 * t193 / 0.2e1 + t102 * t172 / 0.2e1 + t2 * t174 + t3 * t175 + t1 * t176 + t4 * t177 + t179 * t21 + t5 * t157 + t18 * t158 + t78 * t164 + t6 * t135 + t14 * t136 + t7 * t137 + t15 * t138 + t140 * t57 + t141 * t27 + t142 * t25 + t133 * t24 + t23 * t121 + t32 * t122 + t124 * t26 + t105 * t60 + t106 * t62 + t99 * t59 + t95 * t61 - t69 * mrSges(4,2) - pkin(3) * t56 + t283 + (t298 - t320) * t70 + t246 * t53 + t247 * t52 + t248 * t90 + t249 * t89 + (t302 / 0.2e1 + t54 / 0.2e1 - t9 / 0.2e1 - t8 / 0.2e1 + t19 * mrSges(5,3) + (m(5) * t19 + t65) * pkin(10) + (t64 / 0.2e1 - t34 * mrSges(5,3) + t260 * t223 + t261 * t220 + (-m(5) * t34 + m(6) * t32 + t310) * pkin(10)) * qJD(4)) * t224 + (-t169 / 0.2e1 + t107 / 0.2e1 + t108 / 0.2e1) * t101; (t176 + t177) * t98 + (t174 + t175) * t97 + (t121 + t122) * t153 - (t137 + t138) * t231 + (t135 + t136) * t131 + (t157 + t158) * t130 + (-t225 * t164 + (-mrSges(4,2) * t225 + t222 * t298) * qJD(3)) * t215 + m(7) * (t124 * t98 + t130 * t179 + t131 * t95 + t133 * t97 + t140 * t153 - t231 * t99) + m(6) * (-t105 * t231 + t106 * t131 + t141 * t98 + t142 * t97) - t255 * t320 + (t230 * t321 + m(5) * (-t154 * t272 + t230 + t296) / 0.2e1) * t324 + (t296 + t295 + (t153 * t224 - t154 * t221) * qJD(4)) * mrSges(5,3); -0.2e1 * pkin(3) * t164 + 0.2e1 * t105 * t175 + 0.2e1 * t106 * t177 + 0.2e1 * t179 * t121 + 0.2e1 * t124 * t135 + 0.2e1 * t133 * t137 + 0.2e1 * t141 * t136 + 0.2e1 * t142 * t138 + 0.2e1 * t140 * t157 + 0.2e1 * t99 * t174 + 0.2e1 * t95 * t176 + (t105 * t142 + t106 * t141) * t326 + (t124 * t95 + t133 * t99 + t140 * t179) * t325 + (-t107 - t108 + t169 + (t158 * t324 + t220 * t282 + t223 * t281 + t193) * qJD(4)) * t224 + (t122 * t324 + t172 + (t111 + t112) * t223 + (-t109 - t110) * t220 + (-t220 * t281 + t223 * t282) * qJD(5) + (pkin(10) ^ 2 * t224 * t326 + t144 + t145 - t190) * qJD(4)) * t221; t265 + t206 * t21 + t18 * t183 - t185 * t24 + t181 * t26 + t5 * t182 + t23 * t162 + t32 * t163 + t150 * t59 + t151 * t61 + t20 * mrSges(5,1) - t19 * mrSges(5,2) + m(7) * (t1 * t181 + t150 * t7 + t151 * t6 - t185 * t2 + t206 * t5) - t240 * t53 + t241 * t52 + t242 * t87 + t243 * t90 + t244 * t89 + t245 * t101 + (-t4 * mrSges(6,3) - t1 * mrSges(7,3) + (-m(6) * t4 - t27) * pkin(11) + (-t7 * mrSges(7,3) + pkin(5) * t57 - pkin(11) * t60 + t15 * t331 + t23 * t317 + t261) * qJD(5) + t262) * t220 + (t2 * mrSges(7,3) + t3 * mrSges(6,3) + (-t14 * mrSges(6,3) - t6 * mrSges(7,3) + t260) * qJD(5) + (m(6) * (-qJD(5) * t14 + t3) + t25 - qJD(5) * t62) * pkin(11) - t263) * t223 + (-t22 - t316) * pkin(4); -t129 * mrSges(5,2) + (t162 + t163) * t153 + m(7) * (t131 * t151 - t150 * t231 + t153 * t264 + t181 * t98 - t185 * t97) + (t328 + t329) * t130 + (mrSges(7,3) - t331) * (-t98 * t220 + t97 * t223 + (-t131 * t223 + t220 * t231) * qJD(5)); t213 + t206 * t121 - t185 * t137 + t150 * t174 + t151 * t176 + t179 * t162 + t181 * t135 + t140 * t182 - pkin(4) * t122 + t221 * pkin(10) * t163 + m(7) * (t124 * t151 + t133 * t150 + t140 * t206 + t181 * t95 - t185 * t99) - t245 * t224 + ((-Ifges(5,6) + t242) * t221 + (t224 * t329 + t300) * pkin(10)) * qJD(4) + (-t106 * mrSges(6,3) - t95 * mrSges(7,3) - t244 * t221 - t241 * t271 + (-m(6) * t106 - t136) * pkin(11) + (-t133 * mrSges(7,3) + pkin(5) * t157 - pkin(11) * t175 + t142 * t331 + t179 * t317 + t221 * t240 - t247) * qJD(5) + t248) * t220 + (t105 * mrSges(6,3) + t99 * mrSges(7,3) + t243 * t221 - t240 * t271 + (m(6) * t105 + t138) * pkin(11) + (-t141 * mrSges(6,3) - t124 * mrSges(7,3) - t241 * t221 + (-m(6) * t141 - t177) * pkin(11) + t246) * qJD(5) + t249) * t223; 0.2e1 * t206 * t162 + (-t150 * t185 + t151 * t181) * t325 - 0.2e1 * pkin(4) * t163 + (t151 * t322 + t170 + t171 + (0.2e1 * pkin(5) * t328 - t185 * t322 - t188 - t189) * qJD(5)) * t220 + (0.2e1 * t150 * mrSges(7,3) + t167 + t168 + (t181 * t322 + t191 + t192) * qJD(5)) * t223; mrSges(6,1) * t4 + mrSges(7,1) * t1 - mrSges(6,2) * t3 - mrSges(7,2) * t2 + (m(7) * t1 + t26) * pkin(5) + t335; (-mrSges(6,2) - mrSges(7,2)) * t97 + (mrSges(6,1) + t259) * t98; mrSges(6,1) * t106 + mrSges(7,1) * t95 - mrSges(6,2) * t105 - mrSges(7,2) * t99 - t338 * t251 + (m(7) * t95 + t135) * pkin(5) + (-t340 - t341) * t269 + t277 + t278; -mrSges(7,2) * t150 + t211 + t212 + t259 * t151 + ((-mrSges(6,1) * pkin(11) - mrSges(7,3) * pkin(5)) * t223 + (mrSges(6,2) * pkin(11) - t338) * t220) * qJD(5); 0; m(7) * t5 + t21; m(7) * t130; m(7) * t140 + t121; m(7) * t264 + t162; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
