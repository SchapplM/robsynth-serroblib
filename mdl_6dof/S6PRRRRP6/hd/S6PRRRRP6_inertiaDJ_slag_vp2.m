% Calculate time derivative of joint inertia matrix for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:42
% EndTime: 2019-03-09 00:28:56
% DurationCPUTime: 6.23s
% Computational Cost: add. (5813->649), mult. (16902->908), div. (0->0), fcn. (16047->12), ass. (0->260)
t207 = sin(qJ(5));
t331 = (Ifges(7,4) + Ifges(6,5)) * t207;
t208 = sin(qJ(4));
t212 = cos(qJ(4));
t270 = qJD(4) * t212;
t248 = t207 * t270;
t211 = cos(qJ(5));
t267 = qJD(5) * t211;
t218 = t208 * t267 + t248;
t205 = cos(pkin(7));
t203 = sin(pkin(7));
t209 = sin(qJ(3));
t289 = t203 * t209;
t140 = t205 * t208 + t212 * t289;
t213 = cos(qJ(3));
t273 = qJD(3) * t203;
t251 = t213 * t273;
t100 = qJD(4) * t140 + t208 * t251;
t288 = t203 * t213;
t101 = t140 * t207 + t211 * t288;
t272 = qJD(3) * t209;
t252 = t203 * t272;
t139 = -t212 * t205 + t208 * t289;
t99 = -qJD(4) * t139 + t212 * t251;
t50 = -qJD(5) * t101 + t207 * t252 + t211 * t99;
t255 = t207 * t288;
t51 = -qJD(5) * t255 + t140 * t267 + t207 * t99 - t211 * t252;
t12 = Ifges(6,5) * t50 - Ifges(6,6) * t51 + Ifges(6,3) * t100;
t13 = Ifges(7,4) * t50 + Ifges(7,2) * t100 + Ifges(7,6) * t51;
t328 = t12 + t13;
t191 = pkin(9) * t289;
t316 = pkin(2) * t213;
t142 = t205 * t316 - t191;
t214 = cos(qJ(2));
t282 = t213 * t214;
t210 = sin(qJ(2));
t285 = t209 * t210;
t327 = t205 * t282 - t285;
t297 = t208 * mrSges(5,2);
t326 = -m(5) * pkin(3) - t212 * mrSges(5,1) - mrSges(4,1) + t297;
t325 = m(7) * qJ(6) + mrSges(7,3);
t166 = (pkin(4) * t208 - pkin(11) * t212) * qJD(4);
t171 = -pkin(4) * t212 - pkin(11) * t208 - pkin(3);
t266 = qJD(5) * t212;
t269 = qJD(5) * t207;
t271 = qJD(4) * t208;
t76 = pkin(10) * (t207 * t271 - t211 * t266) + t166 * t211 - t171 * t269;
t123 = t191 + (-pkin(3) - t316) * t205;
t66 = pkin(4) * t139 - pkin(11) * t140 + t123;
t143 = t205 * t209 * pkin(2) + pkin(9) * t288;
t124 = pkin(10) * t205 + t143;
t125 = (-pkin(3) * t213 - pkin(10) * t209 - pkin(2)) * t203;
t78 = t212 * t124 + t208 * t125;
t68 = -pkin(11) * t288 + t78;
t308 = t207 * t66 + t211 * t68;
t134 = (pkin(3) * t209 - pkin(10) * t213) * t273;
t135 = t142 * qJD(3);
t35 = -t124 * t271 + t125 * t270 + t208 * t134 + t212 * t135;
t31 = pkin(11) * t252 + t35;
t136 = t143 * qJD(3);
t44 = pkin(4) * t100 - pkin(11) * t99 + t136;
t5 = -qJD(5) * t308 - t207 * t31 + t211 * t44;
t228 = pkin(5) * t211 + qJ(6) * t207;
t265 = qJD(6) * t211;
t324 = qJD(5) * t228 - t265;
t323 = 2 * m(6);
t322 = 2 * m(7);
t321 = 0.2e1 * pkin(10);
t320 = -2 * mrSges(4,3);
t319 = m(6) / 0.2e1;
t36 = -t124 * t270 - t125 * t271 + t134 * t212 - t208 * t135;
t32 = -pkin(4) * t252 - t36;
t317 = m(6) * t32;
t315 = pkin(10) * t212;
t204 = sin(pkin(6));
t206 = cos(pkin(6));
t283 = t210 * t213;
t284 = t209 * t214;
t221 = t205 * t284 + t283;
t59 = t206 * t252 + (t221 * qJD(3) + (t205 * t283 + t284) * qJD(2)) * t204;
t93 = -t204 * t327 - t206 * t288;
t314 = t59 * t93;
t274 = qJD(2) * t204;
t253 = t210 * t274;
t236 = t203 * t253;
t60 = t206 * t251 + (t327 * qJD(3) + (-t205 * t285 + t282) * qJD(2)) * t204;
t138 = -t203 * t204 * t214 + t206 * t205;
t94 = t204 * t221 + t206 * t289;
t65 = t138 * t208 + t212 * t94;
t21 = qJD(4) * t65 + t208 * t60 - t212 * t236;
t224 = t138 * t212 - t208 * t94;
t10 = t224 * t21;
t22 = qJD(4) * t224 + t208 * t236 + t212 * t60;
t34 = t207 * t93 + t211 * t65;
t7 = qJD(5) * t34 + t207 * t22 - t59 * t211;
t313 = t7 * t207;
t33 = t207 * t65 - t93 * t211;
t8 = -qJD(5) * t33 + t207 * t59 + t211 * t22;
t312 = t8 * t211;
t18 = mrSges(6,1) * t51 + mrSges(6,2) * t50;
t81 = mrSges(5,1) * t252 - mrSges(5,3) * t99;
t311 = t18 - t81;
t26 = mrSges(6,1) * t100 - mrSges(6,3) * t50;
t27 = -t100 * mrSges(7,1) + t50 * mrSges(7,2);
t310 = -t26 + t27;
t28 = -mrSges(6,2) * t100 - mrSges(6,3) * t51;
t29 = -mrSges(7,2) * t51 + mrSges(7,3) * t100;
t309 = t28 + t29;
t70 = -mrSges(7,2) * t101 + mrSges(7,3) * t139;
t71 = -mrSges(6,2) * t139 - mrSges(6,3) * t101;
t307 = t70 + t71;
t102 = t140 * t211 - t255;
t72 = mrSges(6,1) * t139 - mrSges(6,3) * t102;
t73 = -mrSges(7,1) * t139 + mrSges(7,2) * t102;
t306 = t73 - t72;
t305 = Ifges(5,4) * t208;
t304 = Ifges(5,4) * t212;
t303 = Ifges(6,4) * t207;
t302 = Ifges(6,4) * t211;
t301 = Ifges(7,5) * t207;
t300 = Ifges(7,5) * t211;
t299 = Ifges(6,6) * t211;
t298 = t136 * t93;
t296 = t21 * t208;
t295 = t22 * t212;
t173 = -t211 * mrSges(6,1) + t207 * mrSges(6,2);
t293 = t173 - mrSges(5,1);
t104 = -mrSges(5,1) * t288 - mrSges(5,3) * t140;
t56 = mrSges(6,1) * t101 + mrSges(6,2) * t102;
t292 = t56 - t104;
t291 = -mrSges(4,1) * t205 + mrSges(5,1) * t139 + mrSges(5,2) * t140 + mrSges(4,3) * t289;
t290 = t171 * t211;
t287 = t207 * t208;
t286 = t208 * t211;
t247 = t211 * t270;
t268 = qJD(5) * t208;
t250 = t207 * t268;
t219 = t247 - t250;
t110 = mrSges(6,1) * t271 - mrSges(6,3) * t219;
t111 = mrSges(7,2) * t247 + (-mrSges(7,1) * qJD(4) - mrSges(7,2) * t269) * t208;
t281 = -t110 + t111;
t112 = -mrSges(6,2) * t271 - mrSges(6,3) * t218;
t113 = -mrSges(7,2) * t218 + mrSges(7,3) * t271;
t280 = t112 + t113;
t229 = Ifges(7,3) * t207 + t300;
t126 = -Ifges(7,6) * t212 + t208 * t229;
t230 = -Ifges(6,2) * t207 + t302;
t129 = -Ifges(6,6) * t212 + t208 * t230;
t279 = t126 - t129;
t231 = Ifges(7,1) * t211 + t301;
t130 = -Ifges(7,4) * t212 + t208 * t231;
t232 = Ifges(6,1) * t211 - t303;
t131 = -Ifges(6,5) * t212 + t208 * t232;
t278 = t130 + t131;
t162 = mrSges(6,2) * t212 - mrSges(6,3) * t287;
t165 = -mrSges(7,2) * t287 - mrSges(7,3) * t212;
t277 = t162 + t165;
t163 = -mrSges(6,1) * t212 - mrSges(6,3) * t286;
t164 = mrSges(7,1) * t212 + mrSges(7,2) * t286;
t276 = -t163 + t164;
t275 = Ifges(6,5) * t247 + Ifges(6,3) * t271;
t121 = t207 * t171 + t211 * t315;
t156 = Ifges(7,4) * t267 + Ifges(7,6) * t269;
t263 = Ifges(5,5) * t99 - Ifges(5,6) * t100 + Ifges(5,3) * t252;
t262 = Ifges(5,6) * t288;
t11 = Ifges(7,5) * t50 + Ifges(7,6) * t100 + Ifges(7,3) * t51;
t14 = Ifges(6,4) * t50 - Ifges(6,2) * t51 + Ifges(6,6) * t100;
t261 = t11 / 0.2e1 - t14 / 0.2e1;
t15 = Ifges(7,1) * t50 + Ifges(7,4) * t100 + Ifges(7,5) * t51;
t16 = Ifges(6,1) * t50 - Ifges(6,4) * t51 + Ifges(6,5) * t100;
t260 = t15 / 0.2e1 + t16 / 0.2e1;
t38 = Ifges(7,5) * t102 + Ifges(7,6) * t139 + Ifges(7,3) * t101;
t41 = Ifges(6,4) * t102 - Ifges(6,2) * t101 + Ifges(6,6) * t139;
t259 = t38 / 0.2e1 - t41 / 0.2e1;
t42 = Ifges(7,1) * t102 + Ifges(7,4) * t139 + Ifges(7,5) * t101;
t43 = Ifges(6,1) * t102 - Ifges(6,4) * t101 + Ifges(6,5) * t139;
t258 = t42 / 0.2e1 + t43 / 0.2e1;
t175 = -Ifges(7,3) * t211 + t301;
t83 = -t175 * t268 + (Ifges(7,6) * t208 + t212 * t229) * qJD(4);
t178 = Ifges(6,2) * t211 + t303;
t86 = -t178 * t268 + (Ifges(6,6) * t208 + t212 * t230) * qJD(4);
t257 = t83 / 0.2e1 - t86 / 0.2e1;
t180 = Ifges(7,1) * t207 - t300;
t87 = -t180 * t268 + (Ifges(7,4) * t208 + t212 * t231) * qJD(4);
t181 = Ifges(6,1) * t207 + t302;
t88 = -t181 * t268 + (Ifges(6,5) * t208 + t212 * t232) * qJD(4);
t256 = t87 / 0.2e1 + t88 / 0.2e1;
t246 = t126 / 0.2e1 - t129 / 0.2e1;
t245 = t130 / 0.2e1 + t131 / 0.2e1;
t154 = t229 * qJD(5);
t157 = t230 * qJD(5);
t244 = t154 / 0.2e1 - t157 / 0.2e1;
t155 = Ifges(6,5) * t267 - Ifges(6,6) * t269;
t243 = t155 / 0.2e1 + t156 / 0.2e1;
t159 = t231 * qJD(5);
t160 = t232 * qJD(5);
t242 = t159 / 0.2e1 + t160 / 0.2e1;
t241 = t175 / 0.2e1 - t178 / 0.2e1;
t240 = t299 / 0.2e1 - Ifges(7,6) * t211 / 0.2e1 + t331 / 0.2e1;
t239 = t180 / 0.2e1 + t181 / 0.2e1;
t77 = -t208 * t124 + t125 * t212;
t237 = Ifges(7,4) * t247 + Ifges(7,2) * t271 + t218 * Ifges(7,6);
t235 = t252 / 0.2e1;
t67 = pkin(4) * t288 - t77;
t234 = mrSges(6,1) * t207 + mrSges(6,2) * t211;
t172 = -t211 * mrSges(7,1) - t207 * mrSges(7,3);
t233 = mrSges(7,1) * t207 - mrSges(7,3) * t211;
t227 = pkin(5) * t207 - qJ(6) * t211;
t23 = -t207 * t68 + t211 * t66;
t223 = pkin(10) + t227;
t222 = -t224 * t270 + t296;
t4 = t207 * t44 + t211 * t31 + t66 * t267 - t269 * t68;
t75 = t207 * t166 + t171 * t267 + (-t207 * t266 - t211 * t271) * pkin(10);
t201 = Ifges(5,5) * t270;
t185 = Ifges(4,5) * t251;
t182 = Ifges(5,1) * t208 + t304;
t179 = Ifges(5,2) * t212 + t305;
t168 = -pkin(4) - t228;
t161 = (Ifges(5,1) * t212 - t305) * qJD(4);
t158 = (-Ifges(5,2) * t208 + t304) * qJD(4);
t153 = (mrSges(5,1) * t208 + mrSges(5,2) * t212) * qJD(4);
t152 = t234 * qJD(5);
t151 = t233 * qJD(5);
t150 = -mrSges(4,2) * t205 + mrSges(4,3) * t288;
t145 = t234 * t208;
t144 = t233 * t208;
t137 = qJD(5) * t227 - qJD(6) * t207;
t133 = (mrSges(4,1) * t209 + mrSges(4,2) * t213) * t273;
t132 = t223 * t208;
t128 = -Ifges(7,2) * t212 + (Ifges(7,4) * t211 + Ifges(7,6) * t207) * t208;
t127 = -Ifges(6,3) * t212 + (Ifges(6,5) * t211 - Ifges(6,6) * t207) * t208;
t120 = -t207 * t315 + t290;
t106 = -t290 + (pkin(10) * t207 + pkin(5)) * t212;
t105 = -qJ(6) * t212 + t121;
t103 = mrSges(5,2) * t288 - mrSges(5,3) * t139;
t91 = mrSges(6,1) * t218 + mrSges(6,2) * t219;
t90 = mrSges(7,1) * t218 - mrSges(7,3) * t219;
t85 = -Ifges(7,4) * t250 + t237;
t84 = -Ifges(6,5) * t250 - Ifges(6,6) * t218 + t275;
t82 = -mrSges(5,2) * t252 - mrSges(5,3) * t100;
t80 = Ifges(5,1) * t140 - Ifges(5,4) * t139 - Ifges(5,5) * t288;
t79 = Ifges(5,4) * t140 - Ifges(5,2) * t139 - t262;
t74 = t208 * t324 + t223 * t270;
t69 = -pkin(5) * t271 - t76;
t62 = qJ(6) * t271 - qJD(6) * t212 + t75;
t55 = mrSges(7,1) * t101 - mrSges(7,3) * t102;
t54 = mrSges(5,1) * t100 + mrSges(5,2) * t99;
t53 = Ifges(5,1) * t99 - Ifges(5,4) * t100 + Ifges(5,5) * t252;
t52 = Ifges(5,4) * t99 - Ifges(5,2) * t100 + Ifges(5,6) * t252;
t40 = Ifges(7,4) * t102 + Ifges(7,2) * t139 + Ifges(7,6) * t101;
t39 = Ifges(6,5) * t102 - Ifges(6,6) * t101 + Ifges(6,3) * t139;
t25 = pkin(5) * t101 - qJ(6) * t102 + t67;
t20 = -pkin(5) * t139 - t23;
t19 = qJ(6) * t139 + t308;
t17 = mrSges(7,1) * t51 - mrSges(7,3) * t50;
t9 = pkin(5) * t51 - qJ(6) * t50 - qJD(6) * t102 + t32;
t6 = pkin(11) * t312;
t3 = -pkin(5) * t100 - t5;
t2 = qJ(6) * t100 + qJD(6) * t139 + t4;
t1 = [0.2e1 * m(5) * (t22 * t65 - t10 + t314) + 0.2e1 * m(4) * (t138 * t236 + t60 * t94 + t314) + 0.2e1 * (m(6) + m(7)) * (t33 * t7 + t34 * t8 - t10); t22 * t103 + t138 * t133 + t60 * t150 + t93 * t54 + t65 * t82 + t307 * t8 + t306 * t7 + t291 * t59 + t309 * t34 + t310 * t33 + (-mrSges(3,1) * t210 - mrSges(3,2) * t214) * t274 - (t17 + t311) * t224 + (t55 + t292) * t21 + ((-mrSges(4,1) * t213 + mrSges(4,2) * t209) * t236 + (-t209 * t94 + t213 * t93) * qJD(3) * mrSges(4,3)) * t203 + m(4) * (-pkin(2) * t203 ^ 2 * t253 + t135 * t94 - t142 * t59 + t143 * t60 + t298) + m(5) * (t123 * t59 - t21 * t77 + t22 * t78 + t224 * t36 + t35 * t65 + t298) + m(6) * (t21 * t67 - t224 * t32 - t23 * t7 + t308 * t8 - t33 * t5 + t34 * t4) + m(7) * (t19 * t8 + t2 * t34 + t20 * t7 + t21 * t25 - t224 * t9 + t3 * t33); (t38 - t41) * t51 + (t42 + t43) * t50 + 0.2e1 * t20 * t27 + 0.2e1 * t19 * t29 + 0.2e1 * t25 * t17 + 0.2e1 * t23 * t26 + (t19 * t2 + t20 * t3 + t25 * t9) * t322 + 0.2e1 * t291 * t136 + 0.2e1 * t308 * t28 + (t23 * t5 + t308 * t4 + t32 * t67) * t323 + (t15 + t16) * t102 + (t11 - t14) * t101 + (t39 + t40 - t79) * t100 + t205 * t185 + (-t52 + t328) * t139 + 0.2e1 * t3 * t73 + 0.2e1 * t77 * t81 + 0.2e1 * t67 * t18 + 0.2e1 * t2 * t70 + 0.2e1 * t4 * t71 + 0.2e1 * t5 * t72 + t99 * t80 + 0.2e1 * t35 * t103 + 0.2e1 * t36 * t104 + 0.2e1 * t123 * t54 + 0.2e1 * t78 * t82 + t140 * t53 + 0.2e1 * t135 * t150 + 0.2e1 * m(4) * (t135 * t143 - t136 * t142) + 0.2e1 * m(5) * (t123 * t136 + t35 * t78 + t36 * t77) + (-t213 * t263 - 0.2e1 * pkin(2) * t133 + ((0.2e1 * Ifges(4,4) * t288 + Ifges(4,5) * t205 + t142 * t320) * t213 + (-0.2e1 * Ifges(4,4) * t289 + t143 * t320 + Ifges(5,5) * t140 - 0.2e1 * Ifges(4,6) * t205 - Ifges(5,6) * t139 + ((2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3)) * t288) * t209) * qJD(3)) * t203 + 0.2e1 * t9 * t55 + 0.2e1 * t32 * t56; -t60 * mrSges(4,2) + t93 * t153 + t277 * t8 + t276 * t7 - (t90 + t91) * t224 + t280 * t34 + t281 * t33 + (t144 + t145) * t21 + m(6) * (-t120 * t7 + t121 * t8 - t33 * t76 + t34 * t75) + m(7) * (t105 * t8 + t106 * t7 + t132 * t21 - t224 * t74 + t33 * t69 + t34 * t62) + (t222 * t319 + m(5) * (-t271 * t65 + t222 + t295) / 0.2e1) * t321 + (t296 + t295 + (-t208 * t65 - t212 * t224) * qJD(4)) * mrSges(5,3) + t326 * t59; (-t36 * mrSges(5,3) + Ifges(5,5) * t235 + t53 / 0.2e1 + t260 * t211 + t261 * t207 + (-t207 * t258 + t211 * t259) * qJD(5) + (-t78 * mrSges(5,3) + t262 / 0.2e1 - t79 / 0.2e1 + t39 / 0.2e1 + t40 / 0.2e1) * qJD(4) + (-qJD(4) * t103 + m(5) * (-qJD(4) * t78 - t36) + t317 + t311) * pkin(10)) * t208 + t185 + (-Ifges(4,6) * t272 - t213 * t201 / 0.2e1) * t203 + t256 * t102 + t257 * t101 + t245 * t50 + t246 * t51 + m(7) * (t105 * t2 + t106 * t3 + t132 * t9 + t19 * t62 + t20 * t69 + t25 * t74) + (t84 / 0.2e1 + t85 / 0.2e1 - t158 / 0.2e1) * t139 + (t127 / 0.2e1 + t128 / 0.2e1 - t179 / 0.2e1) * t100 + m(6) * (t120 * t5 + t121 * t4 + t23 * t76 + t308 * t75) + t308 * t112 + t326 * t136 + t74 * t55 + t75 * t71 + t76 * t72 + t62 * t70 + t69 * t73 + t105 * t29 + t106 * t27 + t23 * t110 + t20 * t111 + t19 * t113 + t120 * t26 + t121 * t28 + t132 * t17 - t135 * mrSges(4,2) + t25 * t90 + t67 * t91 + t9 * t144 + t32 * t145 + t123 * t153 + t140 * t161 / 0.2e1 + t4 * t162 + t5 * t163 + t3 * t164 + t2 * t165 + t99 * t182 / 0.2e1 + (t35 * mrSges(5,3) + Ifges(5,6) * t235 - t12 / 0.2e1 - t13 / 0.2e1 + t52 / 0.2e1 + (m(5) * t35 + t82) * pkin(10) + (t80 / 0.2e1 - t77 * mrSges(5,3) + t258 * t211 + t259 * t207 + (-m(5) * t77 + m(6) * t67 + t292) * pkin(10)) * qJD(4)) * t212 - pkin(3) * t54; -0.2e1 * pkin(3) * t153 + 0.2e1 * t105 * t113 + 0.2e1 * t106 * t111 + 0.2e1 * t120 * t110 + 0.2e1 * t121 * t112 + 0.2e1 * t132 * t90 + 0.2e1 * t74 * t144 + 0.2e1 * t75 * t162 + 0.2e1 * t76 * t163 + 0.2e1 * t69 * t164 + 0.2e1 * t62 * t165 + (t120 * t76 + t121 * t75) * t323 + (t105 * t62 + t106 * t69 + t132 * t74) * t322 + (t158 - t84 - t85 + (t145 * t321 + t207 * t279 + t211 * t278 + t182) * qJD(4)) * t212 + (t91 * t321 + t161 + (t87 + t88) * t211 + (t83 - t86) * t207 + (-t207 * t278 + t211 * t279) * qJD(5) + (pkin(10) ^ 2 * t212 * t323 + t127 + t128 - t179) * qJD(4)) * t208; -t22 * mrSges(5,2) - (t151 + t152) * t224 + (t172 + t293) * t21 + m(6) * (-pkin(4) * t21 + t6) + m(7) * (-t137 * t224 + t168 * t21 + t6) + 0.2e1 * (t319 + m(7) / 0.2e1) * (t267 * t33 - t269 * t34 + t313) * pkin(11) + (mrSges(6,3) + mrSges(7,2)) * (t313 + t312 + (-t207 * t34 + t211 * t33) * qJD(5)); -t35 * mrSges(5,2) + t36 * mrSges(5,1) + ((t20 * mrSges(7,2) - t23 * mrSges(6,3) + t258) * t211 + (-t19 * mrSges(7,2) - mrSges(6,3) * t308 + t259) * t207) * qJD(5) + (t3 * mrSges(7,2) - t5 * mrSges(6,3) + t260) * t207 + (t2 * mrSges(7,2) + t4 * mrSges(6,3) - t261) * t211 + t240 * t100 + t241 * t51 + t242 * t102 + t243 * t139 + t244 * t101 + t239 * t50 + (t309 * t211 + t310 * t207 + (-t207 * t307 + t211 * t306) * qJD(5) + m(6) * (-t207 * t5 + t211 * t4 - t23 * t267 - t269 * t308) + m(7) * (-t19 * t269 + t2 * t211 + t20 * t267 + t207 * t3)) * pkin(11) + t263 + m(7) * (t137 * t25 + t168 * t9) + t137 * t55 + t25 * t151 + t67 * t152 + t168 * t17 + t9 * t172 + t32 * t173 + (-t18 - t317) * pkin(4); t201 + t208 * pkin(10) * t152 - pkin(4) * t91 + t137 * t144 + t132 * t151 + t168 * t90 + t74 * t172 + m(7) * (t132 * t137 + t168 * t74) - t243 * t212 + ((-Ifges(5,6) + t240) * t208 + (t297 + (-m(6) * pkin(4) + t293) * t212) * pkin(10)) * qJD(4) + (t75 * mrSges(6,3) + t62 * mrSges(7,2) + t242 * t208 + t239 * t270 + (t106 * mrSges(7,2) - t120 * mrSges(6,3) + t241 * t208 + t245) * qJD(5) + (t276 * qJD(5) + m(6) * (-qJD(5) * t120 + t75) + m(7) * (qJD(5) * t106 + t62) + t280) * pkin(11) - t257) * t211 + (t69 * mrSges(7,2) - t76 * mrSges(6,3) + t244 * t208 + t241 * t270 + (-t105 * mrSges(7,2) - t121 * mrSges(6,3) - t239 * t208 + t246) * qJD(5) + (-t277 * qJD(5) + m(6) * (-qJD(5) * t121 - t76) + m(7) * (-qJD(5) * t105 + t69) + t281) * pkin(11) + t256) * t207; -0.2e1 * pkin(4) * t152 + 0.2e1 * t151 * t168 + (-t154 + t157) * t211 + (t159 + t160) * t207 + 0.2e1 * (m(7) * t168 + t172) * t137 + ((t180 + t181) * t211 + (t175 - t178) * t207) * qJD(5); m(7) * qJD(6) * t34 + (-mrSges(6,2) + t325) * t8 + (-m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1)) * t7; t2 * mrSges(7,3) + qJD(6) * t70 + qJ(6) * t29 + m(7) * (-pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t19) + t5 * mrSges(6,1) - t3 * mrSges(7,1) - t4 * mrSges(6,2) - pkin(5) * t27 + t328; qJD(6) * t165 + qJ(6) * t113 + m(7) * (-pkin(5) * t69 + qJ(6) * t62 + qJD(6) * t105) + t62 * mrSges(7,3) - Ifges(6,6) * t248 - t69 * mrSges(7,1) - t75 * mrSges(6,2) + t76 * mrSges(6,1) - pkin(5) * t111 + (-t299 - t331) * t268 + t237 + t275; -t324 * mrSges(7,2) + (m(7) * t265 + (-m(7) * t228 + t172 + t173) * qJD(5)) * pkin(11) + t155 + t156; 0.2e1 * t325 * qJD(6); m(7) * t7; m(7) * t3 + t27; m(7) * t69 + t111; (m(7) * pkin(11) + mrSges(7,2)) * t267; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
