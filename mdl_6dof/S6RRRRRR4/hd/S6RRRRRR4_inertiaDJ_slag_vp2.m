% Calculate time derivative of joint inertia matrix for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:46:48
% EndTime: 2019-03-10 03:47:03
% DurationCPUTime: 6.86s
% Computational Cost: add. (19976->642), mult. (46646->953), div. (0->0), fcn. (43901->10), ass. (0->255)
t240 = sin(qJ(3));
t241 = sin(qJ(2));
t245 = cos(qJ(3));
t279 = qJD(3) * t245;
t246 = cos(qJ(2));
t282 = qJD(2) * t246;
t251 = t240 * t282 + t241 * t279;
t239 = sin(qJ(4));
t244 = cos(qJ(4));
t257 = t239 * t240 - t244 * t245;
t184 = t257 * t241;
t303 = -pkin(9) - pkin(8);
t217 = t303 * t240;
t218 = t303 * t245;
t171 = t244 * t217 + t218 * t239;
t201 = t239 * t245 + t240 * t244;
t145 = -pkin(10) * t201 + t171;
t172 = t239 * t217 - t244 * t218;
t146 = -pkin(10) * t257 + t172;
t238 = sin(qJ(5));
t243 = cos(qJ(5));
t99 = t238 * t145 + t243 * t146;
t214 = -pkin(2) * t246 - t241 * pkin(8) - pkin(1);
t199 = t245 * t214;
t287 = t241 * t245;
t301 = pkin(7) * t240;
t160 = -pkin(9) * t287 + t199 + (-pkin(3) - t301) * t246;
t285 = t245 * t246;
t225 = pkin(7) * t285;
t179 = t240 * t214 + t225;
t288 = t240 * t241;
t170 = -pkin(9) * t288 + t179;
t114 = t239 * t160 + t244 * t170;
t270 = t245 * t282;
t283 = qJD(2) * t241;
t313 = -Ifges(4,5) * t270 - Ifges(4,3) * t283;
t312 = qJD(3) + qJD(4);
t166 = t312 * t201;
t122 = -t166 * t241 - t257 * t282;
t123 = t184 * t312 - t201 * t282;
t183 = t201 * t241;
t135 = -t183 * t243 + t184 * t238;
t63 = qJD(5) * t135 + t122 * t243 + t123 * t238;
t136 = -t183 * t238 - t184 * t243;
t64 = -qJD(5) * t136 - t122 * t238 + t123 * t243;
t311 = -Ifges(6,5) * t63 - Ifges(6,6) * t64 - Ifges(6,3) * t283;
t310 = -Ifges(5,5) * t122 - Ifges(5,6) * t123 - Ifges(5,3) * t283;
t309 = 2 * m(4);
t308 = 2 * m(5);
t307 = 2 * m(6);
t306 = 2 * m(7);
t305 = -0.2e1 * pkin(1);
t304 = 0.2e1 * pkin(7);
t302 = -t240 / 0.2e1;
t227 = pkin(3) * t244 + pkin(4);
t291 = t238 * t239;
t187 = -pkin(3) * t291 + t243 * t227;
t185 = pkin(5) + t187;
t289 = t239 * t243;
t189 = pkin(3) * t289 + t227 * t238;
t237 = sin(qJ(6));
t242 = cos(qJ(6));
t140 = t185 * t242 - t189 * t237;
t275 = qJD(5) * t243;
t276 = qJD(5) * t238;
t153 = t227 * t275 + (-t239 * t276 + (t243 * t244 - t291) * qJD(4)) * pkin(3);
t154 = -t227 * t276 + (-t239 * t275 + (-t238 * t244 - t289) * qJD(4)) * pkin(3);
t73 = qJD(6) * t140 + t153 * t242 + t154 * t237;
t300 = t73 * mrSges(7,2);
t100 = -pkin(10) * t183 + t114;
t113 = t244 * t160 - t239 * t170;
t95 = -pkin(4) * t246 + t184 * pkin(10) + t113;
t58 = t243 * t100 + t238 * t95;
t298 = Ifges(4,4) * t240;
t297 = Ifges(4,4) * t245;
t296 = Ifges(4,6) * t240;
t226 = pkin(4) * t243 + pkin(5);
t273 = qJD(6) * t242;
t274 = qJD(6) * t237;
t292 = t237 * t238;
t151 = t226 * t273 + (-t238 * t274 + (t242 * t243 - t292) * qJD(5)) * pkin(4);
t295 = t151 * mrSges(7,2);
t294 = t153 * mrSges(6,2);
t293 = t246 * Ifges(4,6);
t290 = t238 * t242;
t216 = Ifges(4,1) * t240 + t297;
t286 = t245 * t216;
t212 = (pkin(2) * t241 - pkin(8) * t246) * qJD(2);
t284 = t245 * t212 + t283 * t301;
t213 = pkin(3) * t288 + t241 * pkin(7);
t281 = qJD(3) * t240;
t280 = qJD(3) * t241;
t278 = qJD(4) * t239;
t277 = qJD(4) * t244;
t91 = t135 * t242 - t136 * t237;
t24 = qJD(6) * t91 + t237 * t64 + t242 * t63;
t92 = t135 * t237 + t136 * t242;
t25 = -qJD(6) * t92 - t237 * t63 + t242 * t64;
t272 = -Ifges(7,5) * t24 - Ifges(7,6) * t25 - Ifges(7,3) * t283;
t235 = pkin(7) * t282;
t234 = pkin(3) * t281;
t177 = t251 * pkin(3) + t235;
t228 = -pkin(3) * t245 - pkin(2);
t271 = qJD(3) * t303;
t269 = t240 * t280;
t266 = (2 * Ifges(3,4)) + t296;
t147 = pkin(4) * t166 + t234;
t152 = -t226 * t274 + (-t238 * t273 + (-t237 * t243 - t290) * qJD(5)) * pkin(4);
t148 = t152 * mrSges(7,1);
t265 = t148 - t295;
t57 = -t238 * t100 + t243 * t95;
t98 = t243 * t145 - t146 * t238;
t149 = t154 * mrSges(6,1);
t141 = t185 * t237 + t189 * t242;
t74 = -qJD(6) * t141 - t153 * t237 + t154 * t242;
t72 = t74 * mrSges(7,1);
t264 = t149 + t72 - t294;
t163 = pkin(4) * t183 + t213;
t158 = -t201 * t238 - t243 * t257;
t159 = t201 * t243 - t238 * t257;
t104 = t158 * t237 + t159 * t242;
t165 = t312 * t257;
t79 = qJD(5) * t158 - t165 * t243 - t166 * t238;
t80 = -qJD(5) * t159 + t165 * t238 - t166 * t243;
t34 = -qJD(6) * t104 - t237 * t79 + t242 * t80;
t31 = Ifges(7,6) * t34;
t103 = t158 * t242 - t159 * t237;
t33 = qJD(6) * t103 + t237 * t80 + t242 * t79;
t32 = Ifges(7,5) * t33;
t210 = t240 * t271;
t211 = t245 * t271;
t125 = t244 * t210 + t239 * t211 + t217 * t277 + t218 * t278;
t93 = -pkin(10) * t166 + t125;
t126 = -qJD(4) * t172 - t210 * t239 + t244 * t211;
t94 = pkin(10) * t165 + t126;
t38 = t145 * t275 - t146 * t276 + t238 * t94 + t243 * t93;
t26 = pkin(11) * t80 + t38;
t39 = -qJD(5) * t99 - t238 * t93 + t243 * t94;
t27 = -pkin(11) * t79 + t39;
t75 = -pkin(11) * t159 + t98;
t76 = pkin(11) * t158 + t99;
t44 = -t237 * t76 + t242 * t75;
t5 = qJD(6) * t44 + t237 * t27 + t242 * t26;
t45 = t237 * t75 + t242 * t76;
t6 = -qJD(6) * t45 - t237 * t26 + t242 * t27;
t263 = t6 * mrSges(7,1) - t5 * mrSges(7,2) + t31 + t32;
t262 = -mrSges(4,1) * t245 + mrSges(4,2) * t240;
t261 = mrSges(4,1) * t240 + mrSges(4,2) * t245;
t260 = Ifges(4,1) * t245 - t298;
t259 = -Ifges(4,2) * t240 + t297;
t215 = Ifges(4,2) * t245 + t298;
t258 = Ifges(4,5) * t240 + Ifges(4,6) * t245;
t46 = -pkin(5) * t246 - t136 * pkin(11) + t57;
t47 = pkin(11) * t135 + t58;
t18 = -t237 * t47 + t242 * t46;
t19 = t237 * t46 + t242 * t47;
t180 = pkin(4) * t257 + t228;
t101 = -pkin(4) * t123 + t177;
t110 = (pkin(3) * t241 - pkin(9) * t285) * qJD(2) + (-t225 + (pkin(9) * t241 - t214) * t240) * qJD(3) + t284;
t130 = t240 * t212 + t214 * t279 + (-t245 * t283 - t246 * t281) * pkin(7);
t119 = -pkin(9) * t251 + t130;
t50 = -qJD(4) * t114 + t244 * t110 - t119 * t239;
t40 = pkin(4) * t283 - pkin(10) * t122 + t50;
t49 = t239 * t110 + t244 * t119 + t160 * t277 - t170 * t278;
t42 = pkin(10) * t123 + t49;
t14 = -qJD(5) * t58 - t238 * t42 + t243 * t40;
t7 = pkin(5) * t283 - pkin(11) * t63 + t14;
t13 = -t100 * t276 + t238 * t40 + t243 * t42 + t95 * t275;
t8 = pkin(11) * t64 + t13;
t2 = qJD(6) * t18 + t237 * t7 + t242 * t8;
t3 = -qJD(6) * t19 - t237 * t8 + t242 * t7;
t256 = t3 * mrSges(7,1) - t2 * mrSges(7,2) - t272;
t255 = -mrSges(7,1) * t274 - mrSges(7,2) * t273;
t254 = (-mrSges(5,1) * t239 - mrSges(5,2) * t244) * qJD(4) * pkin(3);
t253 = (-mrSges(6,1) * t238 - mrSges(6,2) * t243) * qJD(5) * pkin(4);
t252 = -t269 + t270;
t77 = Ifges(6,6) * t80;
t78 = Ifges(6,5) * t79;
t250 = t39 * mrSges(6,1) - t38 * mrSges(6,2) + t263 + t77 + t78;
t249 = t14 * mrSges(6,1) - t13 * mrSges(6,2) + t256 - t311;
t161 = Ifges(5,6) * t166;
t162 = Ifges(5,5) * t165;
t248 = t126 * mrSges(5,1) - t125 * mrSges(5,2) - t161 - t162 + t250;
t247 = t50 * mrSges(5,1) - t49 * mrSges(5,2) + t249 - t310;
t233 = Ifges(4,5) * t279;
t209 = -mrSges(4,1) * t246 - mrSges(4,3) * t287;
t208 = mrSges(4,2) * t246 - mrSges(4,3) * t288;
t207 = t260 * qJD(3);
t206 = t259 * qJD(3);
t205 = t261 * qJD(3);
t196 = (-t237 * mrSges(7,1) - t242 * mrSges(7,2)) * qJD(6) * pkin(5);
t188 = pkin(4) * t290 + t226 * t237;
t186 = -pkin(4) * t292 + t226 * t242;
t182 = -Ifges(4,5) * t246 + t241 * t260;
t181 = t241 * t259 - t293;
t178 = -t246 * t301 + t199;
t176 = -mrSges(4,2) * t283 - mrSges(4,3) * t251;
t175 = mrSges(4,1) * t283 - mrSges(4,3) * t252;
t174 = -mrSges(5,1) * t246 + t184 * mrSges(5,3);
t173 = mrSges(5,2) * t246 - t183 * mrSges(5,3);
t169 = Ifges(5,1) * t201 - Ifges(5,4) * t257;
t168 = Ifges(5,4) * t201 - Ifges(5,2) * t257;
t167 = mrSges(5,1) * t257 + mrSges(5,2) * t201;
t157 = mrSges(4,1) * t251 + mrSges(4,2) * t252;
t144 = mrSges(5,1) * t183 - mrSges(5,2) * t184;
t143 = -t216 * t280 + (Ifges(4,5) * t241 + t246 * t260) * qJD(2);
t142 = -t215 * t280 + (Ifges(4,6) * t241 + t246 * t259) * qJD(2);
t134 = -Ifges(5,1) * t184 - Ifges(5,4) * t183 - Ifges(5,5) * t246;
t133 = -Ifges(5,4) * t184 - Ifges(5,2) * t183 - Ifges(5,6) * t246;
t131 = -qJD(3) * t179 + t284;
t129 = -mrSges(6,1) * t246 - t136 * mrSges(6,3);
t128 = mrSges(6,2) * t246 + t135 * mrSges(6,3);
t127 = -pkin(5) * t158 + t180;
t118 = -Ifges(5,1) * t165 - Ifges(5,4) * t166;
t117 = -Ifges(5,4) * t165 - Ifges(5,2) * t166;
t116 = mrSges(5,1) * t166 - mrSges(5,2) * t165;
t112 = -mrSges(5,2) * t283 + mrSges(5,3) * t123;
t111 = mrSges(5,1) * t283 - mrSges(5,3) * t122;
t109 = Ifges(6,1) * t159 + Ifges(6,4) * t158;
t108 = Ifges(6,4) * t159 + Ifges(6,2) * t158;
t107 = -mrSges(6,1) * t158 + mrSges(6,2) * t159;
t102 = -pkin(5) * t135 + t163;
t96 = -mrSges(6,1) * t135 + mrSges(6,2) * t136;
t85 = Ifges(6,1) * t136 + Ifges(6,4) * t135 - Ifges(6,5) * t246;
t84 = Ifges(6,4) * t136 + Ifges(6,2) * t135 - Ifges(6,6) * t246;
t82 = -mrSges(7,1) * t246 - t92 * mrSges(7,3);
t81 = mrSges(7,2) * t246 + t91 * mrSges(7,3);
t71 = -mrSges(5,1) * t123 + mrSges(5,2) * t122;
t70 = Ifges(5,1) * t122 + Ifges(5,4) * t123 + Ifges(5,5) * t283;
t69 = Ifges(5,4) * t122 + Ifges(5,2) * t123 + Ifges(5,6) * t283;
t68 = -pkin(5) * t80 + t147;
t67 = Ifges(7,1) * t104 + Ifges(7,4) * t103;
t66 = Ifges(7,4) * t104 + Ifges(7,2) * t103;
t65 = -mrSges(7,1) * t103 + mrSges(7,2) * t104;
t60 = -mrSges(6,2) * t283 + mrSges(6,3) * t64;
t59 = mrSges(6,1) * t283 - mrSges(6,3) * t63;
t56 = -mrSges(7,1) * t91 + mrSges(7,2) * t92;
t55 = Ifges(7,1) * t92 + Ifges(7,4) * t91 - Ifges(7,5) * t246;
t54 = Ifges(7,4) * t92 + Ifges(7,2) * t91 - Ifges(7,6) * t246;
t53 = Ifges(6,1) * t79 + Ifges(6,4) * t80;
t52 = Ifges(6,4) * t79 + Ifges(6,2) * t80;
t51 = -mrSges(6,1) * t80 + mrSges(6,2) * t79;
t43 = -pkin(5) * t64 + t101;
t30 = -mrSges(6,1) * t64 + mrSges(6,2) * t63;
t29 = Ifges(6,1) * t63 + Ifges(6,4) * t64 + Ifges(6,5) * t283;
t28 = Ifges(6,4) * t63 + Ifges(6,2) * t64 + Ifges(6,6) * t283;
t21 = -mrSges(7,2) * t283 + mrSges(7,3) * t25;
t20 = mrSges(7,1) * t283 - mrSges(7,3) * t24;
t17 = Ifges(7,1) * t33 + Ifges(7,4) * t34;
t16 = Ifges(7,4) * t33 + Ifges(7,2) * t34;
t15 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t11 = -mrSges(7,1) * t25 + mrSges(7,2) * t24;
t10 = Ifges(7,1) * t24 + Ifges(7,4) * t25 + Ifges(7,5) * t283;
t9 = Ifges(7,4) * t24 + Ifges(7,2) * t25 + Ifges(7,6) * t283;
t1 = [((mrSges(3,2) * t305 - t240 * t181 + t245 * t182 + t246 * t266) * qJD(2) + t272 + t310 + t311 + t313) * t246 + (t157 * t304 - t240 * t142 + t245 * t143 + (-t245 * t181 - t240 * t182 + t246 * t258) * qJD(3) + (mrSges(3,1) * t305 + Ifges(7,5) * t92 + Ifges(7,6) * t91 + Ifges(6,5) * t136 + Ifges(6,6) * t135 - Ifges(5,5) * t184 - Ifges(5,6) * t183 + (Ifges(4,5) * t245 - t266) * t241 + (pkin(7) ^ 2 * t309 + t261 * t304 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - Ifges(5,3) - Ifges(6,3) - Ifges(7,3)) * t246) * qJD(2)) * t241 + 0.2e1 * t130 * t208 + 0.2e1 * t131 * t209 + 0.2e1 * t213 * t71 - t183 * t69 - t184 * t70 + 0.2e1 * t49 * t173 + 0.2e1 * t50 * t174 + 0.2e1 * t177 * t144 + 0.2e1 * t178 * t175 + 0.2e1 * t179 * t176 + 0.2e1 * t163 * t30 + 0.2e1 * t14 * t129 + t123 * t133 + t122 * t134 + t135 * t28 + t136 * t29 + 0.2e1 * t13 * t128 + 0.2e1 * t114 * t112 + 0.2e1 * t113 * t111 + 0.2e1 * t101 * t96 + 0.2e1 * t102 * t11 + t91 * t9 + t92 * t10 + 0.2e1 * t2 * t81 + 0.2e1 * t3 * t82 + t64 * t84 + t63 * t85 + 0.2e1 * t57 * t59 + 0.2e1 * t58 * t60 + t25 * t54 + t24 * t55 + 0.2e1 * t43 * t56 + 0.2e1 * t19 * t21 + 0.2e1 * t18 * t20 + (t102 * t43 + t18 * t3 + t19 * t2) * t306 + (t101 * t163 + t13 * t58 + t14 * t57) * t307 + (t113 * t50 + t114 * t49 + t177 * t213) * t308 + (t179 * t130 + t178 * t131) * t309; (t113 * t165 - t114 * t166 - t201 * t50 - t257 * t49) * mrSges(5,3) + (t245 * t207 / 0.2e1 + t206 * t302 - Ifges(3,6) * qJD(2) + (-t245 * t215 / 0.2e1 + t216 * t302) * qJD(3) + (mrSges(3,2) * qJD(2) + t205) * pkin(7) + (Ifges(5,5) * t201 + Ifges(6,5) * t159 + Ifges(7,5) * t104 - Ifges(5,6) * t257 + Ifges(6,6) * t158 + Ifges(7,6) * t103 + t258) * qJD(2) / 0.2e1) * t241 - t257 * t69 / 0.2e1 + ((t130 * t245 - t131 * t240 + (-t178 * t245 - t179 * t240) * qJD(3)) * pkin(8) - pkin(2) * t235) * m(4) + (t13 * t158 - t14 * t159 - t57 * t79 + t58 * t80) * mrSges(6,3) + m(7) * (t102 * t68 + t127 * t43 + t18 * t6 + t19 * t5 + t2 * t45 + t3 * t44) + m(6) * (t101 * t180 + t13 * t99 + t14 * t98 + t147 * t163 + t38 * t58 + t39 * t57) + (t103 * t2 - t104 * t3 - t18 * t33 + t19 * t34) * mrSges(7,3) + (t142 / 0.2e1 + t130 * mrSges(4,3) + pkin(8) * t176) * t245 + (t143 / 0.2e1 - t131 * mrSges(4,3) - pkin(8) * t175) * t240 + t228 * t71 + t201 * t70 / 0.2e1 + t213 * t116 + t180 * t30 - t183 * t117 / 0.2e1 - t184 * t118 / 0.2e1 - t166 * t133 / 0.2e1 + t123 * t168 / 0.2e1 + t122 * t169 / 0.2e1 + t171 * t111 + t172 * t112 + t125 * t173 + t126 * t174 + t177 * t167 + t158 * t28 / 0.2e1 + t159 * t29 / 0.2e1 + t163 * t51 - t165 * t134 / 0.2e1 - pkin(2) * t157 + t147 * t96 + t39 * t129 + t135 * t52 / 0.2e1 + t136 * t53 / 0.2e1 + t127 * t11 + t38 * t128 + t103 * t9 / 0.2e1 + t104 * t10 / 0.2e1 + t101 * t107 + t64 * t108 / 0.2e1 + t63 * t109 / 0.2e1 + t98 * t59 + t99 * t60 + t102 * t15 + t91 * t16 / 0.2e1 + t92 * t17 / 0.2e1 + t5 * t81 + t6 * t82 + t80 * t84 / 0.2e1 + t79 * t85 / 0.2e1 + t43 * t65 + t25 * t66 / 0.2e1 + t24 * t67 / 0.2e1 + t68 * t56 + t34 * t54 / 0.2e1 + t33 * t55 / 0.2e1 + t44 * t20 + t45 * t21 + m(5) * (t113 * t126 + t114 * t125 + t171 * t50 + t172 * t49 + t177 * t228 + t213 * t234) + (-t233 / 0.2e1 + t162 / 0.2e1 + t161 / 0.2e1 - t78 / 0.2e1 - t77 / 0.2e1 - t32 / 0.2e1 - t31 / 0.2e1 + (t286 / 0.2e1 + t215 * t302 + Ifges(3,5) + (-mrSges(3,1) + t262) * pkin(7)) * qJD(2)) * t246 + ((-t178 * mrSges(4,3) - pkin(8) * t209 + t182 / 0.2e1) * t245 + (t293 / 0.2e1 - t181 / 0.2e1 + pkin(3) * t144 - t179 * mrSges(4,3) - pkin(8) * t208) * t240) * qJD(3); 0.2e1 * (-t125 * t257 - t126 * t201 + t165 * t171 - t166 * t172) * mrSges(5,3) - t257 * t117 + 0.2e1 * (t103 * t5 - t104 * t6 - t33 * t44 + t34 * t45) * mrSges(7,3) + 0.2e1 * (t158 * t38 - t159 * t39 - t79 * t98 + t80 * t99) * mrSges(6,3) + t245 * t206 + t240 * t207 + 0.2e1 * t228 * t116 + t201 * t118 - 0.2e1 * pkin(2) * t205 + 0.2e1 * t180 * t51 - t166 * t168 - t165 * t169 + t158 * t52 + t159 * t53 + 0.2e1 * t147 * t107 + 0.2e1 * t127 * t15 + t103 * t16 + t104 * t17 + t80 * t108 + t79 * t109 + t34 * t66 + t33 * t67 + 0.2e1 * t68 * t65 + (t127 * t68 + t44 * t6 + t45 * t5) * t306 + (t147 * t180 + t38 * t99 + t39 * t98) * t307 + (t125 * t172 + t126 * t171 + t228 * t234) * t308 + (t286 + (0.2e1 * pkin(3) * t167 - t215) * t240) * qJD(3); (t173 * t277 + m(5) * (-t113 * t278 + t114 * t277 + t239 * t49 + t244 * t50) + t244 * t111 - t174 * t278 + t239 * t112) * pkin(3) - t251 * Ifges(4,6) + t247 + t189 * t60 + t187 * t59 + t153 * t128 + t154 * t129 + t140 * t20 + t141 * t21 - t130 * mrSges(4,2) + t131 * mrSges(4,1) + t73 * t81 + t74 * t82 + m(7) * (t140 * t3 + t141 * t2 + t18 * t74 + t19 * t73) + m(6) * (t13 * t189 + t14 * t187 + t153 * t58 + t154 * t57) - Ifges(4,5) * t269 - t313; (m(5) * (t125 * t239 + t126 * t244 + (-t171 * t239 + t172 * t244) * qJD(4)) + (t244 * t165 - t239 * t166 + (t201 * t239 - t244 * t257) * qJD(4)) * mrSges(5,3)) * pkin(3) + t233 + (pkin(8) * t262 - t296) * qJD(3) + t248 + m(7) * (t140 * t6 + t141 * t5 + t44 * t74 + t45 * t73) + m(6) * (t153 * t99 + t154 * t98 + t187 * t39 + t189 * t38) + (t103 * t73 - t104 * t74 - t140 * t33 + t141 * t34) * mrSges(7,3) + (t153 * t158 - t154 * t159 - t187 * t79 + t189 * t80) * mrSges(6,3); -0.2e1 * t294 - 0.2e1 * t300 + 0.2e1 * t149 + 0.2e1 * t72 + 0.2e1 * t254 + (t140 * t74 + t141 * t73) * t306 + (t153 * t189 + t154 * t187) * t307; (t128 * t275 + t243 * t59 + m(6) * (t13 * t238 + t14 * t243 + t275 * t58 - t276 * t57) - t129 * t276 + t238 * t60) * pkin(4) + m(7) * (t151 * t19 + t152 * t18 + t186 * t3 + t188 * t2) + t247 + t186 * t20 + t188 * t21 + t151 * t81 + t152 * t82; (m(6) * (t238 * t38 + t243 * t39 + (-t238 * t98 + t243 * t99) * qJD(5)) + (t238 * t80 - t243 * t79 + (t158 * t243 + t159 * t238) * qJD(5)) * mrSges(6,3)) * pkin(4) + m(7) * (t151 * t45 + t152 * t44 + t186 * t6 + t188 * t5) + t248 + (t103 * t151 - t104 * t152 - t186 * t33 + t188 * t34) * mrSges(7,3); m(7) * (t140 * t152 + t141 * t151 + t186 * t74 + t188 * t73) + t148 + t254 + (-t73 - t151) * mrSges(7,2) + (m(6) * (t153 * t238 + t154 * t243 - t187 * t276 + t189 * t275) - mrSges(6,2) * t275 - mrSges(6,1) * t276) * pkin(4) + t264; -0.2e1 * t295 + 0.2e1 * t148 + (t151 * t188 + t152 * t186) * t306 + 0.2e1 * t253; (-t82 * t274 + t242 * t20 + m(7) * (-t18 * t274 + t19 * t273 + t2 * t237 + t242 * t3) + t81 * t273 + t237 * t21) * pkin(5) + t249; (m(7) * (t237 * t5 + t242 * t6 + (-t237 * t44 + t242 * t45) * qJD(6)) + (t237 * t34 - t242 * t33 + (t103 * t242 + t104 * t237) * qJD(6)) * mrSges(7,3)) * pkin(5) + t250; -t300 + (m(7) * (-t140 * t274 + t141 * t273 + t237 * t73 + t242 * t74) + t255) * pkin(5) + t264; t253 + (m(7) * (t151 * t237 + t152 * t242 - t186 * t274 + t188 * t273) + t255) * pkin(5) + t265; 0.2e1 * t196; t256; t263; t72 - t300; t265; t196; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
