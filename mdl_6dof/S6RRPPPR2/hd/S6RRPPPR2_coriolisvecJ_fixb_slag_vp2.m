% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2018-11-23 16:43
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPPPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:42:48
% EndTime: 2018-11-23 16:42:58
% DurationCPUTime: 9.75s
% Computational Cost: add. (6599->558), mult. (17021->744), div. (0->0), fcn. (12140->8), ass. (0->263)
t204 = sin(pkin(10));
t206 = cos(pkin(10));
t207 = sin(qJ(6));
t209 = cos(qJ(6));
t331 = -t204 * t207 + t206 * t209;
t165 = t331 * qJD(6);
t205 = sin(pkin(9));
t208 = sin(qJ(2));
t210 = cos(qJ(2));
t265 = cos(pkin(9));
t175 = t205 * t210 + t208 * t265;
t163 = t175 * qJD(1);
t336 = t331 * t163;
t333 = t165 + t336;
t340 = mrSges(5,2) - mrSges(4,1);
t221 = t209 * t204 + t207 * t206;
t166 = t221 * qJD(6);
t213 = t221 * t163;
t347 = t166 + t213;
t197 = -pkin(2) * t210 - pkin(1);
t250 = qJD(1) * t197;
t181 = qJD(3) + t250;
t235 = t265 * t210;
t231 = qJD(1) * t235;
t249 = qJD(1) * t208;
t161 = t205 * t249 - t231;
t136 = -qJD(2) * t204 + t161 * t206;
t273 = t136 * Ifges(6,6);
t135 = -t206 * qJD(2) - t204 * t161;
t274 = t135 * Ifges(6,5);
t212 = -t163 * qJ(4) + t181;
t287 = pkin(3) + qJ(5);
t70 = t161 * t287 + t212;
t285 = -qJ(3) - pkin(7);
t185 = t285 * t210;
t179 = qJD(1) * t185;
t167 = t205 * t179;
t184 = t285 * t208;
t178 = qJD(1) * t184;
t172 = qJD(2) * pkin(2) + t178;
t128 = t172 * t265 + t167;
t217 = qJD(4) - t128;
t295 = t163 * pkin(4);
t74 = -qJD(2) * t287 + t217 + t295;
t32 = -t204 * t70 + t206 * t74;
t33 = t204 * t74 + t206 * t70;
t20 = pkin(5) * t163 + pkin(8) * t135 + t32;
t21 = pkin(8) * t136 + t33;
t5 = t20 * t209 - t207 * t21;
t6 = t20 * t207 + t209 * t21;
t99 = t161 * pkin(3) + t212;
t346 = t32 * mrSges(6,1) + t5 * mrSges(7,1) + t181 * mrSges(4,2) - t33 * mrSges(6,2) - t6 * mrSges(7,2) - t99 * mrSges(5,3) - t274 / 0.2e1 + t273 / 0.2e1;
t162 = t175 * qJD(2);
t151 = qJD(1) * t162;
t246 = qJD(1) * qJD(2);
t240 = t208 * t246;
t152 = qJD(2) * t231 - t205 * t240;
t191 = pkin(2) * t240;
t215 = -qJ(4) * t152 - qJD(4) * t163 + t191;
t44 = qJD(5) * t161 + t151 * t287 + t215;
t237 = qJD(2) * t285;
t157 = qJD(3) * t210 + t208 * t237;
t145 = t157 * qJD(1);
t158 = -t208 * qJD(3) + t210 * t237;
t211 = qJD(1) * t158;
t94 = t145 * t205 - t265 * t211;
t58 = pkin(4) * t152 - qJD(2) * qJD(5) + t94;
t57 = t206 * t58;
t13 = pkin(5) * t152 + t57 + (-pkin(8) * t151 - t44) * t204;
t18 = t204 * t58 + t206 * t44;
t258 = t151 * t206;
t15 = pkin(8) * t258 + t18;
t1 = qJD(6) * t5 + t13 * t207 + t15 * t209;
t2 = -qJD(6) * t6 + t13 * t209 - t15 * t207;
t345 = -t1 * t221 - t2 * t331 - t333 * t6 + t347 * t5;
t233 = t135 * t207 + t209 * t136;
t41 = qJD(6) * t233 + t221 * t151;
t317 = t41 / 0.2e1;
t80 = t135 * t209 - t136 * t207;
t42 = qJD(6) * t80 + t151 * t331;
t316 = t42 / 0.2e1;
t311 = t152 / 0.2e1;
t241 = t265 * pkin(2);
t196 = -t241 - pkin(3);
t190 = -qJ(5) + t196;
t286 = -pkin(8) + t190;
t159 = t286 * t204;
t160 = t286 * t206;
t118 = -t159 * t207 + t160 * t209;
t199 = pkin(2) * t249;
t234 = qJ(4) * t161 + t199;
t76 = t163 * t287 + t234;
t236 = t265 * t179;
t130 = t178 * t205 - t236;
t300 = pkin(4) * t161;
t96 = t130 - t300;
t92 = t206 * t96;
t24 = -pkin(5) * t161 + t92 + (-pkin(8) * t163 - t76) * t204;
t256 = t163 * t206;
t37 = t204 * t96 + t206 * t76;
t28 = pkin(8) * t256 + t37;
t342 = -qJD(5) * t221 + qJD(6) * t118 - t207 * t24 - t209 * t28;
t119 = t159 * t209 + t160 * t207;
t341 = -qJD(5) * t331 - qJD(6) * t119 + t207 * t28 - t209 * t24;
t339 = -Ifges(4,1) - Ifges(6,3);
t338 = -Ifges(4,5) + Ifges(5,4);
t337 = -Ifges(4,6) + Ifges(5,5);
t283 = mrSges(4,3) * t161;
t141 = -qJD(2) * mrSges(4,2) - t283;
t270 = t161 * mrSges(5,1);
t143 = -qJD(2) * mrSges(5,3) + t270;
t335 = t141 - t143;
t282 = mrSges(4,3) * t163;
t284 = mrSges(5,1) * t163;
t334 = -qJD(2) * t340 - t282 - t284;
t259 = t151 * t204;
t108 = mrSges(6,1) * t152 - mrSges(6,3) * t259;
t109 = -mrSges(6,2) * t152 + mrSges(6,3) * t258;
t330 = t206 * t108 + t204 * t109;
t100 = -mrSges(6,2) * t163 + mrSges(6,3) * t136;
t101 = mrSges(6,1) * t163 + mrSges(6,3) * t135;
t329 = -t100 * t204 - t101 * t206;
t17 = -t204 * t44 + t57;
t226 = t17 * t206 + t18 * t204;
t248 = qJD(1) * t210;
t263 = Ifges(3,6) * qJD(2);
t280 = Ifges(3,4) * t208;
t328 = t263 / 0.2e1 + (t210 * Ifges(3,2) + t280) * qJD(1) / 0.2e1 + pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t248);
t34 = -mrSges(7,1) * t233 - mrSges(7,2) * t80;
t88 = -mrSges(6,1) * t136 - mrSges(6,2) * t135;
t245 = t143 - t34 - t88;
t327 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t41 + Ifges(7,6) * t42;
t279 = Ifges(6,4) * t204;
t228 = Ifges(6,2) * t206 + t279;
t278 = Ifges(6,4) * t206;
t281 = Ifges(6,1) * t204;
t229 = t278 + t281;
t230 = -mrSges(6,1) * t206 + mrSges(6,2) * t204;
t129 = t205 * t172 - t236;
t124 = -qJD(2) * qJ(4) - t129;
t86 = qJD(5) - t124 - t300;
t325 = t99 * mrSges(5,2) - t86 * t230 + t135 * t229 / 0.2e1 - t136 * t228 / 0.2e1 - t181 * mrSges(4,1);
t323 = -0.2e1 * pkin(1);
t321 = Ifges(7,4) * t317 + Ifges(7,2) * t316 + Ifges(7,6) * t311;
t320 = Ifges(7,1) * t317 + Ifges(7,4) * t316 + Ifges(7,5) * t311;
t156 = qJD(6) + t163;
t75 = Ifges(7,4) * t233;
t27 = -Ifges(7,1) * t80 + Ifges(7,5) * t156 + t75;
t318 = t27 / 0.2e1;
t315 = -t233 / 0.2e1;
t314 = t233 / 0.2e1;
t313 = t80 / 0.2e1;
t312 = -t80 / 0.2e1;
t310 = -t156 / 0.2e1;
t309 = t156 / 0.2e1;
t308 = -t161 / 0.2e1;
t305 = t163 / 0.2e1;
t304 = t204 / 0.2e1;
t303 = t206 / 0.2e1;
t302 = Ifges(7,4) * t80;
t301 = pkin(2) * t205;
t299 = pkin(7) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t249);
t297 = pkin(8) * t206;
t291 = t233 * Ifges(7,6);
t290 = t80 * Ifges(7,5);
t289 = -qJD(2) / 0.2e1;
t288 = -Ifges(4,4) - Ifges(5,6);
t173 = t205 * t208 - t235;
t247 = qJD(2) * t208;
t164 = qJD(2) * t235 - t205 * t247;
t214 = pkin(2) * t247 - qJ(4) * t164 - qJD(4) * t175;
t49 = qJD(5) * t173 + t162 * t287 + t214;
t111 = t157 * t205 - t265 * t158;
t78 = pkin(4) * t164 + t111;
t23 = t204 * t78 + t206 * t49;
t277 = Ifges(6,5) * t204;
t276 = Ifges(6,6) * t206;
t133 = -t265 * t184 - t185 * t205;
t275 = t133 * t94;
t272 = t152 * mrSges(5,1);
t271 = t156 * Ifges(7,3);
t269 = t163 * Ifges(4,4);
t268 = t163 * Ifges(5,6);
t113 = pkin(4) * t175 + t133;
t219 = -t175 * qJ(4) + t197;
t93 = t173 * t287 + t219;
t46 = t204 * t113 + t206 * t93;
t264 = Ifges(3,5) * qJD(2);
t257 = t163 * t204;
t95 = t265 * t145 + t205 * t211;
t244 = -Ifges(4,4) / 0.2e1 - Ifges(5,6) / 0.2e1;
t243 = qJD(2) * qJD(4) + t95;
t242 = -pkin(5) * t206 - pkin(4);
t14 = -t42 * mrSges(7,1) + t41 * mrSges(7,2);
t239 = t210 * t246;
t192 = qJ(4) + t301;
t98 = -mrSges(6,1) * t258 + mrSges(6,2) * t259;
t227 = t276 + t277;
t225 = -t17 * t204 + t18 * t206;
t224 = t204 * t33 + t206 * t32;
t223 = t204 * t32 - t206 * t33;
t107 = t206 * t113;
t31 = pkin(5) * t175 + t107 + (-pkin(8) * t173 - t93) * t204;
t35 = t173 * t297 + t46;
t9 = -t207 * t35 + t209 * t31;
t10 = t207 * t31 + t209 * t35;
t134 = t205 * t184 - t185 * t265;
t222 = t133 * t152 - t134 * t151;
t218 = t277 / 0.2e1 + t276 / 0.2e1;
t112 = t157 * t265 + t205 * t158;
t131 = t178 * t265 + t167;
t116 = t331 * t173;
t198 = Ifges(3,4) * t248;
t180 = pkin(5) * t204 + t192;
t171 = Ifges(3,1) * t249 + t198 + t264;
t155 = Ifges(4,4) * t161;
t154 = Ifges(5,6) * t161;
t149 = Ifges(7,3) * t152;
t148 = t152 * mrSges(5,3);
t147 = t152 * mrSges(4,2);
t127 = t173 * pkin(3) + t219;
t126 = -mrSges(5,2) * t161 - mrSges(5,3) * t163;
t125 = mrSges(4,1) * t161 + mrSges(4,2) * t163;
t123 = t163 * Ifges(4,1) + Ifges(4,5) * qJD(2) - t155;
t122 = -t161 * Ifges(4,2) + Ifges(4,6) * qJD(2) + t269;
t121 = Ifges(5,4) * qJD(2) - t163 * Ifges(5,2) + t154;
t120 = Ifges(5,5) * qJD(2) + t161 * Ifges(5,3) - t268;
t117 = t221 * t173;
t115 = -qJD(2) * pkin(3) + t217;
t114 = -t173 * pkin(4) + t134;
t110 = pkin(3) * t163 + t234;
t97 = t131 - t295;
t87 = pkin(3) * t162 + t214;
t79 = -t162 * pkin(4) + t112;
t77 = t173 * t242 + t134;
t73 = t206 * t78;
t67 = pkin(3) * t151 + t215;
t65 = t152 * Ifges(6,5) + t151 * t229;
t64 = t152 * Ifges(6,6) + t151 * t228;
t63 = t163 * t242 + t131;
t62 = -t135 * Ifges(6,1) + t136 * Ifges(6,4) + Ifges(6,5) * t163;
t61 = -t135 * Ifges(6,4) + t136 * Ifges(6,2) + Ifges(6,6) * t163;
t60 = t163 * Ifges(6,3) + t273 - t274;
t59 = -pkin(4) * t151 + t243;
t55 = mrSges(7,1) * t156 + mrSges(7,3) * t80;
t54 = -mrSges(7,2) * t156 + mrSges(7,3) * t233;
t53 = t162 * t242 + t112;
t52 = -pkin(5) * t136 + t86;
t51 = t162 * t331 - t166 * t173;
t50 = qJD(6) * t116 + t162 * t221;
t47 = t151 * t242 + t243;
t45 = -t204 * t93 + t107;
t36 = -t204 * t76 + t92;
t30 = -mrSges(7,2) * t152 + mrSges(7,3) * t42;
t29 = mrSges(7,1) * t152 - mrSges(7,3) * t41;
t26 = Ifges(7,2) * t233 + Ifges(7,6) * t156 - t302;
t25 = t271 - t290 + t291;
t22 = -t204 * t49 + t73;
t19 = t162 * t297 + t23;
t16 = pkin(5) * t164 + t73 + (-pkin(8) * t162 - t49) * t204;
t4 = -qJD(6) * t10 + t16 * t209 - t19 * t207;
t3 = qJD(6) * t9 + t16 * t207 + t19 * t209;
t7 = [-t334 * t111 + (t1 * t116 - t117 * t2 - t5 * t50 + t51 * t6) * mrSges(7,3) + m(4) * (-t111 * t128 + t112 * t129 + t134 * t95 + t275) + m(7) * (t1 * t10 + t2 * t9 + t3 * t6 + t4 * t5 + t47 * t77 + t52 * t53) + m(6) * (t114 * t59 + t17 * t45 + t18 * t46 + t22 * t32 + t23 * t33 + t79 * t86) + (-t128 * t164 - t129 * t162 + t222) * mrSges(4,3) + (t115 * t164 + t124 * t162 + t222) * mrSges(5,1) + m(5) * (t111 * t115 - t112 * t124 + t127 * t67 + t134 * t243 + t87 * t99 + t275) + (t59 * t230 - t67 * mrSges(5,2) + t64 * t303 + t65 * t304 - t243 * mrSges(5,1) - t95 * mrSges(4,3) + (t218 + t288) * t152 + t225 * mrSges(6,3) + (Ifges(4,2) + Ifges(5,3) + Ifges(6,2) * t206 ^ 2 / 0.2e1 + (t278 + t281 / 0.2e1) * t204) * t151) * t173 + t197 * (t151 * mrSges(4,1) + t147) + (t149 / 0.2e1 - t67 * mrSges(5,3) - t18 * mrSges(6,2) + t17 * mrSges(6,1) + (mrSges(4,3) + mrSges(5,1)) * t94 + (Ifges(7,3) / 0.2e1 + Ifges(5,2) - t339) * t152 + (t227 + t288) * t151 + t327) * t175 + t87 * t126 + t114 * t98 + t335 * t112 + t47 * (-mrSges(7,1) * t116 + mrSges(7,2) * t117) + t45 * t108 + t46 * t109 + t23 * t100 + (-t263 / 0.2e1 + (mrSges(3,1) * t323 - 0.3e1 / 0.2e1 * t280 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t210) * qJD(1) + (qJD(1) * (mrSges(4,1) * t173 + mrSges(4,2) * t175) + m(4) * (t181 + t250) + t125) * pkin(2) - t328) * t247 + t22 * t101 + t79 * t88 + t77 * t14 + t51 * t26 / 0.2e1 + t52 * (-mrSges(7,1) * t51 + mrSges(7,2) * t50) + t53 * t34 + t3 * t54 + t4 * t55 + (-t122 / 0.2e1 + t120 / 0.2e1 + t61 * t303 + t62 * t304 + (t218 + t244) * t163 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t161 - t223 * mrSges(6,3) - t325) * t162 + t9 * t29 + t10 * t30 + t127 * (-t151 * mrSges(5,2) - t148) + t116 * t321 + (Ifges(7,1) * t50 + Ifges(7,4) * t51) * t312 + (Ifges(7,4) * t50 + Ifges(7,2) * t51) * t314 + (Ifges(7,4) * t117 + Ifges(7,2) * t116) * t316 + (Ifges(7,1) * t117 + Ifges(7,4) * t116) * t317 + t50 * t318 + t117 * t320 + (Ifges(7,5) * t50 + Ifges(7,6) * t51) * t309 + (Ifges(7,5) * t117 + Ifges(7,6) * t116) * t311 + (-t121 / 0.2e1 + t271 / 0.2e1 + t123 / 0.2e1 + t291 / 0.2e1 - t290 / 0.2e1 + t60 / 0.2e1 + t25 / 0.2e1 + t244 * t161 + (Ifges(6,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t163 + t346) * t164 + ((Ifges(4,5) / 0.2e1 - Ifges(5,4) / 0.2e1) * t164 + (-Ifges(4,6) / 0.2e1 + Ifges(5,5) / 0.2e1) * t162 + (t171 / 0.2e1 - t299 + t264 / 0.2e1 + (mrSges(3,2) * t323 + 0.3e1 / 0.2e1 * Ifges(3,4) * t210) * qJD(1)) * t210) * qJD(2); (Ifges(7,4) * t331 - Ifges(7,2) * t221) * t316 + (Ifges(7,1) * t331 - Ifges(7,4) * t221) * t317 + (Ifges(6,5) * t206 + Ifges(7,5) * t331 - Ifges(6,6) * t204 - Ifges(7,6) * t221) * t311 + t47 * (mrSges(7,1) * t221 + mrSges(7,2) * t331) + t331 * t320 + (mrSges(7,1) * t333 - mrSges(7,2) * t347) * t52 + (Ifges(7,5) * t213 + Ifges(7,6) * t336) * t310 + (Ifges(7,4) * t213 + Ifges(7,2) * t336) * t315 + (Ifges(7,1) * t213 + Ifges(7,4) * t336) * t313 - t125 * t199 - t61 * t256 / 0.2e1 - t62 * t257 / 0.2e1 + Ifges(3,5) * t239 + (t268 + t122) * t305 + (-mrSges(3,1) * t239 + mrSges(3,2) * t240) * pkin(7) + (-Ifges(6,2) * t204 + t278) * t258 / 0.2e1 + (Ifges(6,1) * t206 - t279) * t259 / 0.2e1 + t345 * mrSges(7,3) - t128 * t283 - t124 * t284 - Ifges(3,6) * t240 - (Ifges(3,5) * t210 - Ifges(3,6) * t208) * t246 / 0.2e1 - (-Ifges(3,2) * t249 + t171 + t198) * t248 / 0.2e1 + (-Ifges(4,2) * t163 + t123 - t155 + t25 + t60) * t161 / 0.2e1 + (-t208 * (Ifges(3,1) * t210 - t280) / 0.2e1 + pkin(1) * (mrSges(3,1) * t208 + mrSges(3,2) * t210)) * qJD(1) ^ 2 - t213 * t27 / 0.2e1 + (t128 * t130 - t129 * t131 - t181 * t199 + (t205 * t95 - t265 * t94) * pkin(2)) * m(4) + (-t110 * t99 - t115 * t130 + t124 * t131 + t192 * t243 + t196 * t94) * m(5) + t243 * mrSges(5,3) - t221 * t321 + t59 * (mrSges(6,1) * t204 + mrSges(6,2) * t206) - t204 * t64 / 0.2e1 + t192 * t98 + t180 * t14 + t341 * t55 + t342 * t54 + (t1 * t119 + t118 * t2 + t180 * t47 + t341 * t5 + t342 * t6 - t52 * t63) * m(7) - (t339 * t161 + t163 * t227 + t120 - t269) * t163 / 0.2e1 + t340 * t94 + (-t192 * mrSges(5,1) - mrSges(4,3) * t301 + t337) * t151 + (Ifges(5,3) * t308 + t337 * t289 + t325) * t163 + (-mrSges(4,3) * t241 - t338) * t152 + t118 * t29 + t119 * t30 - t110 * t126 - t333 * t26 / 0.2e1 + t334 * t130 - t335 * t131 + (-t256 * t33 + t257 * t32 - t226) * mrSges(6,3) + t329 * qJD(5) + t330 * t190 + (-m(5) * t124 + m(6) * t86 + m(7) * t52 - t245) * qJD(4) - t37 * t100 + t328 * t249 - t36 * t101 - t95 * mrSges(4,2) - t97 * t88 - t63 * t34 - t166 * t318 + t65 * t303 + t248 * t299 + t196 * t272 + t129 * t282 + t115 * t270 + (-Ifges(7,1) * t166 - Ifges(7,4) * t165) * t312 + (-Ifges(7,4) * t166 - Ifges(7,2) * t165) * t314 + (-Ifges(7,5) * t166 - Ifges(7,6) * t165) * t309 + (-Ifges(7,5) * t313 + Ifges(5,2) * t305 - Ifges(7,6) * t315 - Ifges(7,3) * t310 + t338 * t289 + t346) * t161 + (t154 + t121) * t308 + (-qJD(5) * t224 + t190 * t226 + t192 * t59 - t32 * t36 - t33 * t37 - t86 * t97) * m(6); -t204 * t108 + t206 * t109 - t221 * t29 + t331 * t30 + t147 - t148 - t333 * t55 - t347 * t54 - t340 * t151 + (t334 + t329) * t163 + (t141 - t245) * t161 + (t1 * t331 + t161 * t52 - t2 * t221 - t333 * t5 - t347 * t6) * m(7) + (t161 * t86 - t163 * t224 + t225) * m(6) + (-t115 * t163 - t124 * t161 + t67) * m(5) + (t128 * t163 + t129 * t161 + t191) * m(4); t272 + t221 * t30 + t331 * t29 - t347 * t55 + t333 * t54 + (t100 * t206 - t101 * t204 + t126) * t163 + t245 * qJD(2) + (-qJD(2) * t52 - t345) * m(7) + (-qJD(2) * t86 - t163 * t223 + t226) * m(6) + (qJD(2) * t124 + t163 * t99 + t94) * m(5) + t330; -t136 * t100 - t135 * t101 - t233 * t54 - t80 * t55 + t14 + t98 + (-t233 * t6 - t5 * t80 + t47) * m(7) + (-t135 * t32 - t136 * t33 + t59) * m(6); t149 - t52 * (-mrSges(7,1) * t80 + mrSges(7,2) * t233) + (Ifges(7,1) * t233 + t302) * t313 + t26 * t312 + (Ifges(7,5) * t233 + Ifges(7,6) * t80) * t310 - t5 * t54 + t6 * t55 + (t233 * t5 - t6 * t80) * mrSges(7,3) + (Ifges(7,2) * t80 + t27 + t75) * t315 + t327;];
tauc  = t7(:);
