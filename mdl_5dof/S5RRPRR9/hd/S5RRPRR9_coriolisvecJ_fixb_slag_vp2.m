% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:03
% EndTime: 2019-12-31 20:20:21
% DurationCPUTime: 8.11s
% Computational Cost: add. (6623->500), mult. (17007->705), div. (0->0), fcn. (12163->8), ass. (0->247)
t179 = sin(qJ(2));
t182 = cos(qJ(2));
t234 = sin(pkin(9));
t205 = qJD(1) * t234;
t235 = cos(pkin(9));
t206 = qJD(1) * t235;
t144 = -t179 * t205 + t182 * t206;
t266 = -t144 / 0.2e1;
t320 = Ifges(4,2) * t266;
t178 = sin(qJ(4));
t215 = t234 * pkin(2);
t173 = t215 + pkin(7);
t256 = pkin(8) + t173;
t209 = qJD(4) * t256;
t230 = t144 * t178;
t255 = -qJ(3) - pkin(6);
t168 = t255 * t182;
t162 = qJD(1) * t168;
t148 = t234 * t162;
t167 = t255 * t179;
t161 = qJD(1) * t167;
t116 = t161 * t235 + t148;
t181 = cos(qJ(4));
t145 = -t179 * t206 - t182 * t205;
t226 = qJD(1) * t179;
t220 = pkin(2) * t226;
t97 = -pkin(3) * t145 - pkin(7) * t144 + t220;
t51 = t181 * t116 + t178 * t97;
t319 = -pkin(8) * t230 + t178 * t209 + t51;
t229 = t144 * t181;
t50 = -t116 * t178 + t181 * t97;
t318 = pkin(4) * t145 + pkin(8) * t229 - t181 * t209 - t50;
t223 = qJD(4) * t178;
t317 = t223 - t230;
t180 = cos(qJ(5));
t177 = sin(qJ(5));
t124 = qJD(2) * t181 + t145 * t178;
t153 = qJD(2) * pkin(2) + t161;
t208 = t235 * t162;
t114 = t234 * t153 - t208;
t107 = qJD(2) * pkin(7) + t114;
t175 = -pkin(2) * t182 - pkin(1);
t227 = qJD(1) * t175;
t163 = qJD(3) + t227;
t86 = -t144 * pkin(3) + t145 * pkin(7) + t163;
t48 = t107 * t181 + t178 * t86;
t41 = pkin(8) * t124 + t48;
t240 = t177 * t41;
t139 = qJD(4) - t144;
t125 = qJD(2) * t178 - t145 * t181;
t47 = -t107 * t178 + t181 * t86;
t40 = -pkin(8) * t125 + t47;
t34 = pkin(4) * t139 + t40;
t10 = t180 * t34 - t240;
t238 = t180 * t41;
t11 = t177 * t34 + t238;
t157 = t179 * t235 + t182 * t234;
t146 = t157 * qJD(2);
t135 = qJD(1) * t146;
t132 = Ifges(6,3) * t135;
t204 = t180 * t124 - t125 * t177;
t69 = t124 * t177 + t125 * t180;
t262 = Ifges(6,4) * t69;
t134 = qJD(5) + t139;
t270 = -t134 / 0.2e1;
t279 = -t69 / 0.2e1;
t113 = t153 * t235 + t148;
t106 = -qJD(2) * pkin(3) - t113;
t64 = -t124 * pkin(4) + t106;
t316 = t132 + (Ifges(6,5) * t204 - Ifges(6,6) * t69) * t270 + (t10 * t204 + t11 * t69) * mrSges(6,3) - t64 * (mrSges(6,1) * t69 + mrSges(6,2) * t204) + (Ifges(6,1) * t204 - t262) * t279;
t222 = qJD(4) * t181;
t210 = qJD(2) * t255;
t142 = qJD(3) * t182 + t179 * t210;
t129 = t142 * qJD(1);
t143 = -t179 * qJD(3) + t182 * t210;
t183 = qJD(1) * t143;
t79 = t129 * t235 + t183 * t234;
t156 = t179 * t234 - t182 * t235;
t147 = t156 * qJD(2);
t136 = qJD(1) * t147;
t221 = qJD(1) * qJD(2);
t214 = t179 * t221;
t203 = pkin(2) * t214;
t83 = pkin(3) * t135 + pkin(7) * t136 + t203;
t19 = -t107 * t223 + t178 * t83 + t181 * t79 + t86 * t222;
t81 = -qJD(4) * t125 + t136 * t178;
t12 = pkin(8) * t81 + t19;
t20 = -qJD(4) * t48 - t178 * t79 + t181 * t83;
t80 = qJD(4) * t124 - t136 * t181;
t9 = pkin(4) * t135 - pkin(8) * t80 + t20;
t2 = qJD(5) * t10 + t12 * t180 + t177 * t9;
t27 = t204 * qJD(5) + t177 * t81 + t180 * t80;
t28 = -qJD(5) * t69 - t177 * t80 + t180 * t81;
t3 = -qJD(5) * t11 - t12 * t177 + t180 * t9;
t315 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t27 + Ifges(6,6) * t28;
t65 = Ifges(6,4) * t204;
t314 = -Ifges(6,2) * t69 + t65;
t312 = -Ifges(4,6) / 0.2e1;
t287 = t27 / 0.2e1;
t286 = t28 / 0.2e1;
t268 = t135 / 0.2e1;
t154 = t256 * t178;
t155 = t256 * t181;
t110 = -t154 * t177 + t155 * t180;
t311 = -qJD(5) * t110 + t319 * t177 + t318 * t180;
t109 = -t154 * t180 - t155 * t177;
t310 = qJD(5) * t109 + t318 * t177 - t319 * t180;
t305 = t163 * mrSges(4,2);
t190 = t177 * t178 - t180 * t181;
t299 = qJD(4) + qJD(5);
t118 = t299 * t190;
t93 = t190 * t144;
t304 = -t118 + t93;
t160 = t177 * t181 + t178 * t180;
t119 = t299 * t160;
t92 = t160 * t144;
t303 = -t119 + t92;
t254 = mrSges(4,3) * t145;
t236 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t124 + mrSges(5,2) * t125 - t254;
t115 = t161 * t234 - t208;
t302 = t317 * pkin(4) - t115;
t301 = Ifges(4,5) * qJD(2);
t103 = t190 * t157;
t112 = t156 * pkin(3) - t157 * pkin(7) + t175;
t122 = t167 * t234 - t168 * t235;
t117 = t181 * t122;
t63 = t178 * t112 + t117;
t300 = -t178 * t20 + t181 * t19;
t225 = qJD(1) * t182;
t232 = Ifges(3,6) * qJD(2);
t253 = Ifges(3,4) * t179;
t297 = t232 / 0.2e1 + (t182 * Ifges(3,2) + t253) * qJD(1) / 0.2e1 + pkin(6) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t225);
t296 = qJD(2) * t312;
t241 = t145 * Ifges(4,4);
t295 = t241 / 0.2e1 + t320;
t294 = t20 * mrSges(5,1) - t19 * mrSges(5,2) + Ifges(5,5) * t80 + Ifges(5,6) * t81 + t315;
t201 = mrSges(5,1) * t178 + mrSges(5,2) * t181;
t188 = t106 * t201;
t196 = Ifges(5,5) * t181 - Ifges(5,6) * t178;
t250 = Ifges(5,4) * t181;
t198 = -Ifges(5,2) * t178 + t250;
t251 = Ifges(5,4) * t178;
t200 = Ifges(5,1) * t181 - t251;
t263 = t181 / 0.2e1;
t264 = -t178 / 0.2e1;
t271 = t125 / 0.2e1;
t252 = Ifges(5,4) * t125;
t58 = t124 * Ifges(5,2) + t139 * Ifges(5,6) + t252;
t123 = Ifges(5,4) * t124;
t59 = t125 * Ifges(5,1) + t139 * Ifges(5,5) + t123;
t293 = t139 * t196 / 0.2e1 + t200 * t271 + t124 * t198 / 0.2e1 + t188 + t59 * t263 + t58 * t264;
t292 = t163 * mrSges(4,1) + t47 * mrSges(5,1) + t10 * mrSges(6,1) - t48 * mrSges(5,2) - t11 * mrSges(6,2) + t295 + t296;
t291 = -0.2e1 * pkin(1);
t289 = Ifges(6,4) * t287 + Ifges(6,2) * t286 + Ifges(6,6) * t268;
t288 = Ifges(6,1) * t287 + Ifges(6,4) * t286 + Ifges(6,5) * t268;
t32 = Ifges(6,2) * t204 + Ifges(6,6) * t134 + t262;
t285 = -t32 / 0.2e1;
t284 = t32 / 0.2e1;
t33 = Ifges(6,1) * t69 + Ifges(6,5) * t134 + t65;
t283 = -t33 / 0.2e1;
t282 = t33 / 0.2e1;
t281 = -t204 / 0.2e1;
t280 = t204 / 0.2e1;
t278 = t69 / 0.2e1;
t277 = t80 / 0.2e1;
t276 = t81 / 0.2e1;
t273 = -t124 / 0.2e1;
t272 = -t125 / 0.2e1;
t269 = t134 / 0.2e1;
t267 = -t139 / 0.2e1;
t265 = t145 / 0.2e1;
t261 = pkin(6) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t226);
t259 = pkin(8) * t181;
t258 = t204 * Ifges(6,6);
t257 = t69 * Ifges(6,5);
t138 = Ifges(4,4) * t144;
t121 = -t235 * t167 - t168 * t234;
t78 = t129 * t234 - t235 * t183;
t249 = t121 * t78;
t248 = t124 * Ifges(5,6);
t247 = t125 * Ifges(5,5);
t246 = t134 * Ifges(6,3);
t245 = t139 * Ifges(5,3);
t244 = t144 * mrSges(4,3);
t242 = t145 * Ifges(4,1);
t233 = Ifges(3,5) * qJD(2);
t228 = t157 * t178;
t224 = qJD(2) * t179;
t216 = t235 * pkin(2);
t213 = t182 * t221;
t96 = t142 * t235 + t143 * t234;
t98 = pkin(2) * t224 + pkin(3) * t146 + pkin(7) * t147;
t211 = -t178 * t96 + t181 * t98;
t207 = t135 * mrSges(4,1) - t136 * mrSges(4,2);
t62 = t181 * t112 - t122 * t178;
t174 = -t216 - pkin(3);
t202 = mrSges(5,1) * t181 - mrSges(5,2) * t178;
t199 = Ifges(5,1) * t178 + t250;
t197 = Ifges(5,2) * t181 + t251;
t195 = Ifges(5,5) * t178 + Ifges(5,6) * t181;
t44 = pkin(4) * t156 - t157 * t259 + t62;
t49 = -pkin(8) * t228 + t63;
t21 = -t177 * t49 + t180 * t44;
t22 = t177 * t44 + t180 * t49;
t194 = -t178 * t19 - t181 * t20;
t193 = -t178 * t48 - t181 * t47;
t192 = t178 * t47 - t181 * t48;
t84 = -mrSges(5,2) * t139 + mrSges(5,3) * t124;
t85 = mrSges(5,1) * t139 - mrSges(5,3) * t125;
t191 = -t178 * t85 + t181 * t84;
t95 = t142 * t234 - t235 * t143;
t189 = -t147 * t178 + t157 * t222;
t29 = t112 * t222 - t122 * t223 + t178 * t98 + t181 * t96;
t176 = Ifges(3,4) * t225;
t164 = -t181 * pkin(4) + t174;
t152 = Ifges(3,1) * t226 + t176 + t233;
t133 = Ifges(5,3) * t135;
t127 = -qJD(2) * mrSges(4,2) + t244;
t108 = -mrSges(4,1) * t144 - mrSges(4,2) * t145;
t102 = t160 * t157;
t101 = t138 - t242 + t301;
t94 = pkin(4) * t228 + t121;
t61 = -mrSges(5,2) * t135 + mrSges(5,3) * t81;
t60 = mrSges(5,1) * t135 - mrSges(5,3) * t80;
t57 = t245 + t247 + t248;
t54 = mrSges(6,1) * t134 - mrSges(6,3) * t69;
t53 = -mrSges(6,2) * t134 + mrSges(6,3) * t204;
t52 = pkin(4) * t189 + t95;
t45 = -t81 * pkin(4) + t78;
t43 = -mrSges(5,1) * t81 + mrSges(5,2) * t80;
t42 = -mrSges(6,1) * t204 + mrSges(6,2) * t69;
t38 = t80 * Ifges(5,1) + t81 * Ifges(5,4) + t135 * Ifges(5,5);
t37 = t80 * Ifges(5,4) + t81 * Ifges(5,2) + t135 * Ifges(5,6);
t36 = t299 * t103 + t160 * t147;
t35 = -t119 * t157 + t147 * t190;
t31 = t246 + t257 + t258;
t30 = -t63 * qJD(4) + t211;
t24 = -mrSges(6,2) * t135 + mrSges(6,3) * t28;
t23 = mrSges(6,1) * t135 - mrSges(6,3) * t27;
t18 = -pkin(8) * t189 + t29;
t17 = t147 * t259 + pkin(4) * t146 + (-t117 + (pkin(8) * t157 - t112) * t178) * qJD(4) + t211;
t14 = t180 * t40 - t240;
t13 = -t177 * t40 - t238;
t8 = -mrSges(6,1) * t28 + mrSges(6,2) * t27;
t5 = -qJD(5) * t22 + t17 * t180 - t177 * t18;
t4 = qJD(5) * t21 + t17 * t177 + t18 * t180;
t1 = [(t196 * t268 + t200 * t277 + t198 * t276 - Ifges(4,1) * t136 - Ifges(4,4) * t135 + t37 * t264 + t38 * t263 + (mrSges(4,3) + t201) * t78 + t194 * mrSges(5,3) + (t197 * t273 + t199 * t272 + t106 * t202 + t195 * t267 + t59 * t264 - t181 * t58 / 0.2e1 + t192 * mrSges(5,3)) * qJD(4)) * t157 + (t113 * t147 - t114 * t146 - t121 * t136 - t122 * t135) * mrSges(4,3) + (t132 / 0.2e1 + t133 / 0.2e1 + Ifges(4,4) * t136 - t79 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(4,2)) * t135 + t294) * t156 + m(4) * (-t113 * t95 + t114 * t96 + t122 * t79 + t249) + m(5) * (t106 * t95 + t19 * t63 + t20 * t62 + t29 * t48 + t30 * t47 + t249) + t236 * t95 + t175 * t207 + (-t10 * t35 - t102 * t2 + t103 * t3 + t11 * t36) * mrSges(6,3) + (-Ifges(6,4) * t103 - Ifges(6,2) * t102) * t286 + (-Ifges(6,1) * t103 - Ifges(6,4) * t102) * t287 + t45 * (mrSges(6,1) * t102 - mrSges(6,2) * t103) + (-Ifges(6,5) * t103 - Ifges(6,6) * t102) * t268 - (t305 - t242 / 0.2e1 + t138 / 0.2e1 + t101 / 0.2e1 + t193 * mrSges(5,3) + t293) * t147 - t103 * t288 - t102 * t289 + t35 * t282 + t36 * t284 + (Ifges(6,4) * t35 + Ifges(6,2) * t36) * t280 + (-Ifges(4,5) * t147 / 0.2e1 + t146 * t312 + (t152 / 0.2e1 - t261 + t233 / 0.2e1 + (mrSges(3,2) * t291 + 0.3e1 / 0.2e1 * Ifges(3,4) * t182) * qJD(1)) * t182) * qJD(2) + (-t232 / 0.2e1 + (mrSges(3,1) * t291 - 0.3e1 / 0.2e1 * t253 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t182) * qJD(1) + (qJD(1) * (mrSges(4,1) * t156 + mrSges(4,2) * t157) + t108 + m(4) * (t163 + t227)) * pkin(2) - t297) * t224 + m(6) * (t10 * t5 + t11 * t4 + t2 * t22 + t21 * t3 + t45 * t94 + t52 * t64) + (t248 / 0.2e1 + t247 / 0.2e1 + t245 / 0.2e1 + t246 / 0.2e1 + t258 / 0.2e1 + t257 / 0.2e1 + t57 / 0.2e1 + t31 / 0.2e1 + t292 + t295) * t146 + t96 * t127 + t121 * t43 + t94 * t8 + t29 * t84 + t30 * t85 + t62 * t60 + t63 * t61 + t64 * (-mrSges(6,1) * t36 + mrSges(6,2) * t35) + t4 * t53 + t5 * t54 + t52 * t42 + t22 * t24 + t21 * t23 + (Ifges(6,5) * t35 + Ifges(6,6) * t36) * t269 + (Ifges(6,1) * t35 + Ifges(6,4) * t36) * t278; (m(5) * t174 - mrSges(4,1) - t202) * t78 - t114 * t254 + (Ifges(6,1) * t160 - Ifges(6,4) * t190) * t287 + t45 * (mrSges(6,1) * t190 + mrSges(6,2) * t160) - t190 * t289 + (-t304 * t10 + t303 * t11 - t160 * t3 - t190 * t2) * mrSges(6,3) + (-Ifges(6,5) * t118 - Ifges(6,6) * t119) * t269 + (-Ifges(6,1) * t118 - Ifges(6,4) * t119) * t278 + (-Ifges(6,4) * t118 - Ifges(6,2) * t119) * t280 + (-t135 * t215 + t136 * t216) * mrSges(4,3) + (Ifges(6,5) * t160 - Ifges(6,6) * t190 + t195) * t268 + (Ifges(6,4) * t160 - Ifges(6,2) * t190) * t286 - t59 * t229 / 0.2e1 + t58 * t230 / 0.2e1 - (Ifges(3,5) * t182 - Ifges(3,6) * t179) * t221 / 0.2e1 - t108 * t220 - Ifges(3,6) * t214 + (t241 + t57 + t31) * t265 + (-mrSges(3,1) * t213 + mrSges(3,2) * t214) * pkin(6) + (t138 + t101) * t266 + ((t234 * t79 - t235 * t78) * pkin(2) + t113 * t115 - t114 * t116 - t163 * t220) * m(4) + (-Ifges(6,4) * t93 - Ifges(6,2) * t92) * t281 + (-Ifges(6,1) * t93 - Ifges(6,4) * t92) * t279 + (-Ifges(6,5) * t93 - Ifges(6,6) * t92) * t270 - (-Ifges(3,2) * t226 + t152 + t176) * t225 / 0.2e1 + (-t179 * (Ifges(3,1) * t182 - t253) / 0.2e1 + pkin(1) * (mrSges(3,1) * t179 + mrSges(3,2) * t182)) * qJD(1) ^ 2 - Ifges(4,6) * t135 - Ifges(4,5) * t136 + t310 * t53 + (t311 * t10 + t109 * t3 + t310 * t11 + t110 * t2 + t164 * t45 + t302 * t64) * m(6) + t311 * t54 + (-t188 - t301 / 0.2e1 - t305 + Ifges(4,1) * t265 + t196 * t267 + t200 * t272 + t198 * t273) * t144 - t118 * t282 - t93 * t283 - t119 * t284 - t92 * t285 + t160 * t288 + (-t303 * mrSges(6,1) + t304 * mrSges(6,2)) * t64 + t302 * t42 + (-m(5) * t106 - t236) * t115 + (m(5) * t173 * t193 + t293) * qJD(4) + t297 * t226 + (m(5) * t300 - t178 * t60 + t181 * t61 - t85 * t222 - t84 * t223) * t173 + (-Ifges(5,5) * t272 - Ifges(6,5) * t279 - Ifges(5,6) * t273 - Ifges(6,6) * t281 - Ifges(5,3) * t267 - Ifges(6,3) * t270 + t292 + t296 + t320) * t145 + (-t317 * t48 + (-t222 + t229) * t47 + t300) * mrSges(5,3) - t116 * t127 + t178 * t38 / 0.2e1 + t174 * t43 + t164 * t8 + t109 * t23 + t110 * t24 - t51 * t84 - t50 * t85 - t79 * mrSges(4,2) + t113 * t244 + Ifges(3,5) * t213 + t225 * t261 + t37 * t263 + t197 * t276 + t199 * t277 - m(5) * (t47 * t50 + t48 * t51); -t190 * t23 + t160 * t24 + t178 * t61 + t181 * t60 + t303 * t54 + t304 * t53 + t191 * qJD(4) + (t42 + t236) * t145 + (-t127 - t191) * t144 + t207 + (t303 * t10 + t304 * t11 + t145 * t64 + t160 * t2 - t190 * t3) * m(6) + (t106 * t145 - t139 * t192 - t194) * m(5) + (-t113 * t145 - t114 * t144 + t203) * m(4); -t69 * t285 + t294 + (-t125 * t42 + t177 * t24 + t180 * t23 + (-t177 * t54 + t180 * t53) * qJD(5) + (-t125 * t64 + t177 * t2 + t180 * t3 + (-t10 * t177 + t11 * t180) * qJD(5)) * m(6)) * pkin(4) + t204 * t283 + t133 + t314 * t281 - t106 * (mrSges(5,1) * t125 + mrSges(5,2) * t124) - m(6) * (t10 * t13 + t11 * t14) - t47 * t84 + t48 * t85 - t14 * t53 - t13 * t54 + (-Ifges(5,2) * t125 + t123 + t59) * t273 + (Ifges(5,5) * t124 - Ifges(5,6) * t125) * t267 + t58 * t271 + (Ifges(5,1) * t124 - t252) * t272 + (t124 * t47 + t125 * t48) * mrSges(5,3) + t316; t32 * t278 - t10 * t53 + t11 * t54 + (t314 + t33) * t281 + t315 + t316;];
tauc = t1(:);
