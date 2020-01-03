% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:46
% EndTime: 2019-12-31 19:57:02
% DurationCPUTime: 7.34s
% Computational Cost: add. (3592->415), mult. (9526->552), div. (0->0), fcn. (6428->6), ass. (0->198)
t293 = Ifges(5,4) + Ifges(6,4);
t294 = Ifges(5,1) + Ifges(6,1);
t286 = Ifges(6,5) + Ifges(5,5);
t292 = Ifges(5,2) + Ifges(6,2);
t291 = Ifges(5,6) + Ifges(6,6);
t144 = cos(qJ(4));
t296 = t293 * t144;
t142 = sin(qJ(4));
t295 = t293 * t142;
t143 = sin(qJ(2));
t145 = cos(qJ(2));
t219 = sin(pkin(8));
t220 = cos(pkin(8));
t125 = t220 * t143 + t219 * t145;
t113 = t125 * qJD(1);
t93 = qJD(2) * t144 - t113 * t142;
t290 = t293 * t93;
t269 = -t291 * t142 + t286 * t144;
t268 = -t292 * t142 + t296;
t267 = t294 * t144 - t295;
t205 = qJD(4) * t142;
t148 = -t143 * t219 + t145 * t220;
t112 = t148 * qJD(1);
t213 = t112 * t142;
t289 = t205 - t213;
t94 = qJD(2) * t142 + t113 * t144;
t288 = t293 * t94;
t247 = -t112 / 0.2e1;
t287 = Ifges(4,2) * t247 - Ifges(4,6) * qJD(2) / 0.2e1;
t281 = qJD(4) - t112;
t274 = t281 * t291 + t292 * t93 + t288;
t273 = t286 * t281 + t294 * t94 + t290;
t284 = t142 * t286 + t144 * t291;
t283 = t144 * t292 + t295;
t282 = t142 * t294 + t296;
t140 = -pkin(2) * t145 - pkin(1);
t209 = qJD(1) * t140;
t129 = qJD(3) + t209;
t217 = Ifges(4,5) * qJD(2);
t280 = -t217 / 0.2e1 - t129 * mrSges(4,2);
t114 = t125 * qJD(2);
t105 = qJD(1) * t114;
t115 = t148 * qJD(2);
t106 = qJD(1) * t115;
t55 = qJD(4) * t93 + t106 * t144;
t56 = -qJD(4) * t94 - t106 * t142;
t278 = t291 * t105 + t292 * t56 + t293 * t55;
t277 = t286 * t105 + t293 * t56 + t294 * t55;
t195 = t219 * pkin(2);
t138 = t195 + pkin(7);
t210 = qJ(5) + t138;
t180 = qJD(4) * t210;
t212 = t112 * t144;
t208 = qJD(1) * t143;
t201 = pkin(2) * t208;
t72 = pkin(3) * t113 - pkin(7) * t112 + t201;
t234 = -qJ(3) - pkin(6);
t134 = t234 * t145;
t128 = qJD(1) * t134;
t116 = t219 * t128;
t133 = t234 * t143;
t127 = qJD(1) * t133;
t87 = t127 * t220 + t116;
t25 = -t142 * t87 + t144 * t72;
t276 = -pkin(4) * t113 + qJ(5) * t212 - qJD(5) * t142 - t144 * t180 - t25;
t26 = t142 * t72 + t144 * t87;
t275 = qJ(5) * t213 + qJD(5) * t144 - t142 * t180 - t26;
t221 = t113 * mrSges(4,3);
t233 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t93 + mrSges(5,2) * t94 + t221;
t184 = t220 * t128;
t86 = t127 * t219 - t184;
t271 = t289 * pkin(4) - t86;
t172 = mrSges(6,1) * t142 + mrSges(6,2) * t144;
t174 = mrSges(5,1) * t142 + mrSges(5,2) * t144;
t121 = qJD(2) * pkin(2) + t127;
t84 = t121 * t220 + t116;
t79 = -qJD(2) * pkin(3) - t84;
t40 = -t93 * pkin(4) + qJD(5) + t79;
t270 = -t40 * t172 - t79 * t174;
t204 = qJD(4) * t144;
t266 = -t204 + t212;
t185 = qJD(2) * t234;
t110 = qJD(3) * t145 + t143 * t185;
t100 = t110 * qJD(1);
t111 = -t143 * qJD(3) + t145 * t185;
t146 = qJD(1) * t111;
t54 = t100 * t220 + t146 * t219;
t203 = qJD(1) * qJD(2);
t191 = t143 * t203;
t179 = pkin(2) * t191;
t58 = pkin(3) * t105 - pkin(7) * t106 + t179;
t63 = -t112 * pkin(3) - t113 * pkin(7) + t129;
t85 = t219 * t121 - t184;
t80 = qJD(2) * pkin(7) + t85;
t5 = t142 * t58 + t144 * t54 + t63 * t204 - t205 * t80;
t23 = t142 * t63 + t144 * t80;
t6 = -qJD(4) * t23 - t142 * t54 + t144 * t58;
t264 = -t142 * t6 + t144 * t5;
t207 = qJD(1) * t145;
t216 = Ifges(3,6) * qJD(2);
t228 = Ifges(3,4) * t143;
t262 = pkin(6) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t207) + t216 / 0.2e1 + (Ifges(3,2) * t145 + t228) * qJD(1) / 0.2e1;
t227 = Ifges(4,4) * t113;
t255 = -t227 / 0.2e1 + t287;
t251 = t94 / 0.2e1;
t261 = -t270 + t268 * t93 / 0.2e1 + t267 * t251 + t269 * t281 / 0.2e1;
t16 = qJ(5) * t93 + t23;
t22 = -t142 * t80 + t144 * t63;
t15 = -qJ(5) * t94 + t22;
t9 = pkin(4) * t281 + t15;
t260 = t129 * mrSges(4,1) + t22 * mrSges(5,1) + t9 * mrSges(6,1) - t23 * mrSges(5,2) - t16 * mrSges(6,2) + t255;
t259 = -0.2e1 * pkin(1);
t257 = t55 / 0.2e1;
t256 = t56 / 0.2e1;
t254 = -t93 / 0.2e1;
t252 = -t94 / 0.2e1;
t250 = t105 / 0.2e1;
t249 = -t281 / 0.2e1;
t246 = -t113 / 0.2e1;
t245 = -t142 / 0.2e1;
t242 = t144 / 0.2e1;
t241 = pkin(6) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t208);
t53 = t100 * t219 - t220 * t146;
t89 = -t220 * t133 - t134 * t219;
t237 = t53 * t89;
t83 = -pkin(3) * t148 - t125 * pkin(7) + t140;
t90 = t133 * t219 - t134 * t220;
t88 = t144 * t90;
t39 = t142 * t83 + t88;
t230 = mrSges(4,3) * t112;
t229 = Ifges(4,1) * t113;
t108 = Ifges(4,4) * t112;
t218 = Ifges(3,5) * qJD(2);
t214 = qJ(5) * t125;
t211 = t125 * t142;
t206 = qJD(2) * t143;
t71 = t110 * t220 + t111 * t219;
t73 = pkin(2) * t206 + pkin(3) * t114 - pkin(7) * t115;
t202 = t142 * t73 + t144 * t71 + t83 * t204;
t199 = Ifges(5,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t198 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t197 = Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t196 = t220 * pkin(2);
t194 = t125 * t204;
t17 = -t56 * mrSges(6,1) + t55 * mrSges(6,2);
t190 = t145 * t203;
t186 = -t142 * t71 + t144 * t73;
t38 = -t142 * t90 + t144 * t83;
t181 = t105 * mrSges(4,1) + t106 * mrSges(4,2);
t139 = -t196 - pkin(3);
t1 = pkin(4) * t105 - qJ(5) * t55 - qJD(5) * t94 + t6;
t2 = qJ(5) * t56 + qJD(5) * t93 + t5;
t178 = -t1 * t144 - t2 * t142;
t177 = -t142 * t5 - t144 * t6;
t176 = t142 * t9 - t144 * t16;
t175 = mrSges(5,1) * t144 - mrSges(5,2) * t142;
t173 = mrSges(6,1) * t144 - mrSges(6,2) * t142;
t159 = -t142 * t23 - t144 * t22;
t158 = t142 * t22 - t144 * t23;
t70 = t110 * t219 - t220 * t111;
t157 = -qJ(5) * t115 - qJD(5) * t125;
t147 = t6 * mrSges(5,1) + t1 * mrSges(6,1) - t5 * mrSges(5,2) - t2 * mrSges(6,2);
t141 = Ifges(3,4) * t207;
t130 = -t144 * pkin(4) + t139;
t123 = t210 * t144;
t122 = t210 * t142;
t120 = Ifges(3,1) * t208 + t141 + t218;
t104 = Ifges(5,3) * t105;
t103 = Ifges(6,3) * t105;
t98 = -qJD(2) * mrSges(4,2) + t230;
t81 = -mrSges(4,1) * t112 + mrSges(4,2) * t113;
t76 = t108 + t217 + t229;
t69 = pkin(4) * t211 + t89;
t62 = mrSges(5,1) * t281 - mrSges(5,3) * t94;
t61 = mrSges(6,1) * t281 - mrSges(6,3) * t94;
t60 = -mrSges(5,2) * t281 + mrSges(5,3) * t93;
t59 = -mrSges(6,2) * t281 + mrSges(6,3) * t93;
t51 = Ifges(5,5) * t55;
t50 = Ifges(6,5) * t55;
t49 = Ifges(5,6) * t56;
t48 = Ifges(6,6) * t56;
t41 = -mrSges(6,1) * t93 + mrSges(6,2) * t94;
t37 = -mrSges(5,2) * t105 + mrSges(5,3) * t56;
t36 = -mrSges(6,2) * t105 + mrSges(6,3) * t56;
t35 = mrSges(5,1) * t105 - mrSges(5,3) * t55;
t34 = mrSges(6,1) * t105 - mrSges(6,3) * t55;
t29 = Ifges(5,5) * t94 + Ifges(5,6) * t93 + Ifges(5,3) * t281;
t28 = Ifges(6,5) * t94 + Ifges(6,6) * t93 + Ifges(6,3) * t281;
t27 = (t115 * t142 + t194) * pkin(4) + t70;
t24 = -qJ(5) * t211 + t39;
t21 = -t56 * pkin(4) + t53;
t19 = -pkin(4) * t148 - t144 * t214 + t38;
t18 = -mrSges(5,1) * t56 + mrSges(5,2) * t55;
t8 = -qJD(4) * t39 + t186;
t7 = -t205 * t90 + t202;
t4 = -qJ(5) * t194 + (-qJD(4) * t90 + t157) * t142 + t202;
t3 = pkin(4) * t114 + t157 * t144 + (-t88 + (-t83 + t214) * t142) * qJD(4) + t186;
t10 = [t24 * t36 + t38 * t35 + t39 * t37 + t27 * t41 + t19 * t34 + t4 * t59 + t7 * t60 + t3 * t61 + t8 * t62 + t69 * t17 + t89 * t18 + t71 * t98 + t140 * t181 + t233 * t70 + (-t105 * t90 + t106 * t89) * mrSges(4,3) + m(4) * (t54 * t90 - t70 * t84 + t71 * t85 + t237) + m(5) * (t22 * t8 + t23 * t7 + t38 * t6 + t39 * t5 + t70 * t79 + t237) + m(6) * (t1 * t19 + t16 * t4 + t2 * t24 + t21 * t69 + t27 * t40 + t3 * t9) + (-t241 + t120 / 0.2e1 + t218 / 0.2e1 + (mrSges(3,2) * t259 + 0.3e1 / 0.2e1 * Ifges(3,4) * t145) * qJD(1)) * t145 * qJD(2) + (-t85 * mrSges(4,3) + t28 / 0.2e1 + t29 / 0.2e1 + t199 * t94 + t198 * t93 + t197 * t281 + t260 + t255) * t114 + (-t84 * mrSges(4,3) + t261 + (-t142 * t16 - t144 * t9) * mrSges(6,3) + t274 * t245 + t273 * t242 + t159 * mrSges(5,3) + t108 / 0.2e1 + t76 / 0.2e1 + t229 / 0.2e1 - t280) * t115 - (-t54 * mrSges(4,3) + t50 / 0.2e1 + t48 / 0.2e1 + t103 / 0.2e1 + t51 / 0.2e1 + t49 / 0.2e1 + t104 / 0.2e1 - Ifges(4,4) * t106 + t198 * t56 + t199 * t55 + (Ifges(4,2) + t197) * t105 + t147) * t148 + (t21 * t172 - Ifges(4,4) * t105 + Ifges(4,1) * t106 + (mrSges(4,3) + t174) * t53 + t178 * mrSges(6,3) + t177 * mrSges(5,3) + (mrSges(5,3) * t158 + mrSges(6,3) * t176 + t173 * t40 + t175 * t79 + t283 * t254 + t282 * t252 + t284 * t249 - t274 * t144 / 0.2e1) * qJD(4) + t267 * t257 + t268 * t256 + t269 * t250 + t277 * t242 + (qJD(4) * t273 + t278) * t245) * t125 + (-t216 / 0.2e1 + (mrSges(3,1) * t259 - 0.3e1 / 0.2e1 * t228 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t145) * qJD(1) + (t81 + qJD(1) * (-mrSges(4,1) * t148 + mrSges(4,2) * t125) + m(4) * (t129 + t209)) * pkin(2) - t262) * t206; (m(5) * t138 * t159 + t261) * qJD(4) - m(5) * (t22 * t25 + t23 * t26) + (-mrSges(3,1) * t190 + mrSges(3,2) * t191) * pkin(6) + (t204 / 0.2e1 - t212 / 0.2e1) * t273 + (-t227 + t29 + t28) * t246 + t84 * t230 + t207 * t241 + (m(5) * t139 - mrSges(4,1) - t175) * t53 + (t108 + t76) * t247 + (-t105 * t195 - t106 * t196) * mrSges(4,3) + (-t205 / 0.2e1 + t213 / 0.2e1) * t274 + (-t1 * t142 + t144 * t2 - t16 * t289 + t266 * t9) * mrSges(6,3) + (t22 * t266 - t23 * t289 + t264) * mrSges(5,3) - (-t291 * t254 - t286 * t252 + (-Ifges(6,3) - Ifges(5,3)) * t249 + t260 + t287) * t113 + Ifges(3,5) * t190 - t21 * t173 + ((t219 * t54 - t220 * t53) * pkin(2) - t129 * t201 + t84 * t86 - t85 * t87) * m(4) + t85 * t221 + (-m(5) * t79 - t233) * t86 + t275 * t59 + (-t1 * t122 + t123 * t2 + t130 * t21 + t16 * t275 + t271 * t40 + t276 * t9) * m(6) + t276 * t61 + t277 * t142 / 0.2e1 + t278 * t242 + t271 * t41 + (m(5) * t264 - t142 * t35 + t144 * t37 - t62 * t204 - t60 * t205) * t138 + t262 * t208 - t54 * mrSges(4,2) - t26 * t60 - t25 * t62 + (Ifges(4,1) * t246 + t269 * t249 + t267 * t252 + t268 * t254 + t270 + t280) * t112 + t282 * t257 + t283 * t256 + t284 * t250 - t87 * t98 - Ifges(4,6) * t105 + Ifges(4,5) * t106 - t122 * t34 + t123 * t36 + t130 * t17 - (-Ifges(3,2) * t208 + t120 + t141) * t207 / 0.2e1 + t139 * t18 + (pkin(1) * (mrSges(3,1) * t143 + mrSges(3,2) * t145) - t143 * (Ifges(3,1) * t145 - t228) / 0.2e1) * qJD(1) ^ 2 - Ifges(3,6) * t191 - t81 * t201 - (Ifges(3,5) * t145 - Ifges(3,6) * t143) * t203 / 0.2e1; -t112 * t98 - (t41 + t233) * t113 + (t34 + t35 + t281 * (t59 + t60)) * t144 + (t36 + t37 - t281 * (t61 + t62)) * t142 + t181 + (-t40 * t113 - t176 * t281 - t178) * m(6) + (-t113 * t79 - t158 * t281 - t177) * m(5) + (-t112 * t85 + t113 * t84 + t179) * m(4); t147 + (-t41 * t94 + t34) * pkin(4) + (t16 * t94 + t9 * t93) * mrSges(6,3) + (t22 * t93 + t23 * t94) * mrSges(5,3) + (-(t15 - t9) * t16 + (-t40 * t94 + t1) * pkin(4)) * m(6) + t50 + t49 + t48 + t51 - t15 * t59 - t22 * t60 + t16 * t61 + t23 * t62 - t40 * (mrSges(6,1) * t94 + mrSges(6,2) * t93) - t79 * (mrSges(5,1) * t94 + mrSges(5,2) * t93) + t103 + t104 + (t294 * t93 - t288) * t252 + t274 * t251 + (t286 * t93 - t291 * t94) * t249 + (-t292 * t94 + t273 + t290) * t254; -t93 * t59 + t94 * t61 + 0.2e1 * (t21 / 0.2e1 + t16 * t254 + t9 * t251) * m(6) + t17;];
tauc = t10(:);
