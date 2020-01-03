% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR16_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR16_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:44
% EndTime: 2019-12-31 20:45:01
% DurationCPUTime: 7.22s
% Computational Cost: add. (5335->555), mult. (14123->775), div. (0->0), fcn. (9871->8), ass. (0->265)
t179 = sin(qJ(5));
t182 = cos(qJ(5));
t212 = mrSges(6,1) * t179 + mrSges(6,2) * t182;
t181 = sin(qJ(2));
t177 = sin(pkin(5));
t245 = qJD(1) * t177;
t232 = t181 * t245;
t159 = qJD(4) + t232;
t180 = sin(qJ(4));
t183 = cos(qJ(4));
t184 = cos(qJ(2));
t178 = cos(pkin(5));
t244 = qJD(1) * t178;
t236 = pkin(1) * t244;
t168 = t184 * t236;
t170 = qJD(2) + t244;
t185 = -pkin(2) - pkin(8);
t289 = pkin(3) + pkin(7);
t69 = t170 * t185 + t232 * t289 + qJD(3) - t168;
t225 = -qJ(3) * t181 - pkin(1);
t105 = (t184 * t185 + t225) * t177;
t89 = qJD(1) * t105;
t41 = t180 * t69 + t183 * t89;
t32 = pkin(9) * t159 + t41;
t231 = t184 * t245;
t126 = -t170 * t180 - t183 * t231;
t222 = t180 * t231;
t127 = t170 * t183 - t222;
t140 = pkin(7) * t231 + t181 * t236;
t113 = pkin(3) * t231 + t140;
t161 = t170 * qJ(3);
t83 = t161 + t113;
t42 = -pkin(4) * t126 - pkin(9) * t127 + t83;
t10 = t179 * t42 + t182 * t32;
t9 = -t179 * t32 + t182 * t42;
t214 = t10 * t179 + t182 * t9;
t280 = t182 / 0.2e1;
t281 = -t179 / 0.2e1;
t123 = qJD(5) - t126;
t75 = t127 * t182 + t159 * t179;
t278 = Ifges(6,4) * t75;
t74 = -t127 * t179 + t159 * t182;
t29 = Ifges(6,2) * t74 + Ifges(6,6) * t123 + t278;
t71 = Ifges(6,4) * t74;
t30 = Ifges(6,1) * t75 + Ifges(6,5) * t123 + t71;
t40 = -t180 * t89 + t183 * t69;
t31 = -pkin(4) * t159 - t40;
t313 = -t214 * mrSges(6,3) + t212 * t31 + t280 * t30 + t281 * t29;
t309 = -t170 / 0.2e1;
t312 = mrSges(4,1) + mrSges(3,3);
t311 = mrSges(4,2) - mrSges(3,1);
t207 = Ifges(6,5) * t182 - Ifges(6,6) * t179;
t264 = Ifges(6,4) * t182;
t209 = -Ifges(6,2) * t179 + t264;
t265 = Ifges(6,4) * t179;
t211 = Ifges(6,1) * t182 - t265;
t284 = t123 / 0.2e1;
t291 = t75 / 0.2e1;
t293 = t74 / 0.2e1;
t310 = t207 * t284 + t209 * t293 + t211 * t291 + t313;
t308 = -t245 / 0.2e1;
t121 = Ifges(5,4) * t126;
t263 = Ifges(5,5) * t159;
t267 = Ifges(5,1) * t127;
t60 = t121 + t263 + t267;
t307 = t40 * mrSges(5,3) - t60 / 0.2e1 - t263 / 0.2e1 - t83 * mrSges(5,2) - t121 / 0.2e1;
t239 = qJD(4) * t183;
t240 = qJD(4) * t180;
t226 = qJD(2) * t245;
t221 = t181 * t226;
t157 = pkin(2) * t221;
t205 = pkin(8) * t181 - qJ(3) * t184;
t241 = qJD(3) * t181;
t191 = (qJD(2) * t205 - t241) * t177;
t72 = qJD(1) * t191 + t157;
t174 = t178 * t181 * pkin(1);
t250 = t177 * t184;
t114 = (t250 * t289 + t174) * qJD(2);
t92 = qJD(1) * t114;
t13 = t180 * t92 + t183 * t72 + t69 * t239 - t240 * t89;
t220 = t184 * t226;
t11 = pkin(9) * t220 + t13;
t158 = qJD(2) * t168;
t160 = t170 * qJD(3);
t251 = t177 * t181;
t224 = t289 * t251;
t198 = qJD(2) * t224;
t73 = -qJD(1) * t198 + t158 + t160;
t243 = qJD(2) * t181;
t228 = t180 * t243;
t80 = -t170 * t240 + (-t184 * t239 + t228) * t245;
t81 = -qJD(4) * t222 + t170 * t239 - t183 * t221;
t27 = pkin(4) * t81 - pkin(9) * t80 + t73;
t1 = qJD(5) * t9 + t11 * t182 + t179 * t27;
t2 = -qJD(5) * t10 - t11 * t179 + t182 * t27;
t216 = t1 * t182 - t179 * t2;
t139 = pkin(7) * t232 - t168;
t306 = -qJD(3) - t139;
t14 = -qJD(4) * t41 - t180 * t72 + t183 * t92;
t305 = t14 * mrSges(5,1) - t13 * mrSges(5,2) + Ifges(5,5) * t80 - Ifges(5,6) * t81;
t304 = -Ifges(3,4) * t231 / 0.2e1 + Ifges(3,5) * t309;
t171 = pkin(7) * t251;
t275 = pkin(1) * t184;
t233 = -pkin(2) - t275;
t91 = pkin(3) * t251 + t171 + (-pkin(8) + t233) * t178;
t268 = t183 * t105 + t180 * t91;
t230 = t177 * t243;
t165 = pkin(2) * t230;
t87 = t165 + t191;
t22 = -qJD(4) * t268 + t114 * t183 - t180 * t87;
t37 = -qJD(5) * t75 - t179 * t80 + t182 * t220;
t33 = Ifges(6,6) * t37;
t36 = qJD(5) * t74 + t179 * t220 + t182 * t80;
t34 = Ifges(6,5) * t36;
t5 = Ifges(6,3) * t81 + t33 + t34;
t302 = t5 / 0.2e1;
t7 = Ifges(6,1) * t36 + Ifges(6,4) * t37 + Ifges(6,5) * t81;
t301 = t7 / 0.2e1;
t300 = Ifges(3,5) / 0.2e1;
t299 = Ifges(4,5) / 0.2e1;
t298 = Ifges(5,2) / 0.2e1;
t297 = -t29 / 0.2e1;
t296 = t36 / 0.2e1;
t295 = t37 / 0.2e1;
t294 = -t74 / 0.2e1;
t292 = -t75 / 0.2e1;
t290 = t81 / 0.2e1;
t288 = pkin(1) * mrSges(3,1);
t287 = pkin(1) * mrSges(3,2);
t285 = -t123 / 0.2e1;
t143 = t178 * t180 + t183 * t250;
t283 = -t143 / 0.2e1;
t234 = t180 * t250;
t144 = t178 * t183 - t234;
t282 = t144 / 0.2e1;
t65 = mrSges(5,1) * t220 - mrSges(5,3) * t80;
t8 = -mrSges(6,1) * t37 + mrSges(6,2) * t36;
t279 = t65 - t8;
t277 = Ifges(6,5) * t75;
t276 = Ifges(6,6) * t74;
t272 = t80 * Ifges(5,1);
t271 = t80 * Ifges(5,4);
t270 = t81 * Ifges(5,4);
t43 = -mrSges(6,1) * t74 + mrSges(6,2) * t75;
t85 = mrSges(5,1) * t159 - mrSges(5,3) * t127;
t269 = t43 - t85;
t163 = pkin(2) * t232;
t111 = t205 * t245 + t163;
t57 = t183 * t111 + t180 * t113;
t266 = Ifges(5,4) * t127;
t262 = Ifges(5,2) * t126;
t261 = Ifges(4,6) * t184;
t260 = Ifges(5,6) * t159;
t259 = Ifges(6,3) * t123;
t254 = t170 * Ifges(4,5);
t253 = t180 * t31;
t135 = -mrSges(4,1) * t231 - mrSges(4,3) * t170;
t67 = -mrSges(5,1) * t126 + mrSges(5,2) * t127;
t252 = t67 - t135;
t249 = t179 * t181;
t248 = t180 * t185;
t247 = t181 * t182;
t246 = t311 * t170 + t312 * t232;
t146 = pkin(7) * t250 + t174;
t242 = qJD(2) * t184;
t238 = t178 * t275;
t237 = -qJD(2) + t170 / 0.2e1;
t235 = -0.3e1 / 0.2e1 * Ifges(4,6) - 0.3e1 / 0.2e1 * Ifges(3,4);
t130 = -t178 * qJ(3) - t146;
t229 = t177 * t242;
t227 = t185 * t239;
t104 = pkin(3) * t250 - t130;
t223 = t183 * t232;
t219 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t217 = pkin(4) * t183 + pkin(9) * t180;
t215 = t10 * t182 - t179 * t9;
t213 = mrSges(6,1) * t182 - mrSges(6,2) * t179;
t210 = Ifges(6,1) * t179 + t264;
t208 = Ifges(6,2) * t182 + t265;
t206 = Ifges(6,5) * t179 + Ifges(6,6) * t182;
t23 = mrSges(6,1) * t81 - mrSges(6,3) * t36;
t24 = -mrSges(6,2) * t81 + mrSges(6,3) * t37;
t204 = -t179 * t23 + t182 * t24;
t48 = pkin(9) * t251 + t268;
t55 = pkin(4) * t143 - pkin(9) * t144 + t104;
t18 = t179 * t55 + t182 * t48;
t17 = -t179 * t48 + t182 * t55;
t53 = -mrSges(6,2) * t123 + mrSges(6,3) * t74;
t54 = mrSges(6,1) * t123 - mrSges(6,3) * t75;
t203 = -t179 * t53 - t182 * t54;
t202 = t40 * t180 - t41 * t183;
t51 = -t105 * t180 + t183 * t91;
t56 = -t111 * t180 + t113 * t183;
t169 = qJD(2) * t238;
t141 = -pkin(7) * t230 + t169;
t101 = -t144 * t179 + t177 * t247;
t102 = t144 * t182 + t177 * t249;
t150 = pkin(4) * t180 - pkin(9) * t183 + qJ(3);
t125 = t179 * t150 + t182 * t248;
t124 = t182 * t150 - t179 * t248;
t21 = -t105 * t240 + t180 * t114 + t183 * t87 + t91 * t239;
t128 = -pkin(7) * t221 + t158;
t131 = (-pkin(2) * t184 + t225) * t177;
t175 = t178 * qJD(3);
t90 = t169 + t175 - t198;
t195 = (-qJ(3) * t242 - t241) * t177;
t116 = qJD(1) * t131;
t194 = Ifges(3,6) * t309 + (Ifges(3,4) * t181 + Ifges(3,2) * t184) * t308 + t254 / 0.2e1 + (-Ifges(4,6) * t181 - Ifges(4,3) * t184) * t245 / 0.2e1 - t116 * mrSges(4,2) - t140 * mrSges(3,3);
t142 = t146 * qJD(2);
t129 = qJD(1) * t142;
t96 = -t128 - t160;
t193 = -t128 * mrSges(3,2) - t96 * mrSges(4,3) + t311 * t129;
t190 = -t267 / 0.2e1 + t307;
t189 = t139 * mrSges(3,3) + t40 * mrSges(5,1) + Ifges(3,1) * t232 / 0.2e1 + Ifges(4,4) * t309 + (-t181 * Ifges(4,2) - t261) * t308 + t159 * Ifges(5,3) + t127 * Ifges(5,5) + t126 * Ifges(5,6) - t116 * mrSges(4,3) - t41 * mrSges(5,2) - t304;
t28 = t259 + t276 + t277;
t59 = t260 + t262 + t266;
t188 = t260 / 0.2e1 - t276 / 0.2e1 - t277 / 0.2e1 - t259 / 0.2e1 + t10 * mrSges(6,2) - t9 * mrSges(6,1) - t83 * mrSges(5,1) + t41 * mrSges(5,3) - t28 / 0.2e1 + t59 / 0.2e1 + t266 / 0.2e1;
t187 = (-t262 / 0.2e1 - t188) * t183;
t156 = Ifges(3,5) * t220;
t155 = Ifges(4,5) * t221;
t154 = Ifges(5,3) * t220;
t147 = qJD(4) * t217 + qJD(3);
t145 = -t171 + t238;
t138 = -qJ(3) * t231 + t163;
t137 = (mrSges(4,2) * t184 - mrSges(4,3) * t181) * t245;
t134 = -mrSges(3,2) * t170 + mrSges(3,3) * t231;
t132 = t178 * t233 + t171;
t122 = -t141 - t175;
t120 = (t179 * t184 + t180 * t247) * t245;
t119 = (-t180 * t249 + t182 * t184) * t245;
t118 = t170 * t179 - t182 * t223;
t117 = t170 * t182 + t179 * t223;
t115 = t165 + t195;
t112 = -qJD(1) * t224 + t168;
t110 = -t161 - t140;
t103 = -pkin(2) * t170 - t306;
t100 = -qJD(4) * t234 + t178 * t239 - t183 * t230;
t99 = -qJD(4) * t143 + t177 * t228;
t94 = qJD(1) * t195 + t157;
t84 = -mrSges(5,2) * t159 + mrSges(5,3) * t126;
t68 = pkin(4) * t127 - pkin(9) * t126;
t66 = -mrSges(5,2) * t220 - mrSges(5,3) * t81;
t64 = t168 + (-t217 - t289) * t232;
t62 = -qJD(5) * t125 + t182 * t147 - t179 * t227;
t61 = qJD(5) * t124 + t179 * t147 + t182 * t227;
t50 = pkin(9) * t231 + t57;
t49 = -pkin(4) * t231 - t56;
t47 = -pkin(4) * t251 - t51;
t46 = qJD(5) * t101 + t179 * t229 + t182 * t99;
t45 = -qJD(5) * t102 - t179 * t99 + t182 * t229;
t44 = mrSges(5,1) * t81 + mrSges(5,2) * t80;
t39 = Ifges(5,5) * t220 - t270 + t272;
t38 = -t81 * Ifges(5,2) + Ifges(5,6) * t220 + t271;
t35 = pkin(4) * t100 - pkin(9) * t99 + t90;
t26 = t179 * t64 + t182 * t50;
t25 = -t179 * t50 + t182 * t64;
t20 = t179 * t68 + t182 * t40;
t19 = -t179 * t40 + t182 * t68;
t16 = -pkin(4) * t229 - t22;
t15 = pkin(9) * t229 + t21;
t12 = -pkin(4) * t220 - t14;
t6 = Ifges(6,4) * t36 + Ifges(6,2) * t37 + Ifges(6,6) * t81;
t4 = -qJD(5) * t18 - t15 * t179 + t182 * t35;
t3 = qJD(5) * t17 + t15 * t182 + t179 * t35;
t52 = [t246 * t142 + (-t100 * t41 - t13 * t143 - t14 * t144 - t40 * t99) * mrSges(5,3) + t159 * (Ifges(5,5) * t99 - Ifges(5,6) * t100) / 0.2e1 + m(5) * (t104 * t73 + t13 * t268 + t14 * t51 + t21 * t41 + t22 * t40 + t83 * t90) + t268 * t66 + t126 * (Ifges(5,4) * t99 - Ifges(5,2) * t100) / 0.2e1 - t81 * (Ifges(5,4) * t144 - Ifges(5,2) * t143) / 0.2e1 + ((-t145 * mrSges(3,3) + t132 * mrSges(4,1) - t131 * mrSges(4,3) + Ifges(5,5) * t282 + Ifges(5,6) * t283 + (-Ifges(4,4) + t300) * t178 + (-t184 * t235 - 0.2e1 * t287) * t177) * t184 + (t130 * mrSges(4,1) - t131 * mrSges(4,2) - t146 * mrSges(3,3) + (-Ifges(3,6) + t299) * t178 + (t181 * t235 - 0.2e1 * t288) * t177 + (-0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + Ifges(5,3) / 0.2e1) * t250) * t181) * t226 + t127 * (Ifges(5,1) * t99 - Ifges(5,4) * t100) / 0.2e1 + t80 * (Ifges(5,1) * t144 - Ifges(5,4) * t143) / 0.2e1 + m(4) * (t103 * t142 + t110 * t122 + t115 * t116 + t129 * t132 + t130 * t96 + t131 * t94) + m(3) * (t128 * t146 - t129 * t145 + t139 * t142 + t140 * t141) + m(6) * (t1 * t18 + t10 * t3 + t12 * t47 + t16 * t31 + t17 * t2 + t4 * t9) + (t155 / 0.2e1 + t156 / 0.2e1 + t193) * t178 + ((-mrSges(4,1) * t96 + mrSges(4,2) * t94 + mrSges(3,3) * t128) * t184 + (-t94 * mrSges(4,3) + t154 / 0.2e1 + t312 * t129 + t305) * t181 + ((t110 * mrSges(4,1) + (t299 - Ifges(3,6) / 0.2e1) * t170 + t194) * t181 + (t103 * mrSges(4,1) + (-Ifges(4,4) / 0.2e1 + t300) * t170 + t189) * t184) * qJD(2)) * t177 + t39 * t282 + t38 * t283 + (Ifges(6,5) * t46 + Ifges(6,6) * t45 + Ifges(6,3) * t100) * t284 + (Ifges(6,5) * t102 + Ifges(6,6) * t101 + Ifges(6,3) * t143) * t290 + (Ifges(6,1) * t46 + Ifges(6,4) * t45 + Ifges(6,5) * t100) * t291 + t17 * t23 + t18 * t24 + t16 * t43 + t45 * t29 / 0.2e1 + t46 * t30 / 0.2e1 + t31 * (-mrSges(6,1) * t45 + mrSges(6,2) * t46) + t47 * t8 + t3 * t53 + t4 * t54 + t51 * t65 + t21 * t84 + t22 * t85 + t90 * t67 + t99 * t60 / 0.2e1 + t100 * t28 / 0.2e1 + t10 * (-mrSges(6,2) * t100 + mrSges(6,3) * t45) + t9 * (mrSges(6,1) * t100 - mrSges(6,3) * t46) - t100 * t59 / 0.2e1 + t83 * (mrSges(5,1) * t100 + mrSges(5,2) * t99) + t101 * t6 / 0.2e1 + t12 * (-mrSges(6,1) * t101 + mrSges(6,2) * t102) + t104 * t44 + t122 * t135 + t115 * t137 + t141 * t134 + t1 * (-mrSges(6,2) * t143 + mrSges(6,3) * t101) + t2 * (mrSges(6,1) * t143 - mrSges(6,3) * t102) + t73 * (mrSges(5,1) * t143 + mrSges(5,2) * t144) + (Ifges(6,4) * t46 + Ifges(6,2) * t45 + Ifges(6,6) * t100) * t293 + (Ifges(6,4) * t102 + Ifges(6,2) * t101 + Ifges(6,6) * t143) * t295 + (Ifges(6,1) * t102 + Ifges(6,4) * t101 + Ifges(6,5) * t143) * t296 + t102 * t301 + t143 * t302; (t187 + (t207 * t285 + t209 * t294 + t211 * t292 + t190 - t313) * t180 + (-m(5) * t202 + m(6) * t253 + t269 * t180 + t183 * t84) * t185) * qJD(4) + (-t26 + t61) * t53 + t252 * qJD(3) - t246 * t140 + m(5) * (t73 * qJ(3) + t83 * qJD(3)) + m(6) * (t1 * t125 + t10 * t61 + t2 * t124 + t9 * t62) + t156 + t155 + (-t25 + t62) * t54 + t193 + (t34 / 0.2e1 + t33 / 0.2e1 + t73 * mrSges(5,1) - t271 / 0.2e1 + t185 * t66 + t302 - t38 / 0.2e1 + (Ifges(6,3) / 0.2e1 + t298) * t81 + (m(5) * t185 - mrSges(5,3)) * t13 + t219) * t180 + (-t135 + t134) * t139 + (t73 * mrSges(5,2) - t14 * mrSges(5,3) + t272 / 0.2e1 - t270 / 0.2e1 + t39 / 0.2e1 + t12 * t212 + t211 * t296 + t209 * t295 + t207 * t290 + t6 * t281 + t7 * t280 + (-t1 * t179 - t182 * t2) * mrSges(6,3) + (m(5) * t14 - m(6) * t12 + t279) * t185 + (-mrSges(6,3) * t215 + t182 * t297 + t206 * t285 + t208 * t294 + t210 * t292 + t213 * t31 + t281 * t30) * qJD(5)) * t183 - m(5) * (t112 * t83 + t40 * t56 + t41 * t57) - m(6) * (t10 * t26 + t25 * t9 + t31 * t49) + (-pkin(2) * t129 - qJ(3) * t96 - t103 * t140 + t110 * t306 - t116 * t138) * m(4) + (((t287 - t261 / 0.2e1) * t245 + qJD(2) * (Ifges(5,5) * t183 - Ifges(5,6) * t180) / 0.2e1 - t189 + (-pkin(2) * qJD(2) - t103) * mrSges(4,1) + t237 * Ifges(4,4) + t304) * t184 + (-t254 / 0.2e1 + ((t288 + (Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1) * t181) * t177 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 - Ifges(4,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t250) * qJD(1) + t237 * Ifges(3,6) + (-qJ(3) * qJD(2) - t110) * mrSges(4,1) + t190 * t180 + t187 - t194) * t181) * t245 + (Ifges(6,5) * t120 + Ifges(6,6) * t119) * t285 + (Ifges(6,1) * t120 + Ifges(6,4) * t119) * t292 + qJ(3) * t44 - t49 * t43 - t57 * t84 - t56 * t85 + (-t10 * t119 + t9 * t120) * mrSges(6,3) - t112 * t67 - t120 * t30 / 0.2e1 - t31 * (-mrSges(6,1) * t119 + mrSges(6,2) * t120) + t124 * t23 + t125 * t24 - t138 * t137 + (Ifges(6,4) * t120 + Ifges(6,2) * t119) * t294 + t119 * t297; -t117 * t54 - t118 * t53 - t252 * t170 + (mrSges(4,1) * t242 + t137 * t181) * t245 + (t84 * t232 + (-t179 * t54 + t182 * t53 + t84) * qJD(4) + t279) * t183 + (t203 * qJD(5) + t159 * t269 + t204 + t66) * t180 + (-t10 * t118 - t117 * t9 + t232 * t253 + (qJD(4) * t215 - t12) * t183 + (qJD(4) * t31 - qJD(5) * t214 + t216) * t180) * m(6) + (t13 * t180 + t14 * t183 - t159 * t202 - t170 * t83) * m(5) + (t110 * t170 + t116 * t232 + t129) * m(4); (-pkin(4) * t12 - t10 * t20 - t19 * t9 - t31 * t41) * m(6) - t269 * t41 + t305 + t310 * qJD(5) + (m(6) * t216 + t204 + (-m(6) * t214 + t203) * qJD(5)) * pkin(9) + t154 + ((t298 - Ifges(5,1) / 0.2e1) * t127 + t307 - t310) * t126 + t188 * t127 + t6 * t280 + t206 * t290 - pkin(4) * t8 - t20 * t53 - t19 * t54 - t40 * t84 - t12 * t213 + t216 * mrSges(6,3) + t208 * t295 + t210 * t296 + t179 * t301; -t31 * (mrSges(6,1) * t75 + mrSges(6,2) * t74) + (Ifges(6,1) * t74 - t278) * t292 + t29 * t291 + (Ifges(6,5) * t74 - Ifges(6,6) * t75) * t285 - t9 * t53 + t10 * t54 + (t10 * t75 + t74 * t9) * mrSges(6,3) + t219 + t5 + (-Ifges(6,2) * t75 + t30 + t71) * t294;];
tauc = t52(:);
