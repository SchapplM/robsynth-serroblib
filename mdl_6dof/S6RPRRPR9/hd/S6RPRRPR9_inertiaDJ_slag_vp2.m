% Calculate time derivative of joint inertia matrix for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:27:02
% EndTime: 2019-03-09 05:27:19
% DurationCPUTime: 7.21s
% Computational Cost: add. (17228->599), mult. (49872->901), div. (0->0), fcn. (53576->14), ass. (0->259)
t202 = sin(pkin(7));
t205 = cos(pkin(7));
t206 = cos(pkin(6));
t201 = sin(pkin(12));
t203 = sin(pkin(6));
t204 = cos(pkin(12));
t263 = t203 * t204;
t298 = pkin(1) * t206;
t255 = qJ(2) * t263 + t201 * t298;
t140 = (t202 * t206 + t205 * t263) * pkin(9) + t255;
t209 = sin(qJ(3));
t212 = cos(qJ(3));
t192 = t204 * t298;
t267 = t201 * t203;
t146 = pkin(2) * t206 + t192 + (-pkin(9) * t205 - qJ(2)) * t267;
t156 = (-pkin(9) * t201 * t202 - pkin(2) * t204 - pkin(1)) * t203;
t220 = t146 * t205 + t156 * t202;
t87 = -t209 * t140 + t220 * t212;
t211 = cos(qJ(4));
t200 = sin(pkin(13));
t208 = sin(qJ(4));
t268 = t200 * t208;
t275 = cos(pkin(13));
t216 = t211 * t275 - t268;
t167 = t216 * qJD(4);
t228 = t275 * t208;
t171 = t200 * t211 + t228;
t210 = cos(qJ(6));
t207 = sin(qJ(6));
t249 = qJD(6) * t207;
t217 = -t210 * t167 + t171 * t249;
t261 = t205 * t212;
t264 = t202 * t212;
t326 = t203 * (-t201 * t209 + t204 * t261) + t206 * t264;
t325 = -2 * Ifges(4,4);
t301 = t207 / 0.2e1;
t300 = t210 / 0.2e1;
t231 = -t249 / 0.2e1;
t297 = pkin(4) * t200;
t194 = pkin(11) + t297;
t324 = m(7) * t194;
t166 = t171 * qJD(4);
t250 = qJD(4) * t211;
t322 = Ifges(5,5) * t250 + Ifges(6,5) * t167 - Ifges(6,6) * t166;
t265 = t202 * t209;
t168 = t205 * t211 - t208 * t265;
t252 = qJD(3) * t212;
t236 = t202 * t252;
t151 = qJD(4) * t168 + t211 * t236;
t169 = t205 * t208 + t211 * t265;
t152 = -qJD(4) * t169 - t208 * t236;
t107 = t151 * t275 + t200 * t152;
t125 = t200 * t168 + t169 * t275;
t114 = -t207 * t125 - t210 * t264;
t253 = qJD(3) * t209;
t237 = t202 * t253;
t73 = qJD(6) * t114 + t210 * t107 + t207 * t237;
t219 = -t210 * t125 + t207 * t264;
t74 = qJD(6) * t219 - t207 * t107 + t210 * t237;
t321 = -t74 * t207 + t73 * t210;
t196 = -pkin(4) * t211 - pkin(3);
t139 = -pkin(5) * t216 - pkin(11) * t171 + t196;
t294 = -qJ(5) - pkin(10);
t183 = t294 * t211;
t154 = -t183 * t275 + t268 * t294;
t104 = t139 * t210 - t154 * t207;
t229 = qJD(4) * t294;
t164 = qJD(5) * t211 + t208 * t229;
t214 = -qJD(5) * t208 + t211 * t229;
t122 = t164 * t275 + t200 * t214;
t251 = qJD(4) * t208;
t244 = pkin(4) * t251;
t123 = pkin(5) * t166 - pkin(11) * t167 + t244;
t62 = qJD(6) * t104 + t122 * t210 + t123 * t207;
t105 = t139 * t207 + t154 * t210;
t63 = -qJD(6) * t105 - t122 * t207 + t123 * t210;
t320 = -t207 * t63 + t210 * t62;
t128 = t212 * t140;
t254 = qJD(2) * t203;
t76 = (t201 * t261 + t204 * t209) * t254 + (t209 * t220 + t128) * qJD(3);
t262 = t205 * t209;
t145 = t206 * t265 + (t201 * t212 + t204 * t262) * t203;
t165 = -t202 * t263 + t205 * t206;
t113 = t145 * t211 + t165 * t208;
t137 = t326 * qJD(3);
t95 = -qJD(4) * t113 - t137 * t208;
t53 = -t95 * pkin(4) + t76;
t112 = -t145 * t208 + t165 * t211;
t96 = qJD(4) * t112 + t137 * t211;
t58 = t200 * t96 - t275 * t95;
t59 = t200 * t95 + t275 * t96;
t21 = t58 * pkin(5) - t59 * pkin(11) + t53;
t138 = t145 * qJD(3);
t238 = t201 * t254;
t227 = t202 * t238;
t101 = pkin(3) * t138 - pkin(10) * t137 + t227;
t108 = -t146 * t202 + t205 * t156;
t80 = -pkin(3) * t326 - pkin(10) * t145 + t108;
t88 = t146 * t262 + t156 * t265 + t128;
t83 = pkin(10) * t165 + t88;
t47 = t208 * t80 + t211 * t83;
t75 = (-t201 * t262 + t204 * t212) * t254 + t87 * qJD(3);
t25 = -qJD(4) * t47 + t211 * t101 - t208 * t75;
t16 = pkin(4) * t138 - qJ(5) * t96 - qJD(5) * t113 + t25;
t24 = t208 * t101 + t211 * t75 + t80 * t250 - t251 * t83;
t20 = qJ(5) * t95 + qJD(5) * t112 + t24;
t6 = t200 * t16 + t275 * t20;
t4 = pkin(11) * t138 + t6;
t46 = -t208 * t83 + t211 * t80;
t38 = -pkin(4) * t326 - qJ(5) * t113 + t46;
t41 = qJ(5) * t112 + t47;
t18 = t200 * t38 + t275 * t41;
t14 = -pkin(11) * t326 + t18;
t82 = -t165 * pkin(3) - t87;
t61 = -t112 * pkin(4) + t82;
t84 = -t112 * t275 + t113 * t200;
t85 = t200 * t112 + t113 * t275;
t31 = t84 * pkin(5) - t85 * pkin(11) + t61;
t7 = -t14 * t207 + t210 * t31;
t1 = qJD(6) * t7 + t207 * t21 + t210 * t4;
t8 = t14 * t210 + t207 * t31;
t2 = -qJD(6) * t8 - t207 * t4 + t21 * t210;
t319 = t1 * t210 - t2 * t207;
t310 = m(6) * pkin(4);
t318 = t200 * t310 - mrSges(6,2);
t181 = -mrSges(7,1) * t210 + mrSges(7,2) * t207;
t239 = t275 * pkin(4);
t195 = -t239 - pkin(5);
t317 = m(7) * t195 - t275 * t310 - mrSges(6,1) + t181;
t316 = 0.2e1 * m(6);
t315 = 0.2e1 * m(7);
t314 = -2 * mrSges(6,3);
t121 = t164 * t200 - t214 * t275;
t313 = 0.2e1 * t121;
t153 = -t183 * t200 - t228 * t294;
t312 = 0.2e1 * t153;
t311 = m(5) * pkin(3);
t64 = -t207 * t85 - t210 * t326;
t35 = qJD(6) * t64 + t138 * t207 + t210 * t59;
t309 = t35 / 0.2e1;
t65 = -t207 * t326 + t210 * t85;
t36 = -qJD(6) * t65 + t138 * t210 - t207 * t59;
t308 = t36 / 0.2e1;
t307 = t64 / 0.2e1;
t248 = qJD(6) * t210;
t197 = Ifges(7,5) * t248;
t305 = Ifges(7,6) * t231 + t197 / 0.2e1;
t289 = Ifges(7,4) * t207;
t224 = Ifges(7,1) * t210 - t289;
t178 = t224 * qJD(6);
t304 = t178 / 0.2e1;
t303 = Ifges(7,5) * t301 + Ifges(7,6) * t300;
t302 = -t207 / 0.2e1;
t299 = t211 / 0.2e1;
t12 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t50 = mrSges(6,1) * t138 - mrSges(6,3) * t59;
t293 = t12 - t50;
t39 = -mrSges(7,1) * t64 + mrSges(7,2) * t65;
t67 = -mrSges(6,1) * t326 - mrSges(6,3) * t85;
t292 = t39 - t67;
t291 = Ifges(5,4) * t208;
t290 = Ifges(5,4) * t211;
t288 = Ifges(7,4) * t210;
t287 = Ifges(7,6) * t207;
t286 = t137 * mrSges(4,3);
t285 = t138 * mrSges(4,3);
t284 = t138 * Ifges(6,5);
t283 = t138 * Ifges(6,6);
t282 = t326 * Ifges(5,6);
t279 = t212 * t76;
t276 = -mrSges(5,1) * t211 + mrSges(5,2) * t208 - mrSges(4,1);
t274 = t121 * t153;
t106 = t151 * t200 - t152 * t275;
t124 = -t168 * t275 + t169 * t200;
t273 = t124 * t106;
t272 = t171 * t207;
t271 = t171 * t210;
t270 = t194 * t207;
t269 = t194 * t210;
t185 = Ifges(7,2) * t210 + t289;
t260 = t207 * t185;
t187 = Ifges(7,1) * t207 + t288;
t258 = t210 * t187;
t257 = Ifges(4,5) * t137 - Ifges(4,6) * t138;
t9 = Ifges(7,5) * t35 + Ifges(7,6) * t36 + Ifges(7,3) * t58;
t246 = Ifges(6,5) * t59 - Ifges(6,6) * t58 + Ifges(6,3) * t138;
t245 = Ifges(5,5) * t96 + Ifges(5,6) * t95 + Ifges(5,3) * t138;
t243 = mrSges(7,3) * t249;
t242 = mrSges(7,3) * t248;
t235 = t194 * t249;
t234 = t194 * t248;
t232 = t171 * t248;
t32 = t58 * mrSges(6,1) + t59 * mrSges(6,2);
t230 = t248 / 0.2e1;
t133 = t166 * mrSges(6,1) + t167 * mrSges(6,2);
t226 = t202 ^ 2 * t209 * t252;
t225 = mrSges(7,1) * t207 + mrSges(7,2) * t210;
t223 = -Ifges(7,2) * t207 + t288;
t222 = -t208 * t25 + t211 * t24;
t221 = t153 * t106 + t121 * t124;
t5 = t16 * t275 - t200 * t20;
t17 = -t200 * t41 + t275 * t38;
t218 = t207 * t167 + t232;
t92 = -Ifges(7,5) * t217 - Ifges(7,6) * t218 + Ifges(7,3) * t166;
t188 = Ifges(5,1) * t208 + t290;
t186 = Ifges(5,2) * t211 + t291;
t179 = (Ifges(5,1) * t211 - t291) * qJD(4);
t177 = (-Ifges(5,2) * t208 + t290) * qJD(4);
t176 = t223 * qJD(6);
t174 = (mrSges(5,1) * t208 + mrSges(5,2) * t211) * qJD(4);
t173 = t225 * qJD(6);
t149 = Ifges(6,1) * t171 + Ifges(6,4) * t216;
t148 = Ifges(6,4) * t171 + Ifges(6,2) * t216;
t147 = -mrSges(6,1) * t216 + mrSges(6,2) * t171;
t142 = -mrSges(7,1) * t216 - mrSges(7,3) * t271;
t141 = mrSges(7,2) * t216 - mrSges(7,3) * t272;
t136 = t225 * t171;
t135 = Ifges(6,1) * t167 - Ifges(6,4) * t166;
t134 = Ifges(6,4) * t167 - Ifges(6,2) * t166;
t120 = -Ifges(7,5) * t216 + t171 * t224;
t119 = -Ifges(7,6) * t216 + t171 * t223;
t118 = -Ifges(7,3) * t216 + (Ifges(7,5) * t210 - t287) * t171;
t117 = mrSges(4,1) * t165 - mrSges(4,3) * t145;
t116 = -mrSges(4,2) * t165 + mrSges(4,3) * t326;
t111 = -mrSges(7,2) * t166 - mrSges(7,3) * t218;
t110 = mrSges(7,1) * t166 + mrSges(7,3) * t217;
t103 = mrSges(4,1) * t138 + mrSges(4,2) * t137;
t102 = mrSges(7,1) * t218 - mrSges(7,2) * t217;
t98 = -mrSges(5,1) * t326 - mrSges(5,3) * t113;
t97 = mrSges(5,2) * t326 + mrSges(5,3) * t112;
t94 = -Ifges(7,1) * t217 - Ifges(7,4) * t218 + Ifges(7,5) * t166;
t93 = -Ifges(7,4) * t217 - Ifges(7,2) * t218 + Ifges(7,6) * t166;
t86 = -mrSges(5,1) * t112 + mrSges(5,2) * t113;
t71 = mrSges(5,1) * t138 - mrSges(5,3) * t96;
t70 = -mrSges(5,2) * t138 + mrSges(5,3) * t95;
t69 = Ifges(5,1) * t113 + Ifges(5,4) * t112 - Ifges(5,5) * t326;
t68 = Ifges(5,4) * t113 + Ifges(5,2) * t112 - t282;
t66 = mrSges(6,2) * t326 - mrSges(6,3) * t84;
t60 = -mrSges(5,1) * t95 + mrSges(5,2) * t96;
t52 = Ifges(5,1) * t96 + Ifges(5,4) * t95 + t138 * Ifges(5,5);
t51 = Ifges(5,4) * t96 + Ifges(5,2) * t95 + t138 * Ifges(5,6);
t49 = -mrSges(6,2) * t138 - mrSges(6,3) * t58;
t48 = mrSges(6,1) * t84 + mrSges(6,2) * t85;
t45 = Ifges(6,1) * t85 - Ifges(6,4) * t84 - Ifges(6,5) * t326;
t44 = Ifges(6,4) * t85 - Ifges(6,2) * t84 - Ifges(6,6) * t326;
t43 = mrSges(7,1) * t84 - mrSges(7,3) * t65;
t42 = -mrSges(7,2) * t84 + mrSges(7,3) * t64;
t30 = Ifges(7,1) * t65 + Ifges(7,4) * t64 + Ifges(7,5) * t84;
t29 = Ifges(7,4) * t65 + Ifges(7,2) * t64 + Ifges(7,6) * t84;
t28 = Ifges(7,5) * t65 + Ifges(7,6) * t64 + Ifges(7,3) * t84;
t27 = Ifges(6,1) * t59 - Ifges(6,4) * t58 + t284;
t26 = Ifges(6,4) * t59 - Ifges(6,2) * t58 + t283;
t23 = -mrSges(7,2) * t58 + mrSges(7,3) * t36;
t22 = mrSges(7,1) * t58 - mrSges(7,3) * t35;
t13 = pkin(5) * t326 - t17;
t11 = Ifges(7,1) * t35 + Ifges(7,4) * t36 + Ifges(7,5) * t58;
t10 = Ifges(7,4) * t35 + Ifges(7,2) * t36 + Ifges(7,6) * t58;
t3 = -t138 * pkin(5) - t5;
t15 = [0.2e1 * m(5) * (t24 * t47 + t25 * t46 + t76 * t82) + (-t44 + t28) * t58 + (t9 - t26) * t84 + (t1 * t8 + t13 * t3 + t2 * t7) * t315 + (t17 * t5 + t18 * t6 + t53 * t61) * t316 + (t145 * t325 + Ifges(5,5) * t113 + Ifges(6,5) * t85 - Ifges(4,6) * t165 + Ifges(5,6) * t112 - Ifges(6,6) * t84 - ((2 * Ifges(4,2)) + Ifges(5,3) + Ifges(6,3)) * t326) * t138 + (0.2e1 * Ifges(4,1) * t145 + Ifges(4,5) * t165 - t325 * t326) * t137 + 0.2e1 * (-mrSges(4,1) * t326 + mrSges(4,2) * t145) * t227 - t326 * t245 - t326 * t246 + 0.2e1 * (t204 * (-mrSges(3,2) * t206 + mrSges(3,3) * t263) + m(3) * (t255 * t204 + (qJ(2) * t267 - t192) * t201)) * t254 + t113 * t52 + 0.2e1 * t75 * t116 - 0.2e1 * t76 * t117 + 0.2e1 * t108 * t103 + t112 * t51 + t95 * t68 + t96 * t69 + 0.2e1 * t24 * t97 + 0.2e1 * t25 * t98 + t85 * t27 + 0.2e1 * t76 * t86 + 0.2e1 * t82 * t60 + 0.2e1 * t47 * t70 + 0.2e1 * t46 * t71 + t64 * t10 + t65 * t11 + 0.2e1 * t6 * t66 + 0.2e1 * t5 * t67 + t59 * t45 + 0.2e1 * t61 * t32 + 0.2e1 * t18 * t49 + 0.2e1 * t17 * t50 + 0.2e1 * t53 * t48 + 0.2e1 * t1 * t42 + 0.2e1 * t2 * t43 + t36 * t29 + 0.2e1 * t3 * t39 + t35 * t30 + 0.2e1 * t7 * t22 + 0.2e1 * t8 * t23 + 0.2e1 * t13 * t12 + 0.2e1 * m(4) * (t108 * t227 + t75 * t88 - t76 * t87) + t165 * t257 - 0.2e1 * (mrSges(3,1) * t206 - mrSges(3,3) * t267) * t238 - 0.2e1 * t88 * t285 - 0.2e1 * t87 * t286; t205 * t103 + t107 * t66 + t114 * t22 - t219 * t23 + t125 * t49 + t151 * t97 + t152 * t98 + t168 * t71 + t169 * t70 + t73 * t42 + t74 * t43 + t293 * t124 + t292 * t106 + m(5) * (t151 * t47 + t152 * t46 + t168 * t25 + t169 * t24) + m(6) * (-t106 * t17 + t107 * t18 - t124 * t5 + t125 * t6) + m(7) * (-t1 * t219 + t106 * t13 + t114 * t2 + t124 * t3 + t7 * t74 + t73 * t8) + (-t209 * t285 + (-t32 - t60 - t286) * t212 + (t212 * t116 + (-t117 + t48 + t86) * t209) * qJD(3) + m(4) * (t205 * t238 + t209 * t75 + t252 * t88 - t253 * t87 - t279) + m(5) * (t253 * t82 - t279) + m(6) * (-t212 * t53 + t253 * t61)) * t202; 0.2e1 * m(6) * (t125 * t107 - t226 + t273) + 0.2e1 * m(7) * (t114 * t74 - t219 * t73 + t273) + 0.2e1 * m(5) * (t169 * t151 + t168 * t152 - t226); t257 + (t69 * t299 + (pkin(4) * t48 + t282 / 0.2e1 - t68 / 0.2e1) * t208) * qJD(4) + ((-t208 * t47 - t211 * t46) * qJD(4) + t222) * mrSges(5,3) + t93 * t307 + t119 * t308 + t120 * t309 - t322 * t326 / 0.2e1 + (t28 / 0.2e1 - t44 / 0.2e1 - t18 * mrSges(6,3)) * t166 + t138 * (Ifges(5,5) * t208 + Ifges(5,6) * t211) / 0.2e1 + t208 * t52 / 0.2e1 + t95 * t186 / 0.2e1 + t96 * t188 / 0.2e1 + t196 * t32 + t82 * t174 + t112 * t177 / 0.2e1 + t113 * t179 / 0.2e1 + t154 * t49 + t53 * t147 + t59 * t149 / 0.2e1 + t1 * t141 + t2 * t142 + t61 * t133 + t85 * t135 / 0.2e1 + t3 * t136 + t122 * t66 + t7 * t110 + t8 * t111 + t13 * t102 + t104 * t22 + t105 * t23 + t65 * t94 / 0.2e1 - t75 * mrSges(4,2) + t63 * t43 - pkin(3) * t60 + t62 * t42 + m(7) * (t1 * t105 + t104 * t2 + t121 * t13 + t153 * t3 + t62 * t8 + t63 * t7) + (-t134 / 0.2e1 + t92 / 0.2e1) * t84 + (-t148 / 0.2e1 + t118 / 0.2e1) * t58 + t51 * t299 + m(6) * (-t121 * t17 + t122 * t18 - t153 * t5 + t154 * t6 + t196 * t53 + t244 * t61) + (-t97 * t251 + m(5) * (-t250 * t46 - t251 * t47 + t222) + t211 * t70 - t208 * t71 - t98 * t250) * pkin(10) + (t276 - t311) * t76 + t292 * t121 + t293 * t153 + (t45 / 0.2e1 + t30 * t300 + t29 * t302 - t17 * mrSges(6,3)) * t167 + (t284 / 0.2e1 + t27 / 0.2e1 + t11 * t300 + t10 * t302 - t5 * mrSges(6,3) + (t30 * t302 - t210 * t29 / 0.2e1) * qJD(6)) * t171 - (t9 / 0.2e1 - t26 / 0.2e1 - t283 / 0.2e1 - t6 * mrSges(6,3)) * t216; t124 * t102 + t106 * t136 + t114 * t110 - t219 * t111 + t73 * t141 + t74 * t142 + ((-t133 - t174) * t212 + (-mrSges(4,2) * t212 + (t147 + t276) * t209) * qJD(3)) * t202 + m(6) * (t154 * t107 + t122 * t125 + (t196 * t253 - t212 * t244) * t202 + t221) + m(7) * (t104 * t74 + t105 * t73 + t114 * t63 - t219 * t62 + t221) - t237 * t311 + (t106 * t171 + t107 * t216 + t124 * t167 - t125 * t166) * mrSges(6,3) + (m(5) * pkin(10) + mrSges(5,3)) * (t151 * t211 - t152 * t208 + (-t168 * t211 - t169 * t208) * qJD(4)); -0.2e1 * pkin(3) * t174 + t102 * t312 + 0.2e1 * t104 * t110 + 0.2e1 * t105 * t111 + t136 * t313 + 0.2e1 * t196 * t133 + 0.2e1 * t62 * t141 + 0.2e1 * t63 * t142 + t211 * t177 + t208 * t179 + (t211 * t188 + (0.2e1 * pkin(4) * t147 - t186) * t208) * qJD(4) + (t122 * t154 + t196 * t244 + t274) * t316 + (t104 * t63 + t105 * t62 + t274) * t315 - (t122 * t314 - t134 + t92) * t216 + (t154 * t314 + t118 - t148) * t166 + (mrSges(6,3) * t312 - t207 * t119 + t210 * t120 + t149) * t167 + (mrSges(6,3) * t313 - t207 * t93 + t210 * t94 + t135 + (-t119 * t210 - t120 * t207) * qJD(6)) * t171; t30 * t230 + t29 * t231 + t58 * t303 + t65 * t304 + t84 * t305 + t176 * t307 + t185 * t308 + t187 * t309 + (t200 * t6 + t275 * t5) * t310 + t245 + t246 + t195 * t12 + t13 * t173 + t3 * t181 - t24 * mrSges(5,2) + t25 * mrSges(5,1) - t6 * mrSges(6,2) + t5 * mrSges(6,1) + m(7) * (t195 * t3 + ((-t207 * t8 - t210 * t7) * qJD(6) + t319) * t194) + t319 * mrSges(7,3) + t50 * t239 + t23 * t269 + t49 * t297 + t10 * t300 + t11 * t301 - t43 * t234 - t42 * t235 - t7 * t242 - t8 * t243 - t22 * t270; ((-t114 * t210 + t207 * t219) * qJD(6) + t321) * t324 + t219 * t243 - t114 * t242 + t124 * t173 - t151 * mrSges(5,2) + t152 * mrSges(5,1) + t318 * t107 + t317 * t106 + t321 * mrSges(7,3); t322 + t120 * t230 + (t171 * t187 + t119) * t231 + (-mrSges(6,3) * t239 + t258 / 0.2e1 - t260 / 0.2e1) * t167 + (-mrSges(6,3) * t297 + t303) * t166 + ((-t104 * t210 - t105 * t207) * qJD(6) + t320) * t324 + (-mrSges(5,1) * t250 + mrSges(5,2) * t251) * pkin(10) + t271 * t304 + t195 * t102 + t153 * t173 + t320 * mrSges(7,3) + t318 * t122 + t317 * t121 + t111 * t269 + t93 * t300 + t94 * t301 - t185 * t232 / 0.2e1 - t142 * t234 - t141 * t235 - t104 * t242 - t105 * t243 - Ifges(5,6) * t251 - t110 * t270 - t176 * t272 / 0.2e1 - t216 * t305; 0.2e1 * t173 * t195 + t176 * t210 + t178 * t207 + (t258 - t260) * qJD(6); t207 * t23 + t210 * t22 + (-t207 * t43 + t210 * t42) * qJD(6) + m(7) * (t1 * t207 + t2 * t210 + (-t207 * t7 + t210 * t8) * qJD(6)) + m(6) * t53 + t32; m(6) * t237 + m(7) * (t207 * t73 + t210 * t74 + (-t114 * t207 - t210 * t219) * qJD(6)); m(7) * (t207 * t62 + t210 * t63 + (-t104 * t207 + t105 * t210) * qJD(6)) + t141 * t248 + t207 * t111 - t142 * t249 + t210 * t110 + m(6) * t244 + t133; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t9; mrSges(7,1) * t74 - mrSges(7,2) * t73; mrSges(7,1) * t63 - mrSges(7,2) * t62 + t92; t197 + (t181 * t194 - t287) * qJD(6); -t173; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
