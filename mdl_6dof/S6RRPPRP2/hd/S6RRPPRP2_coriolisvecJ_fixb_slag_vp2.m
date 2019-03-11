% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:29:41
% EndTime: 2019-03-09 08:30:00
% DurationCPUTime: 10.11s
% Computational Cost: add. (5318->534), mult. (13635->665), div. (0->0), fcn. (9232->6), ass. (0->256)
t339 = Ifges(6,4) + Ifges(7,4);
t340 = Ifges(6,1) + Ifges(7,1);
t331 = Ifges(6,5) + Ifges(7,5);
t338 = Ifges(6,2) + Ifges(7,2);
t337 = Ifges(6,6) + Ifges(7,6);
t324 = mrSges(5,2) - mrSges(4,1);
t175 = sin(pkin(9));
t179 = cos(qJ(2));
t250 = cos(pkin(9));
t218 = t250 * t179;
t215 = qJD(1) * t218;
t177 = sin(qJ(2));
t241 = qJD(1) * t177;
t137 = t175 * t241 - t215;
t176 = sin(qJ(5));
t178 = cos(qJ(5));
t112 = -qJD(2) * t176 + t137 * t178;
t336 = t339 * t112;
t113 = qJD(2) * t178 + t137 * t176;
t335 = t339 * t113;
t334 = t339 * t178;
t333 = t339 * t176;
t332 = mrSges(4,3) + mrSges(5,1);
t150 = t175 * t179 + t177 * t250;
t139 = t150 * qJD(1);
t134 = qJD(5) + t139;
t319 = t112 * t338 + t134 * t337 + t335;
t318 = t113 * t340 + t331 * t134 + t336;
t208 = mrSges(7,1) * t178 - mrSges(7,2) * t176;
t210 = mrSges(6,1) * t178 - mrSges(6,2) * t176;
t270 = -qJ(3) - pkin(7);
t158 = t270 * t177;
t152 = qJD(1) * t158;
t148 = qJD(2) * pkin(2) + t152;
t159 = t270 * t179;
t153 = qJD(1) * t159;
t219 = t250 * t153;
t105 = t175 * t148 - t219;
t100 = -qJD(2) * qJ(4) - t105;
t281 = pkin(4) * t137;
t63 = -t100 - t281;
t32 = -pkin(5) * t112 + qJD(6) + t63;
t329 = t32 * t208 + t63 * t210;
t328 = -t176 * t337 + t178 * t331;
t327 = -t176 * t338 + t334;
t326 = t178 * t340 - t333;
t141 = t175 * t153;
t107 = t152 * t250 + t141;
t325 = -qJD(4) + t107;
t323 = Ifges(5,4) - Ifges(4,5);
t322 = Ifges(5,5) - Ifges(4,6);
t235 = qJD(1) * qJD(2);
t225 = t177 * t235;
t130 = qJD(2) * t215 - t175 * t225;
t138 = t150 * qJD(2);
t129 = qJD(1) * t138;
t80 = qJD(5) * t112 + t129 * t176;
t81 = -qJD(5) * t113 + t129 * t178;
t321 = t130 * t337 + t338 * t81 + t339 * t80;
t320 = t130 * t331 + t339 * t81 + t340 * t80;
t236 = qJD(6) * t178;
t238 = qJD(5) * t176;
t228 = t250 * pkin(2);
t168 = -t228 - pkin(3);
t164 = -pkin(8) + t168;
t244 = qJ(6) - t164;
t171 = pkin(2) * t241;
t216 = qJ(4) * t137 + t171;
t297 = pkin(3) + pkin(8);
t58 = t139 * t297 + t216;
t106 = t152 * t175 - t219;
t78 = t106 - t281;
t69 = t178 * t78;
t317 = pkin(5) * t137 - t69 - (-qJ(6) * t139 - t58) * t176 + t238 * t244 - t236;
t145 = t244 * t178;
t24 = t176 * t78 + t178 * t58;
t245 = t139 * t178;
t316 = -qJ(6) * t245 - qJD(5) * t145 - qJD(6) * t176 - t24;
t229 = -pkin(5) * t178 - pkin(4);
t237 = qJD(5) * t178;
t315 = pkin(5) * t237 - t139 * t229 - t325;
t265 = mrSges(5,1) * t137;
t121 = -qJD(2) * mrSges(5,3) + t265;
t67 = -mrSges(6,1) * t112 + mrSges(6,2) * t113;
t314 = t67 - t121;
t255 = t137 * mrSges(4,3);
t119 = -qJD(2) * mrSges(4,2) - t255;
t313 = t119 - t121;
t253 = t139 * mrSges(4,3);
t254 = t139 * mrSges(5,1);
t312 = -qJD(2) * t324 - t253 - t254;
t311 = t176 * t331 + t178 * t337;
t310 = t178 * t338 + t333;
t309 = t176 * t340 + t334;
t308 = -t237 - t245;
t246 = t139 * t176;
t307 = t238 + t246;
t165 = pkin(2) * t225;
t186 = -qJ(4) * t130 - qJD(4) * t139 + t165;
t31 = t129 * t297 + t186;
t220 = qJD(2) * t270;
t135 = qJD(3) * t179 + t177 * t220;
t123 = t135 * qJD(1);
t136 = -t177 * qJD(3) + t179 * t220;
t181 = qJD(1) * t136;
t76 = t123 * t175 - t250 * t181;
t48 = pkin(4) * t130 + t76;
t169 = -pkin(2) * t179 - pkin(1);
t242 = qJD(1) * t169;
t155 = qJD(3) + t242;
t182 = -t139 * qJ(4) + t155;
t52 = t137 * t297 + t182;
t104 = t148 * t250 + t141;
t190 = qJD(4) - t104;
t277 = t139 * pkin(4);
t57 = -qJD(2) * t297 + t190 + t277;
t5 = t176 * t48 + t178 * t31 + t57 * t237 - t238 * t52;
t16 = t176 * t57 + t178 * t52;
t6 = -qJD(5) * t16 - t176 * t31 + t178 * t48;
t306 = -t176 * t5 - t178 * t6;
t1 = pkin(5) * t130 - qJ(6) * t80 - qJD(6) * t113 + t6;
t2 = qJ(6) * t81 + qJD(6) * t112 + t5;
t305 = -t1 * t178 - t176 * t2;
t240 = qJD(1) * t179;
t248 = Ifges(3,6) * qJD(2);
t264 = Ifges(3,4) * t177;
t304 = t248 / 0.2e1 + (t179 * Ifges(3,2) + t264) * qJD(1) / 0.2e1 + pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t240);
t87 = t137 * pkin(3) + t182;
t303 = -t155 * mrSges(4,1) + t87 * mrSges(5,2) + t329;
t11 = qJ(6) * t112 + t16;
t15 = -t176 * t52 + t178 * t57;
t10 = -qJ(6) * t113 + t15;
t9 = pkin(5) * t134 + t10;
t302 = t15 * mrSges(6,1) + t9 * mrSges(7,1) + t155 * mrSges(4,2) - t16 * mrSges(6,2) - t11 * mrSges(7,2) - t87 * mrSges(5,3);
t301 = -0.2e1 * pkin(1);
t299 = t80 / 0.2e1;
t298 = t81 / 0.2e1;
t296 = -t112 / 0.2e1;
t295 = t112 / 0.2e1;
t294 = -t113 / 0.2e1;
t293 = t113 / 0.2e1;
t291 = -t134 / 0.2e1;
t290 = t134 / 0.2e1;
t289 = -t137 / 0.2e1;
t288 = t137 / 0.2e1;
t287 = -t139 / 0.2e1;
t286 = t139 / 0.2e1;
t285 = -t176 / 0.2e1;
t284 = t176 / 0.2e1;
t283 = t178 / 0.2e1;
t282 = pkin(2) * t175;
t280 = pkin(7) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t241);
t273 = -qJD(2) / 0.2e1;
t271 = -Ifges(4,4) - Ifges(5,6);
t44 = mrSges(7,1) * t130 - mrSges(7,3) * t80;
t45 = mrSges(6,1) * t130 - mrSges(6,3) * t80;
t269 = -t44 - t45;
t46 = -mrSges(7,2) * t130 + mrSges(7,3) * t81;
t47 = -mrSges(6,2) * t130 + mrSges(6,3) * t81;
t268 = t46 + t47;
t149 = t175 * t177 - t218;
t191 = -t150 * qJ(4) + t169;
t75 = t149 * t297 + t191;
t108 = -t250 * t158 - t159 * t175;
t93 = pkin(4) * t150 + t108;
t28 = t176 * t93 + t178 * t75;
t82 = -mrSges(7,2) * t134 + mrSges(7,3) * t112;
t83 = -mrSges(6,2) * t134 + mrSges(6,3) * t112;
t267 = t82 + t83;
t84 = mrSges(7,1) * t134 - mrSges(7,3) * t113;
t85 = mrSges(6,1) * t134 - mrSges(6,3) * t113;
t266 = -t84 - t85;
t259 = t108 * t76;
t256 = t130 * mrSges(5,1);
t252 = t139 * Ifges(4,4);
t251 = t139 * Ifges(5,6);
t249 = Ifges(3,5) * qJD(2);
t247 = qJ(6) * t149;
t77 = t250 * t123 + t175 * t181;
t239 = qJD(2) * t177;
t66 = -mrSges(7,1) * t112 + mrSges(7,2) * t113;
t234 = -t66 - t314;
t233 = -Ifges(4,4) / 0.2e1 - Ifges(5,6) / 0.2e1;
t232 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t231 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t230 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t65 = -qJD(2) * qJD(4) - t77;
t25 = -t81 * mrSges(7,1) + t80 * mrSges(7,2);
t224 = t179 * t235;
t166 = qJ(4) + t282;
t217 = -t75 - t247;
t91 = t135 * t175 - t250 * t136;
t214 = -t1 * t176 + t178 * t2;
t213 = -t176 * t6 + t178 * t5;
t212 = t11 * t178 - t176 * t9;
t211 = t11 * t176 + t178 * t9;
t209 = mrSges(6,1) * t176 + mrSges(6,2) * t178;
t207 = mrSges(7,1) * t176 + mrSges(7,2) * t178;
t194 = t15 * t178 + t16 * t176;
t193 = t15 * t176 - t16 * t178;
t35 = -pkin(4) * t129 - t65;
t140 = qJD(2) * t218 - t175 * t239;
t185 = pkin(2) * t239 - qJ(4) * t140 - qJD(4) * t150;
t34 = t138 * t297 + t185;
t59 = pkin(4) * t140 + t91;
t7 = t176 * t59 + t178 * t34 + t93 * t237 - t238 * t75;
t187 = -t138 * t178 + t149 * t238;
t92 = t135 * t250 + t175 * t136;
t109 = t175 * t158 - t159 * t250;
t184 = -t176 * t267 + t178 * t266;
t183 = t6 * mrSges(6,1) + t1 * mrSges(7,1) - t5 * mrSges(6,2) - t2 * mrSges(7,2);
t60 = -t138 * pkin(4) + t92;
t180 = -qJD(5) * t193 - t306;
t170 = Ifges(3,4) * t240;
t154 = pkin(5) * t176 + t166;
t147 = Ifges(3,1) * t241 + t170 + t249;
t144 = t244 * t176;
t133 = Ifges(4,4) * t137;
t132 = Ifges(5,6) * t137;
t128 = Ifges(6,3) * t130;
t127 = Ifges(7,3) * t130;
t126 = t130 * mrSges(5,3);
t125 = t130 * mrSges(4,2);
t103 = t149 * pkin(3) + t191;
t102 = -mrSges(5,2) * t137 - mrSges(5,3) * t139;
t101 = mrSges(4,1) * t137 + mrSges(4,2) * t139;
t99 = t139 * Ifges(4,1) + Ifges(4,5) * qJD(2) - t133;
t98 = -t137 * Ifges(4,2) + Ifges(4,6) * qJD(2) + t252;
t97 = Ifges(5,4) * qJD(2) - t139 * Ifges(5,2) + t132;
t96 = Ifges(5,5) * qJD(2) + t137 * Ifges(5,3) - t251;
t95 = -qJD(2) * pkin(3) + t190;
t94 = -t149 * pkin(4) + t109;
t90 = pkin(3) * t139 + t216;
t89 = t178 * t93;
t79 = t107 - t277;
t74 = Ifges(6,5) * t80;
t73 = Ifges(7,5) * t80;
t72 = Ifges(6,6) * t81;
t71 = Ifges(7,6) * t81;
t62 = pkin(3) * t138 + t185;
t61 = t149 * t229 + t109;
t56 = t178 * t59;
t50 = pkin(3) * t129 + t186;
t39 = t113 * Ifges(6,5) + t112 * Ifges(6,6) + t134 * Ifges(6,3);
t38 = t113 * Ifges(7,5) + t112 * Ifges(7,6) + t134 * Ifges(7,3);
t29 = pkin(5) * t187 + t60;
t27 = -t176 * t75 + t89;
t26 = -mrSges(6,1) * t81 + mrSges(6,2) * t80;
t23 = -t176 * t58 + t69;
t22 = t178 * t247 + t28;
t17 = -pkin(5) * t81 + t35;
t14 = pkin(5) * t150 + t176 * t217 + t89;
t8 = -qJD(5) * t28 - t176 * t34 + t56;
t4 = -qJ(6) * t187 + t149 * t236 + t7;
t3 = pkin(5) * t140 + t56 + t217 * t237 + (-qJ(6) * t138 - qJD(5) * t93 - qJD(6) * t149 - t34) * t176;
t12 = [(t230 * t134 + t231 * t112 + t232 * t113 + t233 * t137 + t99 / 0.2e1 - t97 / 0.2e1 + t302 + t39 / 0.2e1 + t38 / 0.2e1 + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t139 + (Ifges(4,5) / 0.2e1 - Ifges(5,4) / 0.2e1) * qJD(2) - t104 * mrSges(4,3) + t95 * mrSges(5,1)) * t140 - t312 * t91 + m(7) * (t1 * t14 + t11 * t4 + t17 * t61 + t2 * t22 + t29 * t32 + t3 * t9) + m(6) * (t15 * t8 + t16 * t7 + t27 * t6 + t28 * t5 + t35 * t94 + t60 * t63) + (-t98 / 0.2e1 + t96 / 0.2e1 - t105 * mrSges(4,3) + t100 * mrSges(5,1) + t233 * t139 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t137 + (-Ifges(4,6) / 0.2e1 + Ifges(5,5) / 0.2e1) * qJD(2) + t212 * mrSges(7,3) - t193 * mrSges(6,3) + t310 * t295 + t309 * t293 + t311 * t290 + t318 * t284 + t319 * t283 - t303) * t138 + (-t248 / 0.2e1 + (mrSges(3,1) * t301 - 0.3e1 / 0.2e1 * t264 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t179) * qJD(1) + (m(4) * (t155 + t242) + qJD(1) * (mrSges(4,1) * t149 + mrSges(4,2) * t150) + t101) * pkin(2) - t304) * t239 + t103 * (-t129 * mrSges(5,2) - t126) + (-t77 * mrSges(4,3) + t65 * mrSges(5,1) - t50 * mrSges(5,2) - t17 * t208 - t35 * t210 + (t176 * t232 + t178 * t231 + t271) * t130 + (Ifges(5,3) + Ifges(4,2)) * t129 + t214 * mrSges(7,3) + t213 * mrSges(6,3) + (-mrSges(6,3) * t194 - mrSges(7,3) * t211 + t207 * t32 + t209 * t63 + t285 * t319 + t290 * t328 + t293 * t326 + t295 * t327) * qJD(5) + t309 * t299 + t310 * t298 + t320 * t284 + (qJD(5) * t318 + t321) * t283) * t149 + t332 * (t108 * t130 - t109 * t129) + (-t50 * mrSges(5,3) + t127 / 0.2e1 + t128 / 0.2e1 + t74 / 0.2e1 + t73 / 0.2e1 + t72 / 0.2e1 + t71 / 0.2e1 + t231 * t81 + t232 * t80 + t332 * t76 + t271 * t129 + (Ifges(4,1) + Ifges(5,2) + t230) * t130 + t183) * t150 + t313 * t92 + t169 * (t129 * mrSges(4,1) + t125) + (t147 / 0.2e1 - t280 + t249 / 0.2e1 + (mrSges(3,2) * t301 + 0.3e1 / 0.2e1 * Ifges(3,4) * t179) * qJD(1)) * t179 * qJD(2) + m(5) * (-t100 * t92 + t103 * t50 - t109 * t65 + t62 * t87 + t91 * t95 + t259) + m(4) * (-t104 * t91 + t105 * t92 + t109 * t77 + t259) + t14 * t44 + t27 * t45 + t22 * t46 + t28 * t47 + t61 * t25 + t29 * t66 + t60 * t67 + t4 * t82 + t7 * t83 + t3 * t84 + t8 * t85 + t94 * t26 + t62 * t102; t327 * t298 + t329 * qJD(5) - t101 * t171 + t320 * t283 + t321 * t285 + t304 * t241 + ((t175 * t77 - t250 * t76) * pkin(2) + t104 * t106 - t105 * t107 - t155 * t171) * m(4) + (-t133 + t99 + t39 + t38) * t288 + (-mrSges(5,1) * t166 - mrSges(4,3) * t282 + t322) * t129 + (-Ifges(4,2) * t288 + Ifges(5,3) * t289 + t273 * t322 + t291 * t311 + t294 * t309 + t296 * t310 + t303) * t139 + t324 * t76 - (Ifges(3,5) * t179 - Ifges(3,6) * t177) * t235 / 0.2e1 + (t132 + t97) * t289 - Ifges(3,6) * t225 + (-Ifges(4,1) * t287 + Ifges(5,2) * t286 - t337 * t296 - t331 * t294 + (-Ifges(6,3) - Ifges(7,3)) * t291 + t323 * t273 + t302) * t137 + (-mrSges(3,1) * t224 + mrSges(3,2) * t225) * pkin(7) + (-t15 * t23 - t16 * t24 + t166 * t35 + (qJD(4) - t79) * t63) * m(6) - (-Ifges(3,2) * t241 + t147 + t170) * t240 / 0.2e1 + (-t177 * (Ifges(3,1) * t179 - t264) / 0.2e1 + pkin(1) * (mrSges(3,1) * t177 + mrSges(3,2) * t179)) * qJD(1) ^ 2 + t314 * qJD(4) + t315 * t66 + t316 * t82 + t317 * t84 + (-t1 * t145 + t11 * t316 - t144 * t2 + t154 * t17 + t315 * t32 + t317 * t9) * m(7) + (t251 + t98) * t286 + (t11 * t308 + t307 * t9 + t305) * mrSges(7,3) + (t15 * t307 + t16 * t308 + t306) * mrSges(6,3) - (t112 * t310 + t113 * t309 + t134 * t311) * qJD(5) / 0.2e1 + t312 * t106 - t313 * t107 + (-t252 + t96) * t287 + (t100 * t325 - t106 * t95 - t166 * t65 + t168 * t76 - t87 * t90) * m(5) + t326 * t299 + t17 * t207 + t35 * t209 + (-t238 / 0.2e1 - t246 / 0.2e1) * t318 + (t328 / 0.2e1 - mrSges(4,3) * t228 - t323) * t130 + t166 * t26 + t154 * t25 - t144 * t46 - t145 * t44 + t240 * t280 + t105 * t253 + t168 * t256 + t95 * t265 + (-t237 / 0.2e1 - t245 / 0.2e1) * t319 - t100 * t254 - t104 * t255 + Ifges(3,5) * t224 + (m(6) * t180 + t176 * t47 + t178 * t45 + t237 * t83 - t238 * t85) * t164 - t65 * mrSges(5,3) - t77 * mrSges(4,2) - t79 * t67 - t24 * t83 - t23 * t85 - t90 * t102; t125 - t126 + t268 * t178 + t269 * t176 - t324 * t129 + (t119 - t234) * t137 + t184 * qJD(5) + (t184 + t312) * t139 + (-t134 * t211 + t137 * t32 + t214) * m(7) + (-t134 * t194 + t137 * t63 + t213) * m(6) + (-t100 * t137 - t139 * t95 + t50) * m(5) + (t104 * t139 + t105 * t137 + t165) * m(4); t256 + t139 * t102 + t234 * qJD(2) + (t134 * t267 - t269) * t178 + (t134 * t266 + t268) * t176 + (-qJD(2) * t32 + t134 * t212 - t305) * m(7) + (-qJD(2) * t63 - t139 * t193 + t180) * m(6) + (qJD(2) * t100 + t139 * t87 + t76) * m(5); (t11 * t113 + t112 * t9) * mrSges(7,3) + (t112 * t15 + t113 * t16) * mrSges(6,3) + (-t113 * t66 + t44) * pkin(5) + t127 + t128 + t183 + t74 + t73 + t72 + t71 + (-(t10 - t9) * t11 + (-t113 * t32 + t1) * pkin(5)) * m(7) - t10 * t82 - t15 * t83 + t11 * t84 + t16 * t85 - t32 * (mrSges(7,1) * t113 + mrSges(7,2) * t112) - t63 * (mrSges(6,1) * t113 + mrSges(6,2) * t112) + (t112 * t340 - t335) * t294 + t319 * t293 + (t112 * t331 - t113 * t337) * t291 + (-t113 * t338 + t318 + t336) * t296; -t112 * t82 + t113 * t84 + 0.2e1 * (t17 / 0.2e1 + t11 * t296 + t9 * t293) * m(7) + t25;];
tauc  = t12(:);
