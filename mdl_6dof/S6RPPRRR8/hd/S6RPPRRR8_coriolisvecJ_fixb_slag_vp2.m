% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:35:25
% EndTime: 2019-03-09 02:35:35
% DurationCPUTime: 5.84s
% Computational Cost: add. (8481->478), mult. (19189->662), div. (0->0), fcn. (13649->8), ass. (0->218)
t179 = sin(pkin(10));
t180 = cos(pkin(10));
t267 = sin(qJ(4));
t268 = cos(qJ(4));
t148 = t267 * t179 - t268 * t180;
t183 = sin(qJ(5));
t280 = -pkin(9) - pkin(8);
t229 = qJD(5) * t280;
t223 = qJD(1) * t267;
t224 = qJD(1) * t268;
t141 = -t179 * t224 - t180 * t223;
t236 = t183 * t141;
t142 = -t179 * t223 + t180 * t224;
t116 = pkin(4) * t142 - pkin(8) * t141;
t185 = cos(qJ(5));
t181 = -pkin(1) - qJ(3);
t162 = qJD(1) * t181 + qJD(2);
t219 = -pkin(7) * qJD(1) + t162;
t137 = t219 * t179;
t138 = t219 * t180;
t98 = -t267 * t137 + t268 * t138;
t54 = t183 * t116 + t185 * t98;
t312 = pkin(9) * t236 + t183 * t229 - t54;
t264 = pkin(9) * t185;
t53 = t185 * t116 - t183 * t98;
t311 = -pkin(5) * t142 + t141 * t264 + t185 * t229 - t53;
t99 = t137 * t268 + t138 * t267;
t91 = qJD(4) * pkin(8) + t99;
t178 = qJD(1) * qJ(2);
t174 = qJD(3) + t178;
t175 = t179 * pkin(3);
t157 = qJD(1) * t175 + t174;
t92 = -pkin(4) * t141 - pkin(8) * t142 + t157;
t49 = -t183 * t91 + t185 * t92;
t50 = t183 * t92 + t185 * t91;
t200 = t50 * t183 + t49 * t185;
t205 = Ifges(6,5) * t185 - Ifges(6,6) * t183;
t257 = Ifges(6,4) * t185;
t207 = -Ifges(6,2) * t183 + t257;
t258 = Ifges(6,4) * t183;
t209 = Ifges(6,1) * t185 - t258;
t210 = mrSges(6,1) * t183 + mrSges(6,2) * t185;
t269 = t185 / 0.2e1;
t140 = qJD(5) - t141;
t271 = t140 / 0.2e1;
t127 = qJD(4) * t183 + t142 * t185;
t275 = t127 / 0.2e1;
t126 = qJD(4) * t185 - t142 * t183;
t276 = t126 / 0.2e1;
t259 = Ifges(6,4) * t127;
t61 = Ifges(6,2) * t126 + Ifges(6,6) * t140 + t259;
t123 = Ifges(6,4) * t126;
t62 = Ifges(6,1) * t127 + Ifges(6,5) * t140 + t123;
t90 = -qJD(4) * pkin(4) - t98;
t310 = t62 * t269 - t183 * t61 / 0.2e1 + t90 * t210 + t207 * t276 + t209 * t275 + t205 * t271 - t200 * mrSges(6,3);
t184 = cos(qJ(6));
t182 = sin(qJ(6));
t36 = pkin(9) * t126 + t50;
t246 = t182 * t36;
t35 = -pkin(9) * t127 + t49;
t31 = pkin(5) * t140 + t35;
t10 = t184 * t31 - t246;
t244 = t184 * t36;
t11 = t182 * t31 + t244;
t136 = t142 * qJD(4);
t132 = Ifges(7,3) * t136;
t217 = t184 * t126 - t127 * t182;
t75 = t126 * t182 + t127 * t184;
t265 = Ifges(7,4) * t75;
t134 = qJD(6) + t140;
t274 = -t134 / 0.2e1;
t282 = -t75 / 0.2e1;
t284 = -t217 / 0.2e1;
t71 = Ifges(7,4) * t217;
t34 = Ifges(7,1) * t75 + Ifges(7,5) * t134 + t71;
t59 = -t126 * pkin(5) + t90;
t309 = t132 + (Ifges(7,5) * t217 - Ifges(7,6) * t75) * t274 + (t10 * t217 + t11 * t75) * mrSges(7,3) + (-Ifges(7,2) * t75 + t34 + t71) * t284 - t59 * (mrSges(7,1) * t75 + mrSges(7,2) * t217) + (Ifges(7,1) * t217 - t265) * t282;
t308 = t148 * qJD(3);
t139 = Ifges(5,4) * t141;
t291 = t139 / 0.2e1 + t142 * Ifges(5,1) / 0.2e1;
t307 = t157 * mrSges(5,2) + Ifges(5,5) * qJD(4) + t291 + t310;
t232 = qJD(5) * t185;
t233 = qJD(5) * t183;
t149 = t179 * t268 + t180 * t267;
t190 = t149 * qJD(3);
t66 = -qJD(1) * t190 + t98 * qJD(4);
t135 = t149 * qJD(4) * qJD(1);
t231 = qJD(1) * qJD(2);
t95 = pkin(4) * t136 + pkin(8) * t135 + t231;
t18 = t183 * t95 + t185 * t66 + t92 * t232 - t233 * t91;
t83 = -qJD(5) * t127 + t135 * t183;
t12 = pkin(9) * t83 + t18;
t19 = -qJD(5) * t50 - t183 * t66 + t185 * t95;
t82 = qJD(5) * t126 - t135 * t185;
t9 = pkin(5) * t136 - pkin(9) * t82 + t19;
t2 = qJD(6) * t10 + t12 * t184 + t182 * t9;
t29 = t217 * qJD(6) + t182 * t83 + t184 * t82;
t3 = -qJD(6) * t11 - t12 * t182 + t184 * t9;
t30 = -qJD(6) * t75 - t182 * t82 + t184 * t83;
t306 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t29 + Ifges(7,6) * t30;
t286 = t29 / 0.2e1;
t285 = t30 / 0.2e1;
t33 = Ifges(7,2) * t217 + Ifges(7,6) * t134 + t265;
t304 = t33 / 0.2e1;
t272 = t136 / 0.2e1;
t158 = t280 * t183;
t159 = t280 * t185;
t125 = t158 * t182 - t159 * t184;
t303 = -qJD(6) * t125 - t312 * t182 + t311 * t184;
t124 = t158 * t184 + t159 * t182;
t302 = qJD(6) * t124 + t311 * t182 + t312 * t184;
t234 = t179 ^ 2 + t180 ^ 2;
t300 = mrSges(4,3) * t234;
t152 = t182 * t185 + t183 * t184;
t293 = qJD(5) + qJD(6);
t121 = t293 * t152;
t221 = qJD(4) * t267;
t222 = qJD(4) * t268;
t144 = -t179 * t221 + t180 * t222;
t195 = t182 * t183 - t184 * t185;
t296 = -t152 * qJD(1) - t121 * t149 - t144 * t195;
t111 = t195 * t149;
t295 = t195 * qJD(1) + t293 * t111 - t152 * t144;
t44 = -mrSges(7,1) * t217 + mrSges(7,2) * t75;
t294 = m(7) * t59 + t44;
t96 = t152 * t141;
t242 = t96 - t121;
t120 = t293 * t195;
t97 = t195 * t141;
t241 = t97 - t120;
t110 = t152 * t148;
t170 = qJ(2) + t175;
t115 = pkin(4) * t149 + pkin(8) * t148 + t170;
t261 = -pkin(7) + t181;
t155 = t261 * t179;
t156 = t261 * t180;
t119 = t155 * t268 + t156 * t267;
t117 = t185 * t119;
t58 = t183 * t115 + t117;
t118 = t267 * t155 - t268 * t156;
t202 = t18 * t185 - t183 * t19;
t93 = -mrSges(6,2) * t140 + mrSges(6,3) * t126;
t94 = mrSges(6,1) * t140 - mrSges(6,3) * t127;
t196 = -t183 * t93 - t185 * t94;
t290 = -m(6) * t200 + t196;
t289 = t19 * mrSges(6,1) - t18 * mrSges(6,2) + Ifges(6,5) * t82 + Ifges(6,6) * t83 + t306;
t288 = Ifges(7,4) * t286 + Ifges(7,2) * t285 + Ifges(7,6) * t272;
t287 = Ifges(7,1) * t286 + Ifges(7,4) * t285 + Ifges(7,5) * t272;
t283 = t217 / 0.2e1;
t281 = t75 / 0.2e1;
t273 = t134 / 0.2e1;
t270 = t183 / 0.2e1;
t266 = m(3) * qJ(2);
t260 = m(4) * qJD(3);
t67 = -t308 * qJD(1) + qJD(4) * t99;
t256 = t118 * t67;
t251 = t141 * Ifges(5,2);
t248 = t148 * t67;
t243 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t126 + mrSges(6,2) * t127 + t142 * mrSges(5,3);
t237 = t148 * t183;
t212 = mrSges(4,1) * t179 + mrSges(4,2) * t180;
t235 = -mrSges(5,1) * t141 + mrSges(5,2) * t142 + qJD(1) * t212;
t230 = -t44 - t243;
t220 = t136 * mrSges(5,1) - t135 * mrSges(5,2);
t143 = -t179 * t222 - t180 * t221;
t113 = pkin(4) * t144 - pkin(8) * t143 + qJD(2);
t84 = -t118 * qJD(4) - t190;
t218 = t185 * t113 - t183 * t84;
t57 = t185 * t115 - t119 * t183;
t216 = qJD(1) * t234;
t211 = -mrSges(6,1) * t185 + mrSges(6,2) * t183;
t208 = Ifges(6,1) * t183 + t257;
t206 = Ifges(6,2) * t185 + t258;
t204 = Ifges(6,5) * t183 + Ifges(6,6) * t185;
t203 = -t143 * t98 - t144 * t99;
t201 = t18 * t183 + t185 * t19;
t47 = pkin(5) * t149 + t148 * t264 + t57;
t51 = pkin(9) * t237 + t58;
t21 = -t182 * t51 + t184 * t47;
t22 = t182 * t47 + t184 * t51;
t199 = -t183 * t49 + t185 * t50;
t63 = mrSges(6,1) * t136 - mrSges(6,3) * t82;
t64 = -mrSges(6,2) * t136 + mrSges(6,3) * t83;
t198 = -t183 * t63 + t185 * t64;
t197 = -t183 * t94 + t185 * t93;
t129 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t141;
t194 = -t129 - t197;
t192 = -t143 * t183 + t148 * t232;
t25 = t183 * t113 + t115 * t232 - t119 * t233 + t185 * t84;
t85 = qJD(4) * t119 - t308;
t188 = t11 * mrSges(7,2) + t50 * mrSges(6,2) + Ifges(5,6) * qJD(4) + t142 * Ifges(5,4) + t251 / 0.2e1 - t134 * Ifges(7,3) - t75 * Ifges(7,5) - t217 * Ifges(7,6) - t140 * Ifges(6,3) - t127 * Ifges(6,5) - t126 * Ifges(6,6) - t10 * mrSges(7,1) - t157 * mrSges(5,1) - t49 * mrSges(6,1);
t186 = qJD(1) ^ 2;
t173 = -pkin(5) * t185 - pkin(4);
t133 = Ifges(6,3) * t136;
t112 = t195 * t148;
t109 = t152 * t149;
t88 = -pkin(5) * t237 + t118;
t68 = pkin(5) * t236 + t99;
t56 = mrSges(7,1) * t134 - mrSges(7,3) * t75;
t55 = -mrSges(7,2) * t134 + mrSges(7,3) * t217;
t52 = -pkin(5) * t192 + t85;
t48 = -mrSges(6,1) * t83 + mrSges(6,2) * t82;
t45 = -t83 * pkin(5) + t67;
t43 = t82 * Ifges(6,1) + t83 * Ifges(6,4) + t136 * Ifges(6,5);
t42 = t82 * Ifges(6,4) + t83 * Ifges(6,2) + t136 * Ifges(6,6);
t40 = -t120 * t148 - t143 * t152;
t38 = t293 * t110 - t195 * t143;
t26 = -t58 * qJD(5) + t218;
t24 = -mrSges(7,2) * t136 + mrSges(7,3) * t30;
t23 = mrSges(7,1) * t136 - mrSges(7,3) * t29;
t20 = pkin(9) * t192 + t25;
t17 = -t143 * t264 + pkin(5) * t144 + (-t117 + (-pkin(9) * t148 - t115) * t183) * qJD(5) + t218;
t14 = t184 * t35 - t246;
t13 = -t182 * t35 - t244;
t8 = -mrSges(7,1) * t30 + mrSges(7,2) * t29;
t5 = -qJD(6) * t22 + t17 * t184 - t182 * t20;
t4 = qJD(6) * t21 + t17 * t182 + t184 * t20;
t1 = [t40 * t304 + (m(4) * (t174 + t178) + m(5) * (qJD(1) * t170 + t157) + t235 + ((2 * mrSges(3,3)) + t212 + 0.2e1 * t266) * qJD(1)) * qJD(2) + (-t162 * t234 - t181 * t216) * t260 + (-t118 * t135 - t119 * t136 + t203) * mrSges(5,3) + (-t185 * t43 / 0.2e1 + t42 * t270 - mrSges(5,2) * t231 + Ifges(5,1) * t135 - t82 * t209 / 0.2e1 - t83 * t207 / 0.2e1 + (-mrSges(5,3) - t210) * t67 + t201 * mrSges(6,3) + (mrSges(6,3) * t199 + t204 * t271 + t206 * t276 + t208 * t275 + t211 * t90 + t269 * t61 + t270 * t62) * qJD(5) + (Ifges(5,4) - t205 / 0.2e1) * t136) * t148 + (mrSges(5,1) * t231 - t66 * mrSges(5,3) + t132 / 0.2e1 + t133 / 0.2e1 + Ifges(5,4) * t135 + (Ifges(7,3) / 0.2e1 + Ifges(5,2) + Ifges(6,3) / 0.2e1) * t136 + t289) * t149 + (-t10 * t38 + t11 * t40 + t110 * t2 - t112 * t3) * mrSges(7,3) + 0.2e1 * qJD(3) * t300 * qJD(1) + (Ifges(7,4) * t112 + Ifges(7,2) * t110) * t285 + (Ifges(7,1) * t112 + Ifges(7,4) * t110) * t286 + t112 * t287 + t110 * t288 + (Ifges(7,1) * t38 + Ifges(7,4) * t40) * t281 + (Ifges(7,4) * t38 + Ifges(7,2) * t40) * t283 + (t291 + t307) * t143 + t21 * t23 + t22 * t24 + t38 * t34 / 0.2e1 + t170 * t220 + t52 * t44 + t4 * t55 + t5 * t56 + t59 * (-mrSges(7,1) * t40 + mrSges(7,2) * t38) + t57 * t63 + t58 * t64 + (Ifges(7,5) * t112 + Ifges(7,6) * t110) * t272 + (Ifges(7,5) * t38 + Ifges(7,6) * t40) * t273 + t88 * t8 + t25 * t93 + t26 * t94 + t243 * t85 + t45 * (-mrSges(7,1) * t110 + mrSges(7,2) * t112) + t118 * t48 + (-t251 / 0.2e1 - t188) * t144 + m(5) * (t119 * t66 + t84 * t99 - t85 * t98 + t256) + m(6) * (t18 * t58 + t19 * t57 + t25 * t50 + t26 * t49 + t85 * t90 + t256) + t84 * t129 + m(7) * (t10 * t5 + t11 * t4 + t2 * t22 + t21 * t3 + t45 * t88 + t52 * t59); -t109 * t23 - t111 * t24 + t295 * t56 + t296 * t55 + (-mrSges(3,3) - t266) * t186 + (-mrSges(5,3) * t135 + t48 + t8) * t148 - t194 * t144 + t230 * t143 + (-t136 * mrSges(5,3) + t196 * qJD(5) + t198) * t149 + m(6) * (-t143 * t90 + t248 + t199 * t144 + (-t200 * qJD(5) + t202) * t149) + m(5) * (t149 * t66 - t203 + t248) + (-m(4) * t174 - m(5) * t157 - t234 * t260 - t235 + t290) * qJD(1) + (t295 * t10 - t109 * t3 + t296 * t11 - t111 * t2 - t143 * t59 + t148 * t45) * m(7); -t195 * t23 + t152 * t24 + t183 * t64 + t185 * t63 + t242 * t56 + t241 * t55 + t197 * qJD(5) + (m(4) + m(5)) * t231 - t186 * t300 + t230 * t142 + t194 * t141 - m(5) * (t141 * t99 - t142 * t98) + m(4) * t162 * t216 + t220 + (t242 * t10 + t241 * t11 - t142 * t59 + t152 * t2 - t195 * t3) * m(7) + (t140 * t199 - t142 * t90 + t201) * m(6); (t294 * t183 * pkin(5) + t310) * qJD(5) + (-pkin(4) * t67 - t49 * t53 - t50 * t54 - t90 * t99) * m(6) + (Ifges(7,5) * t152 - Ifges(7,6) * t195 + t204) * t272 + (Ifges(7,1) * t152 - Ifges(7,4) * t195) * t286 + (Ifges(7,4) * t152 - Ifges(7,2) * t195) * t285 + (-t10 * t241 + t11 * t242 - t152 * t3 - t195 * t2) * mrSges(7,3) + t45 * (mrSges(7,1) * t195 + mrSges(7,2) * t152) - t195 * t288 + (t97 / 0.2e1 - t120 / 0.2e1) * t34 + (t96 / 0.2e1 - t121 / 0.2e1) * t33 + (-Ifges(7,1) * t120 - Ifges(7,4) * t121) * t281 + (-Ifges(7,4) * t120 - Ifges(7,2) * t121) * t283 + (-Ifges(7,5) * t120 - Ifges(7,6) * t121) * t273 + (-Ifges(7,5) * t97 - Ifges(7,6) * t96) * t274 + (-Ifges(7,1) * t97 - Ifges(7,4) * t96) * t282 + (-Ifges(7,4) * t97 - Ifges(7,2) * t96) * t284 + t302 * t55 + (t303 * t10 + t302 * t11 + t124 * t3 + t125 * t2 + t173 * t45 - t59 * t68) * m(7) + t303 * t56 + (m(6) * t202 + t290 * qJD(5) + t198) * pkin(8) + (t99 * mrSges(5,3) + t188) * t142 + t152 * t287 + (t98 * mrSges(5,3) - t139 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t142 - t307) * t141 + t202 * mrSges(6,3) + t83 * t206 / 0.2e1 + t82 * t208 / 0.2e1 + (-mrSges(5,1) + t211) * t67 - pkin(4) * t48 - t66 * mrSges(5,2) + t42 * t269 + t43 * t270 - t68 * t44 - t54 * t93 - t53 * t94 + (-mrSges(7,1) * t242 + mrSges(7,2) * t241) * t59 - t243 * t99 + t124 * t23 + t125 * t24 - t98 * t129 - Ifges(5,5) * t135 - Ifges(5,6) * t136 + t173 * t8; -(-Ifges(6,2) * t127 + t123 + t62) * t126 / 0.2e1 + t75 * t304 + (t182 * t24 + t184 * t23 + m(7) * (t182 * t2 + t184 * t3) - t294 * t127 + (-t182 * t56 + t184 * t55 + m(7) * (-t10 * t182 + t11 * t184)) * qJD(6)) * pkin(5) + t289 + (t126 * t49 + t127 * t50) * mrSges(6,3) + t133 - m(7) * (t10 * t13 + t11 * t14) + t61 * t275 - t14 * t55 - t13 * t56 - t49 * t93 + t50 * t94 - t90 * (mrSges(6,1) * t127 + mrSges(6,2) * t126) - t127 * (Ifges(6,1) * t126 - t259) / 0.2e1 - t140 * (Ifges(6,5) * t126 - Ifges(6,6) * t127) / 0.2e1 + t309; -t10 * t55 + t11 * t56 + t33 * t281 + t306 + t309;];
tauc  = t1(:);
