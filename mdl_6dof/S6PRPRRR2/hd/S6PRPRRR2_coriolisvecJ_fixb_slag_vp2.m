% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:46
% EndTime: 2019-03-08 20:27:02
% DurationCPUTime: 6.88s
% Computational Cost: add. (6200->498), mult. (15472->705), div. (0->0), fcn. (11365->12), ass. (0->245)
t177 = sin(qJ(4));
t181 = cos(qJ(4));
t192 = -pkin(4) * t181 - pkin(9) * t177 - pkin(3);
t173 = cos(pkin(12));
t263 = pkin(2) * t173;
t142 = t192 - t263;
t206 = pkin(4) * t177 - pkin(9) * t181;
t154 = t206 * qJD(4);
t171 = sin(pkin(12));
t167 = pkin(2) * t171 + pkin(8);
t176 = sin(qJ(5));
t180 = cos(qJ(5));
t219 = qJD(5) * t180;
t220 = qJD(5) * t176;
t222 = qJD(4) * t180;
t54 = t142 * t219 + t176 * t154 + (-t177 * t222 - t181 * t220) * t167;
t172 = sin(pkin(6));
t178 = sin(qJ(2));
t182 = cos(qJ(2));
t126 = (t171 * t182 + t173 * t178) * t172;
t119 = qJD(1) * t126;
t226 = qJD(1) * t172;
t212 = t178 * t226;
t156 = t171 * t212;
t211 = t182 * t226;
t122 = t173 * t211 - t156;
t230 = t180 * t181;
t68 = t119 * t176 + t122 * t230;
t320 = t54 - t68;
t150 = t167 * t230;
t191 = pkin(5) * t177 - pkin(10) * t230;
t223 = qJD(4) * t176;
t228 = t177 * t167 * t223 + t180 * t154;
t261 = pkin(10) * t177;
t231 = t176 * t181;
t67 = t119 * t180 - t122 * t231;
t319 = t191 * qJD(4) + (-t150 + (-t142 + t261) * t176) * qJD(5) + t228 - t67;
t221 = qJD(4) * t181;
t187 = t176 * t221 + t177 * t219;
t318 = -pkin(10) * t187 + t320;
t281 = -pkin(10) - pkin(9);
t213 = qJD(5) * t281;
t224 = qJD(2) * t181;
t151 = t206 * qJD(2);
t155 = qJD(2) * pkin(2) + t211;
t108 = t171 * t155 + t173 * t212;
t105 = qJD(2) * pkin(8) + t108;
t174 = cos(pkin(6));
t162 = qJD(1) * t174 + qJD(3);
t77 = -t177 * t105 + t162 * t181;
t53 = t176 * t151 + t180 * t77;
t317 = -t53 + (pkin(10) * t224 + t213) * t176;
t52 = t180 * t151 - t176 * t77;
t316 = -qJD(2) * t191 + t180 * t213 - t52;
t315 = -Ifges(5,1) / 0.2e1;
t170 = Ifges(5,4) * t224;
t314 = -t170 / 0.2e1;
t313 = qJD(4) / 0.2e1;
t225 = qJD(2) * t177;
t145 = -t176 * t225 + t222;
t146 = t180 * t225 + t223;
t215 = mrSges(5,3) * t225;
t227 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t145 - mrSges(6,2) * t146 - t215;
t175 = sin(qJ(6));
t179 = cos(qJ(6));
t207 = t179 * t145 - t146 * t175;
t95 = t145 * t175 + t146 * t179;
t51 = -mrSges(7,1) * t207 + mrSges(7,2) * t95;
t216 = -t51 + t227;
t72 = -qJD(4) * pkin(4) - t77;
t56 = -pkin(5) * t145 + t72;
t278 = m(7) * t56;
t294 = -m(5) * t77 + m(6) * t72;
t312 = -t294 + t216 - t278;
t233 = t162 * t177;
t78 = t105 * t181 + t233;
t73 = qJD(4) * pkin(9) + t78;
t107 = t155 * t173 - t156;
t83 = qJD(2) * t192 - t107;
t34 = t176 * t83 + t180 * t73;
t29 = pkin(10) * t145 + t34;
t241 = t175 * t29;
t165 = qJD(5) - t224;
t33 = -t176 * t73 + t180 * t83;
t28 = -pkin(10) * t146 + t33;
t25 = pkin(5) * t165 + t28;
t10 = t179 * t25 - t241;
t240 = t179 * t29;
t11 = t175 * t25 + t240;
t218 = qJD(2) * qJD(4);
t208 = t177 * t218;
t163 = Ifges(7,3) * t208;
t264 = Ifges(7,4) * t95;
t161 = qJD(6) + t165;
t269 = -t161 / 0.2e1;
t283 = -t95 / 0.2e1;
t285 = -t207 / 0.2e1;
t89 = Ifges(7,4) * t207;
t44 = Ifges(7,1) * t95 + Ifges(7,5) * t161 + t89;
t311 = t163 + (Ifges(7,5) * t207 - Ifges(7,6) * t95) * t269 + (t10 * t207 + t11 * t95) * mrSges(7,3) + (-Ifges(7,2) * t95 + t44 + t89) * t285 - t56 * (mrSges(7,1) * t95 + mrSges(7,2) * t207) + (Ifges(7,1) * t207 - t264) * t283;
t210 = Ifges(5,5) * t313;
t131 = t180 * t142;
t82 = -t180 * t261 + t131 + (-t167 * t176 - pkin(5)) * t181;
t98 = t176 * t142 + t150;
t88 = -t176 * t261 + t98;
t39 = -t175 * t88 + t179 * t82;
t310 = qJD(6) * t39 + t319 * t175 + t318 * t179;
t40 = t175 * t82 + t179 * t88;
t309 = -qJD(6) * t40 - t318 * t175 + t319 * t179;
t236 = qJD(4) * t78;
t194 = t171 * t178 - t173 * t182;
t293 = qJD(2) * t172;
t121 = t194 * t293;
t114 = qJD(1) * t121;
t49 = -t114 * t177 + t236;
t308 = m(5) * (-t49 + t236);
t159 = t281 * t176;
t160 = t281 * t180;
t112 = t159 * t175 - t160 * t179;
t307 = -qJD(6) * t112 - t317 * t175 + t316 * t179;
t111 = t159 * t179 + t160 * t175;
t306 = qJD(6) * t111 + t316 * t175 + t317 * t179;
t104 = -qJD(2) * pkin(3) - t107;
t197 = t176 * t34 + t180 * t33;
t252 = Ifges(6,4) * t180;
t201 = -Ifges(6,2) * t176 + t252;
t253 = Ifges(6,4) * t176;
t203 = Ifges(6,1) * t180 - t253;
t204 = mrSges(6,1) * t176 + mrSges(6,2) * t180;
t250 = Ifges(6,6) * t176;
t251 = Ifges(6,5) * t180;
t265 = t180 / 0.2e1;
t266 = -t176 / 0.2e1;
t270 = t146 / 0.2e1;
t254 = Ifges(6,4) * t146;
t86 = t145 * Ifges(6,2) + t165 * Ifges(6,6) + t254;
t143 = Ifges(6,4) * t145;
t87 = t146 * Ifges(6,1) + t165 * Ifges(6,5) + t143;
t184 = -t197 * mrSges(6,3) + t145 * t201 / 0.2e1 + t203 * t270 + t72 * t204 + t165 * (-t250 + t251) / 0.2e1 + t86 * t266 + t87 * t265;
t305 = -t104 * mrSges(5,2) + t77 * mrSges(5,3) + t225 * t315 - t184 - t210 + t314;
t217 = qJD(4) * qJD(5);
t115 = t180 * t217 + (-t177 * t220 + t180 * t221) * qJD(2);
t48 = qJD(4) * t77 - t114 * t181;
t84 = (t154 + t119) * qJD(2);
t15 = -qJD(5) * t34 - t176 * t48 + t180 * t84;
t8 = pkin(5) * t208 - pkin(10) * t115 + t15;
t116 = -qJD(2) * t187 - t176 * t217;
t14 = t176 * t84 + t180 * t48 + t83 * t219 - t220 * t73;
t9 = pkin(10) * t116 + t14;
t2 = qJD(6) * t10 + t175 * t8 + t179 * t9;
t3 = -qJD(6) * t11 - t175 * t9 + t179 * t8;
t37 = qJD(6) * t207 + t115 * t179 + t116 * t175;
t38 = -qJD(6) * t95 - t115 * t175 + t116 * t179;
t304 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t37 + Ifges(7,6) * t38;
t43 = Ifges(7,2) * t207 + Ifges(7,6) * t161 + t264;
t302 = t43 / 0.2e1;
t209 = -Ifges(5,6) * qJD(4) / 0.2e1;
t296 = t51 + t278;
t193 = t175 * t176 - t179 * t180;
t129 = t193 * t177;
t198 = t14 * t180 - t15 * t176;
t292 = qJD(5) + qJD(6);
t289 = -t15 * mrSges(6,1) + t14 * mrSges(6,2) - Ifges(6,5) * t115 - Ifges(6,6) * t116 - t304;
t288 = 0.2e1 * m(5);
t287 = t37 / 0.2e1;
t286 = t38 / 0.2e1;
t284 = t207 / 0.2e1;
t282 = t95 / 0.2e1;
t276 = t115 / 0.2e1;
t275 = t116 / 0.2e1;
t148 = t175 * t180 + t176 * t179;
t128 = t148 * t177;
t274 = -t128 / 0.2e1;
t273 = -t129 / 0.2e1;
t272 = -t145 / 0.2e1;
t271 = -t146 / 0.2e1;
t268 = t161 / 0.2e1;
t267 = -t165 / 0.2e1;
t262 = pkin(5) * t176;
t18 = -mrSges(7,1) * t38 + mrSges(7,2) * t37;
t69 = -mrSges(6,1) * t116 + mrSges(6,2) * t115;
t256 = t18 + t69;
t255 = Ifges(5,4) * t177;
t102 = t126 * t177 - t174 * t181;
t249 = t102 * t49;
t189 = t193 * t181;
t124 = qJD(2) * t189;
t99 = t292 * t193;
t239 = -t124 + t99;
t120 = qJD(2) * t126;
t113 = qJD(1) * t120;
t125 = t194 * t172;
t235 = t113 * t125;
t100 = t292 * t148;
t190 = t148 * t181;
t123 = qJD(2) * t190;
t229 = -t100 + t123;
t214 = mrSges(5,3) * t224;
t205 = mrSges(6,1) * t180 - mrSges(6,2) * t176;
t202 = Ifges(6,1) * t176 + t252;
t200 = Ifges(6,2) * t180 + t253;
t199 = Ifges(6,5) * t176 + Ifges(6,6) * t180;
t103 = t126 * t181 + t174 * t177;
t65 = -t103 * t176 + t125 * t180;
t66 = t103 * t180 + t125 * t176;
t23 = -t175 * t66 + t179 * t65;
t24 = t175 * t65 + t179 * t66;
t90 = mrSges(6,1) * t208 - mrSges(6,3) * t115;
t91 = -mrSges(6,2) * t208 + mrSges(6,3) * t116;
t196 = -t176 * t90 + t180 * t91;
t117 = -mrSges(6,2) * t165 + mrSges(6,3) * t145;
t118 = mrSges(6,1) * t165 - mrSges(6,3) * t146;
t195 = -t176 * t117 - t180 * t118;
t185 = t10 * mrSges(7,1) + t104 * mrSges(5,1) + t33 * mrSges(6,1) + t209 - (t181 * Ifges(5,2) + t255) * qJD(2) / 0.2e1 + t161 * Ifges(7,3) + t95 * Ifges(7,5) + t207 * Ifges(7,6) + t165 * Ifges(6,3) + t146 * Ifges(6,5) + t145 * Ifges(6,6) - t11 * mrSges(7,2) - t34 * mrSges(6,2) - t78 * mrSges(5,3);
t169 = -pkin(5) * t180 - pkin(4);
t168 = -pkin(3) - t263;
t164 = Ifges(6,3) * t208;
t158 = -qJD(4) * mrSges(5,2) + t214;
t149 = (-mrSges(5,1) * t181 + mrSges(5,2) * t177) * qJD(2);
t140 = (mrSges(5,1) * t177 + mrSges(5,2) * t181) * t218;
t132 = (t167 + t262) * t177;
t106 = pkin(5) * t187 + t167 * t221;
t97 = -t167 * t231 + t131;
t75 = mrSges(7,1) * t161 - mrSges(7,3) * t95;
t74 = -mrSges(7,2) * t161 + mrSges(7,3) * t207;
t70 = t233 + (qJD(2) * t262 + t105) * t181;
t64 = -qJD(4) * t102 - t121 * t181;
t63 = qJD(4) * t103 - t121 * t177;
t60 = t115 * Ifges(6,1) + t116 * Ifges(6,4) + Ifges(6,5) * t208;
t59 = t115 * Ifges(6,4) + t116 * Ifges(6,2) + Ifges(6,6) * t208;
t58 = -qJD(4) * t190 + t129 * t292;
t57 = -qJD(4) * t189 - t100 * t177;
t55 = -qJD(5) * t98 + t228;
t32 = -mrSges(7,2) * t208 + mrSges(7,3) * t38;
t31 = mrSges(7,1) * t208 - mrSges(7,3) * t37;
t30 = -pkin(5) * t116 + t49;
t22 = qJD(5) * t65 + t120 * t176 + t180 * t64;
t21 = -qJD(5) * t66 + t120 * t180 - t176 * t64;
t17 = t37 * Ifges(7,1) + t38 * Ifges(7,4) + Ifges(7,5) * t208;
t16 = t37 * Ifges(7,4) + t38 * Ifges(7,2) + Ifges(7,6) * t208;
t13 = t179 * t28 - t241;
t12 = -t175 * t28 - t240;
t5 = -qJD(6) * t24 - t175 * t22 + t179 * t21;
t4 = qJD(6) * t23 + t175 * t21 + t179 * t22;
t1 = [t22 * t117 + t21 * t118 + t120 * t149 + t125 * t140 + t64 * t158 + t23 * t31 + t24 * t32 + t4 * t74 + t5 * t75 + t65 * t90 + t66 * t91 + t256 * t102 - t216 * t63 + m(7) * (t10 * t5 + t102 * t30 + t11 * t4 + t2 * t24 + t23 * t3 + t56 * t63) + m(6) * (t14 * t66 + t15 * t65 + t21 * t33 + t22 * t34 + t63 * t72 + t249) + m(5) * (t103 * t48 + t104 * t120 - t63 * t77 + t64 * t78 + t235 + t249) + m(4) * (-t107 * t120 - t108 * t121 - t114 * t126 + t235) + (-t120 * mrSges(4,1) + t121 * mrSges(4,2) + (t102 * t181 - t103 * t177) * qJD(4) * mrSges(5,3) + (-mrSges(3,1) * t178 - mrSges(3,2) * t182) * t293) * qJD(2); t320 * t117 + (qJD(2) * t119 - t113) * mrSges(4,1) + t58 * t302 + (t201 * t275 + t113 * mrSges(5,2) + t203 * t276 + t59 * t266 + t60 * t265 + (mrSges(5,3) + t204) * t49 + (-t14 * t176 - t15 * t180) * mrSges(6,3) + (-t180 * t86 / 0.2e1 + t87 * t266 + t199 * t267 + t200 * t272 + t202 * t271 + t72 * t205 + (t33 * t176 - t34 * t180) * mrSges(6,3)) * qJD(5) + (t185 + t209 + (Ifges(7,5) * t273 + Ifges(7,6) * t274 + (t251 / 0.2e1 - t250 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,4)) * t177 + (-Ifges(7,3) / 0.2e1 - Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(5,1)) * t181) * qJD(2)) * qJD(4) + (m(6) * t49 - qJD(4) * t158 - t308 + t69) * t167 + t312 * t122) * t177 + m(6) * (t14 * t98 + t15 * t97 + t33 * t55 + t34 * t54) + (-t163 / 0.2e1 - t164 / 0.2e1 - t113 * mrSges(5,1) + t48 * mrSges(5,3) - t122 * t158 + (-t122 * t78 / 0.2e1 + t167 * t48 / 0.2e1) * t288 + (0.3e1 / 0.2e1 * t170 + t210 + (-t227 + t294) * t167 - t305) * qJD(4) + t289) * t181 + (-t10 * t57 + t11 * t58 - t128 * t2 + t129 * t3) * mrSges(7,3) + t30 * (mrSges(7,1) * t128 - mrSges(7,2) * t129) + (-Ifges(7,4) * t129 - Ifges(7,2) * t128) * t286 + (-Ifges(7,1) * t129 - Ifges(7,4) * t128) * t287 + (t55 - t67) * t118 + t168 * t140 - t119 * t149 - m(6) * (t33 * t67 + t34 * t68) + t132 * t18 + t309 * t75 + t310 * t74 + (t10 * t309 + t106 * t56 + t11 * t310 + t132 * t30 + t2 * t40 + t3 * t39) * m(7) + (Ifges(7,5) * t57 + Ifges(7,6) * t58) * t268 + t17 * t273 + t16 * t274 + (Ifges(7,1) * t57 + Ifges(7,4) * t58) * t282 + (Ifges(7,4) * t57 + Ifges(7,2) * t58) * t284 + (t113 * t168 / 0.2e1 - t104 * t119 / 0.2e1) * t288 + (qJD(2) * t122 + t114) * mrSges(4,2) + (t107 * t119 - t108 * t122 + (-t113 * t173 - t114 * t171) * pkin(2)) * m(4) + t39 * t31 + t40 * t32 + t57 * t44 / 0.2e1 + t56 * (-mrSges(7,1) * t58 + mrSges(7,2) * t57) + t97 * t90 + t98 * t91 + t106 * t51; t57 * t74 + t58 * t75 - t129 * t32 - t128 * t31 + m(7) * (t10 * t58 + t11 * t57 - t128 * t3 - t129 * t2) + ((t117 * t180 - t118 * t176 + t158 - t214) * qJD(4) + m(6) * (t222 * t34 - t223 * t33 - t49) - m(7) * t30 + t308 - t256) * t181 + (t195 * qJD(5) + m(6) * (-t219 * t33 - t220 * t34 + t198) + m(5) * t48 + t196 + (-t215 - t312) * qJD(4)) * t177; t196 * pkin(9) + m(6) * (-pkin(4) * t49 + pkin(9) * t198) + (t123 / 0.2e1 - t100 / 0.2e1) * t43 + (-Ifges(7,4) * t124 - Ifges(7,2) * t123) * t285 + (t124 / 0.2e1 - t99 / 0.2e1) * t44 + (-Ifges(7,5) * t124 - Ifges(7,6) * t123) * t269 + (-Ifges(7,1) * t124 - Ifges(7,4) * t123) * t283 + ((t210 + t314 + t305) * t181 + (-t185 + (t255 / 0.2e1 + (t315 + Ifges(5,2) / 0.2e1) * t181) * qJD(2) + t209 + (Ifges(7,5) * t148 - Ifges(7,6) * t193 + t199) * t313) * t177) * qJD(2) + (-mrSges(5,1) - t205) * t49 + t306 * t74 + (t10 * t307 + t11 * t306 + t111 * t3 + t112 * t2 + t169 * t30 - t56 * t70) * m(7) + t307 * t75 + t198 * mrSges(6,3) + (t296 * t262 + (-m(6) * t197 + t195) * pkin(9) + t184) * qJD(5) - m(6) * (t33 * t52 + t34 * t53 + t72 * t78) + t176 * t60 / 0.2e1 + t169 * t18 - t77 * t158 + t148 * t17 / 0.2e1 + t59 * t265 + t200 * t275 + t202 * t276 + t227 * t78 + (-mrSges(7,1) * t229 - mrSges(7,2) * t239) * t56 + t30 * (mrSges(7,1) * t193 + mrSges(7,2) * t148) + (Ifges(7,4) * t148 - Ifges(7,2) * t193) * t286 + (Ifges(7,1) * t148 - Ifges(7,4) * t193) * t287 + (t10 * t239 + t11 * t229 - t148 * t3 - t193 * t2) * mrSges(7,3) - t193 * t16 / 0.2e1 + (-Ifges(7,5) * t99 - Ifges(7,6) * t100) * t268 + (-Ifges(7,1) * t99 - Ifges(7,4) * t100) * t282 + (-Ifges(7,4) * t99 - Ifges(7,2) * t100) * t284 - t48 * mrSges(5,2) - pkin(4) * t69 - t70 * t51 + t111 * t31 + t112 * t32 - t53 * t117 - t52 * t118; t164 - m(7) * (t10 * t12 + t11 * t13) - t289 + (t175 * t32 + t179 * t31 + m(7) * (t175 * t2 + t179 * t3) - t296 * t146 + (-t175 * t75 + t179 * t74 + m(7) * (-t10 * t175 + t11 * t179)) * qJD(6)) * pkin(5) + (-Ifges(6,2) * t146 + t143 + t87) * t272 - t72 * (mrSges(6,1) * t146 + mrSges(6,2) * t145) + (t145 * t33 + t146 * t34) * mrSges(6,3) + (Ifges(6,5) * t145 - Ifges(6,6) * t146) * t267 + t86 * t270 + (Ifges(6,1) * t145 - t254) * t271 + t95 * t302 - t13 * t74 - t12 * t75 - t33 * t117 + t34 * t118 + t311; -t10 * t74 + t11 * t75 + t43 * t282 + t304 + t311;];
tauc  = t1(:);
