% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2018-11-23 14:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:55:52
% EndTime: 2018-11-23 14:55:57
% DurationCPUTime: 5.03s
% Computational Cost: add. (4221->453), mult. (11003->659), div. (0->0), fcn. (8166->12), ass. (0->217)
t165 = sin(pkin(11));
t166 = sin(pkin(6));
t168 = cos(pkin(11));
t172 = sin(qJ(2));
t175 = cos(qJ(2));
t120 = (t165 * t175 + t168 * t172) * t166;
t107 = qJD(1) * t120;
t216 = qJD(1) * t166;
t205 = t172 * t216;
t149 = t165 * t205;
t204 = t175 * t216;
t109 = t168 * t204 - t149;
t167 = cos(pkin(12));
t164 = sin(pkin(12));
t174 = cos(qJ(4));
t222 = t164 * t174;
t55 = t107 * t167 - t109 * t222;
t171 = sin(qJ(4));
t192 = pkin(4) * t171 - qJ(5) * t174;
t128 = qJD(4) * t192 - qJD(5) * t171;
t160 = pkin(2) * t165 + pkin(8);
t213 = qJD(4) * t171;
t203 = t160 * t213;
t84 = t167 * t128 + t164 * t203;
t274 = t84 - t55;
t220 = t167 * t174;
t185 = pkin(5) * t171 - pkin(9) * t220;
t181 = t185 * qJD(4);
t273 = t181 + t274;
t116 = t164 * t128;
t210 = pkin(9) * t222;
t221 = t167 * t171;
t56 = t107 * t164 + t109 * t220;
t272 = t56 - t116 - (-t160 * t221 - t210) * qJD(4);
t242 = pkin(9) + qJ(5);
t150 = t242 * t164;
t151 = t242 * t167;
t170 = sin(qJ(6));
t173 = cos(qJ(6));
t102 = -t150 * t170 + t151 * t173;
t143 = t164 * t173 + t167 * t170;
t146 = t192 * qJD(2);
t169 = cos(pkin(6));
t156 = qJD(1) * t169 + qJD(3);
t148 = qJD(2) * pkin(2) + t204;
t104 = t165 * t148 + t168 * t205;
t97 = qJD(2) * pkin(8) + t104;
t90 = t171 * t97;
t72 = t156 * t174 - t90;
t43 = t167 * t146 - t164 * t72;
t36 = qJD(2) * t185 + t43;
t44 = t164 * t146 + t167 * t72;
t39 = -qJD(2) * t210 + t44;
t271 = -qJD(5) * t143 - qJD(6) * t102 + t170 * t39 - t173 * t36;
t270 = qJD(4) / 0.2e1;
t101 = -t150 * t173 - t151 * t170;
t187 = t164 * t170 - t167 * t173;
t269 = -qJD(5) * t187 + qJD(6) * t101 - t170 * t36 - t173 * t39;
t215 = qJD(2) * t171;
t268 = t215 / 0.2e1;
t202 = -Ifges(5,6) * qJD(4) / 0.2e1;
t267 = Ifges(5,5) * t270;
t184 = -pkin(4) * t174 - qJ(5) * t171 - pkin(3);
t248 = pkin(2) * t168;
t138 = t184 - t248;
t124 = t167 * t138;
t75 = -pkin(9) * t221 + t124 + (-t160 * t164 - pkin(5)) * t174;
t92 = t164 * t138 + t160 * t220;
t79 = -pkin(9) * t164 * t171 + t92;
t27 = t170 * t75 + t173 * t79;
t266 = -qJD(6) * t27 + t272 * t170 + t273 * t173;
t26 = -t170 * t79 + t173 * t75;
t265 = qJD(6) * t26 + t273 * t170 - t272 * t173;
t140 = qJD(4) * t167 - t164 * t215;
t141 = t164 * qJD(4) + t167 * t215;
t198 = t173 * t140 - t141 * t170;
t264 = qJD(2) * t166;
t130 = t187 * qJD(6);
t186 = t165 * t172 - t168 * t175;
t110 = t186 * t264;
t106 = qJD(1) * t110;
t212 = qJD(4) * t174;
t219 = -t174 * t106 + t156 * t212;
t35 = (qJD(5) - t90) * qJD(4) + t219;
t67 = (t107 + t128) * qJD(2);
t15 = -t164 * t35 + t167 * t67;
t11 = qJD(2) * t181 + t15;
t16 = t164 * t67 + t167 * t35;
t211 = qJD(2) * qJD(4);
t200 = t174 * t211;
t197 = t164 * t200;
t12 = -pkin(9) * t197 + t16;
t214 = qJD(2) * t174;
t223 = t156 * t171;
t73 = t174 * t97 + t223;
t64 = qJD(4) * qJ(5) + t73;
t103 = t148 * t168 - t149;
t78 = qJD(2) * t184 - t103;
t24 = -t164 * t64 + t167 * t78;
t20 = -pkin(5) * t214 - pkin(9) * t141 + t24;
t25 = t164 * t78 + t167 * t64;
t23 = pkin(9) * t140 + t25;
t5 = -t170 * t23 + t173 * t20;
t1 = qJD(6) * t5 + t11 * t170 + t12 * t173;
t6 = t170 * t20 + t173 * t23;
t2 = -qJD(6) * t6 + t11 * t173 - t12 * t170;
t182 = t187 * t174;
t179 = qJD(4) * t182;
t49 = -qJD(2) * t179 + qJD(6) * t198;
t183 = t143 * t174;
t180 = qJD(4) * t183;
t89 = t140 * t170 + t141 * t173;
t50 = -qJD(2) * t180 - qJD(6) * t89;
t263 = -t2 * mrSges(7,1) + t1 * mrSges(7,2) - Ifges(7,5) * t49 - Ifges(7,6) * t50;
t163 = Ifges(5,4) * t214;
t237 = Ifges(6,2) * t164;
t239 = Ifges(6,4) * t167;
t193 = -t237 + t239;
t240 = Ifges(6,4) * t164;
t194 = Ifges(6,1) * t167 - t240;
t241 = mrSges(6,2) * t167;
t195 = mrSges(6,1) * t164 + t241;
t250 = t167 / 0.2e1;
t251 = -t164 / 0.2e1;
t62 = -qJD(4) * pkin(4) + qJD(5) - t72;
t96 = -qJD(2) * pkin(3) - t103;
t262 = -(t164 * t25 + t167 * t24) * mrSges(6,3) + t62 * t195 + t96 * mrSges(5,2) + Ifges(5,1) * t268 + t163 / 0.2e1 + t267 + t140 * t193 / 0.2e1 + t141 * t194 / 0.2e1 + (t141 * Ifges(6,4) + t140 * Ifges(6,2) - Ifges(6,6) * t214) * t251 + (t141 * Ifges(6,1) + t140 * Ifges(6,4) - Ifges(6,5) * t214) * t250 - t72 * mrSges(5,3);
t261 = t49 / 0.2e1;
t260 = t50 / 0.2e1;
t259 = -t198 / 0.2e1;
t258 = t198 / 0.2e1;
t257 = -t89 / 0.2e1;
t256 = t89 / 0.2e1;
t121 = t143 * t171;
t255 = -t121 / 0.2e1;
t122 = t187 * t171;
t254 = -t122 / 0.2e1;
t158 = qJD(6) - t214;
t253 = -t158 / 0.2e1;
t252 = t158 / 0.2e1;
t249 = Ifges(7,4) * t89;
t247 = pkin(5) * t164;
t226 = t106 * t171;
t38 = qJD(4) * t73 - t226;
t94 = t120 * t171 - t169 * t174;
t246 = t38 * t94;
t238 = Ifges(6,5) * t167;
t236 = Ifges(6,6) * t164;
t232 = t174 * t38;
t115 = mrSges(6,1) * t197 + t200 * t241;
t17 = -t50 * mrSges(7,1) + t49 * mrSges(7,2);
t231 = t115 + t17;
t207 = mrSges(5,3) * t215;
t230 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t140 + mrSges(6,2) * t141 + t207;
t108 = qJD(2) * t120;
t105 = qJD(1) * t108;
t119 = t186 * t166;
t227 = t105 * t119;
t225 = t109 * t171;
t113 = qJD(2) * t183;
t131 = t143 * qJD(6);
t218 = -t113 + t131;
t114 = qJD(2) * t182;
t217 = -t114 + t130;
t40 = -mrSges(7,1) * t198 + mrSges(7,2) * t89;
t209 = -t40 - t230;
t208 = t160 * t171 * t38;
t206 = mrSges(5,3) * t214;
t201 = t171 * t211;
t199 = t160 + t247;
t191 = -t15 * t164 + t16 * t167;
t189 = -t164 * t24 + t167 * t25;
t95 = t120 * t174 + t169 * t171;
t53 = t119 * t167 - t164 * t95;
t54 = t119 * t164 + t167 * t95;
t18 = -t170 * t54 + t173 * t53;
t19 = t170 * t53 + t173 * t54;
t188 = -t171 * t72 + t174 * t73;
t57 = t223 + (qJD(2) * t247 + t97) * t174;
t178 = t24 * mrSges(6,1) + t5 * mrSges(7,1) + t96 * mrSges(5,1) + t202 - (Ifges(5,4) * t171 + t174 * Ifges(5,2)) * qJD(2) / 0.2e1 + t158 * Ifges(7,3) + t89 * Ifges(7,5) + t198 * Ifges(7,6) - Ifges(6,3) * t214 / 0.2e1 + t141 * Ifges(6,5) + t140 * Ifges(6,6) - t25 * mrSges(6,2) - t6 * mrSges(7,2) - t73 * mrSges(5,3);
t162 = -pkin(5) * t167 - pkin(4);
t161 = -pkin(3) - t248;
t157 = Ifges(7,3) * t201;
t155 = -qJD(4) * mrSges(5,2) + t206;
t145 = (-mrSges(5,1) * t174 + mrSges(5,2) * t171) * qJD(2);
t136 = (mrSges(5,1) * t171 + mrSges(5,2) * t174) * t211;
t127 = t199 * t171;
t126 = (mrSges(6,1) * t171 - mrSges(6,3) * t220) * t211;
t125 = (-mrSges(6,2) * t171 - mrSges(6,3) * t222) * t211;
t118 = t199 * t212;
t112 = -mrSges(6,1) * t214 - mrSges(6,3) * t141;
t111 = mrSges(6,2) * t214 + mrSges(6,3) * t140;
t99 = (Ifges(6,5) * t171 + t174 * t194) * t211;
t98 = (Ifges(6,6) * t171 + t174 * t193) * t211;
t91 = -t160 * t222 + t124;
t85 = -t167 * t203 + t116;
t83 = Ifges(7,4) * t198;
t77 = t130 * t171 - t180;
t76 = -t131 * t171 - t179;
t69 = mrSges(7,1) * t158 - mrSges(7,3) * t89;
t68 = -mrSges(7,2) * t158 + mrSges(7,3) * t198;
t52 = qJD(4) * t95 - t110 * t171;
t51 = -qJD(4) * t94 - t110 * t174;
t45 = -pkin(5) * t140 + t62;
t42 = -mrSges(7,2) * t201 + mrSges(7,3) * t50;
t41 = mrSges(7,1) * t201 - mrSges(7,3) * t49;
t37 = -t213 * t97 + t219;
t34 = Ifges(7,1) * t89 + Ifges(7,5) * t158 + t83;
t33 = Ifges(7,2) * t198 + Ifges(7,6) * t158 + t249;
t30 = qJD(4) * t57 - t226;
t29 = t108 * t164 + t167 * t51;
t28 = t108 * t167 - t164 * t51;
t14 = Ifges(7,1) * t49 + Ifges(7,4) * t50 + Ifges(7,5) * t201;
t13 = Ifges(7,4) * t49 + Ifges(7,2) * t50 + Ifges(7,6) * t201;
t4 = -qJD(6) * t19 - t170 * t29 + t173 * t28;
t3 = qJD(6) * t18 + t170 * t28 + t173 * t29;
t7 = [t108 * t145 + t29 * t111 + t28 * t112 + t119 * t136 + t54 * t125 + t53 * t126 + t51 * t155 + t18 * t41 + t19 * t42 + t3 * t68 + t4 * t69 + t231 * t94 - t209 * t52 + m(7) * (t1 * t19 + t18 * t2 + t3 * t6 + t30 * t94 + t4 * t5 + t45 * t52) + m(6) * (t15 * t53 + t16 * t54 + t24 * t28 + t25 * t29 + t52 * t62 + t246) + m(5) * (t108 * t96 + t37 * t95 + t51 * t73 - t52 * t72 + t227 + t246) + m(4) * (-t103 * t108 - t104 * t110 - t106 * t120 + t227) + (-t108 * mrSges(4,1) + t110 * mrSges(4,2) + (-t171 * t95 + t174 * t94) * qJD(4) * mrSges(5,3) + (-mrSges(3,1) * t172 - mrSges(3,2) * t175) * t264) * qJD(2); (-t1 * t121 + t122 * t2 - t5 * t76 + t6 * t77) * mrSges(7,3) + (-t157 / 0.2e1 - t105 * mrSges(5,1) + t16 * mrSges(6,2) - t15 * mrSges(6,1) + t37 * mrSges(5,3) - t109 * t155 + t263) * t174 - m(5) * (t107 * t96 + t109 * t188) + ((-t105 * t168 - t106 * t165) * pkin(2) + t103 * t107 - t104 * t109) * m(4) + m(6) * (t15 * t91 + t16 * t92 + t24 * t84 + t25 * t85 + t208) + m(5) * (t160 * t174 * t37 + t105 * t161 + t208) + (t105 * mrSges(5,2) + t160 * t115 + t98 * t251 + t99 * t250 + (mrSges(5,3) + t195) * t38 + (-t15 * t167 - t16 * t164) * mrSges(6,3) + t209 * t109) * t171 + (t107 * mrSges(4,1) + t109 * mrSges(4,2) + (Ifges(7,5) * t254 + Ifges(7,6) * t255 + (t238 / 0.2e1 - t236 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,4)) * t171) * t213 + ((0.3e1 / 0.2e1 * Ifges(5,4) - 0.3e1 / 0.2e1 * t238 + 0.3e1 / 0.2e1 * t236) * t174 + (-0.3e1 / 0.2e1 * Ifges(6,3) - Ifges(7,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2) + Ifges(6,1) * t167 ^ 2 / 0.2e1 + (-t239 + t237 / 0.2e1) * t164) * t171) * t212) * qJD(2) + t30 * (mrSges(7,1) * t121 - mrSges(7,2) * t122) + (-Ifges(7,4) * t122 - Ifges(7,2) * t121) * t260 + (-Ifges(7,1) * t122 - Ifges(7,4) * t121) * t261 + (((-m(5) * t73 - t155) * t160 + t178 + t202) * t171 + (t267 + (-m(5) * t72 + m(6) * t62 + t230) * t160 + t262) * t174) * qJD(4) + t265 * t68 + t266 * t69 + (t1 * t27 + t127 * t30 + t2 * t26 + t265 * t6 + t266 * t5 + (t118 - t225) * t45) * m(7) - m(6) * (t225 * t62 + t24 * t55 + t25 * t56) + (t85 - t56) * t111 - t105 * mrSges(4,1) + t106 * mrSges(4,2) + t118 * t40 + t92 * t125 + t91 * t126 + t127 * t17 + t76 * t34 / 0.2e1 + t45 * (-mrSges(7,1) * t77 + mrSges(7,2) * t76) + t77 * t33 / 0.2e1 + t26 * t41 + t27 * t42 - t107 * t145 + t161 * t136 + (Ifges(7,5) * t76 + Ifges(7,6) * t77) * t252 + t14 * t254 + t13 * t255 + (Ifges(7,1) * t76 + Ifges(7,4) * t77) * t256 + (Ifges(7,4) * t76 + Ifges(7,2) * t77) * t258 + t274 * t112; -t121 * t41 - t122 * t42 + t76 * t68 + t77 * t69 - t231 * t174 + (t125 * t167 - t126 * t164) * t171 + ((t111 * t167 - t112 * t164 + t155 - t206) * t174 + (-t207 - t209) * t171) * qJD(4) + m(6) * (-t232 + t191 * t171 + (t171 * t62 + t174 * t189) * qJD(4)) + m(7) * (-t1 * t122 - t121 * t2 - t174 * t30 + t213 * t45 + t5 * t77 + t6 * t76) + m(5) * (qJD(4) * t188 + t171 * t37 - t232); (mrSges(7,1) * t218 - mrSges(7,2) * t217) * t45 + (-pkin(4) * t38 + qJ(5) * t191 + qJD(5) * t189 - t24 * t43 - t25 * t44 - t62 * t73) * m(6) - t187 * t13 / 0.2e1 + (t113 / 0.2e1 - t131 / 0.2e1) * t33 + (-Ifges(7,5) * t114 - Ifges(7,6) * t113) * t253 + (-Ifges(7,1) * t114 - Ifges(7,4) * t113) * t257 + (-Ifges(7,4) * t114 - Ifges(7,2) * t113) * t259 + (t114 / 0.2e1 - t130 / 0.2e1) * t34 + (t1 * t102 + t101 * t2 + t162 * t30 + t269 * t6 + t271 * t5 - t45 * t57) * m(7) + t271 * t69 + t269 * t68 + ((-t163 / 0.2e1 + (-t236 + t238) * t214 / 0.2e1 + (Ifges(5,5) / 0.2e1 + (Ifges(6,1) * t164 + t239) * t250 + (Ifges(6,2) * t167 + t240) * t251) * qJD(4) - t262) * t174 + (Ifges(5,4) * t268 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t214 - t178 + t202 + (Ifges(6,5) * t164 + Ifges(7,5) * t143 + Ifges(6,6) * t167 - Ifges(7,6) * t187) * t270) * t171) * qJD(2) + (t16 * mrSges(6,3) + qJD(5) * t111 + qJ(5) * t125 + t98 / 0.2e1 - t38 * mrSges(6,1)) * t167 + (-t15 * mrSges(6,3) - qJD(5) * t112 - qJ(5) * t126 + t99 / 0.2e1 + t38 * mrSges(6,2)) * t164 + (-Ifges(7,5) * t130 - Ifges(7,6) * t131) * t252 + (-Ifges(7,1) * t130 - Ifges(7,4) * t131) * t256 + (-Ifges(7,4) * t130 - Ifges(7,2) * t131) * t258 + (-t1 * t187 - t143 * t2 + t217 * t5 - t218 * t6) * mrSges(7,3) + t30 * (mrSges(7,1) * t187 + t143 * mrSges(7,2)) + (Ifges(7,4) * t143 - Ifges(7,2) * t187) * t260 + (Ifges(7,1) * t143 - Ifges(7,4) * t187) * t261 - t230 * t73 + t101 * t41 + t102 * t42 - t44 * t111 - t43 * t112 - pkin(4) * t115 - t57 * t40 - t37 * mrSges(5,2) - t38 * mrSges(5,1) + t143 * t14 / 0.2e1 - t72 * t155 + t162 * t17; -t140 * t111 + t141 * t112 - t198 * t68 + t89 * t69 + t231 + (-t198 * t6 + t5 * t89 + t30) * m(7) + (-t140 * t25 + t141 * t24 + t38) * m(6); t157 - t45 * (mrSges(7,1) * t89 + mrSges(7,2) * t198) + (Ifges(7,1) * t198 - t249) * t257 + t33 * t256 + (Ifges(7,5) * t198 - Ifges(7,6) * t89) * t253 - t5 * t68 + t6 * t69 + (t198 * t5 + t6 * t89) * mrSges(7,3) + (-Ifges(7,2) * t89 + t34 + t83) * t259 - t263;];
tauc  = t7(:);
