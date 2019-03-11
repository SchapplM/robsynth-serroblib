% Calculate time derivative of joint inertia matrix for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR11_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR11_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:29:17
% EndTime: 2019-03-09 14:29:27
% DurationCPUTime: 4.59s
% Computational Cost: add. (8076->480), mult. (16877->714), div. (0->0), fcn. (14746->8), ass. (0->214)
t182 = sin(qJ(5));
t183 = sin(qJ(4));
t230 = qJD(5) * t182;
t233 = qJD(4) * t183;
t186 = cos(qJ(5));
t187 = cos(qJ(4));
t239 = t186 * t187;
t274 = qJD(4) + qJD(5);
t103 = -t182 * t233 - t183 * t230 + t239 * t274;
t201 = t182 * t187 + t186 * t183;
t104 = t274 * t201;
t181 = sin(qJ(6));
t185 = cos(qJ(6));
t200 = t182 * t183 - t239;
t276 = -t181 * t201 - t185 * t200;
t34 = qJD(6) * t276 + t103 * t185 - t104 * t181;
t95 = t181 * t200 - t185 * t201;
t286 = t34 * t95;
t188 = cos(qJ(2));
t231 = qJD(4) * t188;
t223 = t183 * t231;
t184 = sin(qJ(2));
t235 = qJD(2) * t184;
t196 = t187 * t235 + t223;
t121 = t201 * t188;
t259 = pkin(3) + pkin(7);
t160 = t259 * t184;
t146 = t187 * t160;
t189 = -pkin(2) - pkin(8);
t245 = qJ(3) * t184;
t135 = t188 * t189 - pkin(1) - t245;
t218 = pkin(9) * t188 - t135;
t81 = pkin(4) * t184 + t183 * t218 + t146;
t145 = t183 * t160;
t102 = t187 * t135 + t145;
t236 = t187 * t188;
t87 = -pkin(9) * t236 + t102;
t51 = -t182 * t87 + t186 * t81;
t29 = pkin(5) * t184 + pkin(10) * t121 + t51;
t120 = t200 * t188;
t52 = t182 * t81 + t186 * t87;
t30 = pkin(10) * t120 + t52;
t15 = -t181 * t30 + t185 * t29;
t16 = t181 * t29 + t185 * t30;
t191 = qJD(6) * t95 - t103 * t181 - t185 * t104;
t215 = pkin(2) * t235 - qJD(3) * t184;
t116 = (pkin(8) * t184 - qJ(3) * t188) * qJD(2) + t215;
t234 = qJD(2) * t188;
t153 = t259 * t234;
t217 = -t116 * t183 + t187 * t153;
t40 = (-pkin(9) * t183 * t184 + pkin(4) * t188) * qJD(2) + (t187 * t218 - t145) * qJD(4) + t217;
t232 = qJD(4) * t187;
t53 = t187 * t116 - t135 * t233 + t183 * t153 + t160 * t232;
t48 = pkin(9) * t196 + t53;
t14 = -qJD(5) * t52 - t182 * t48 + t186 * t40;
t271 = t188 * t274;
t68 = t200 * t271 + t201 * t235;
t6 = pkin(5) * t234 - pkin(10) * t68 + t14;
t229 = qJD(5) * t186;
t13 = t182 * t40 + t186 * t48 + t81 * t229 - t230 * t87;
t69 = -t200 * t235 + t201 * t271;
t8 = pkin(10) * t69 + t13;
t2 = qJD(6) * t15 + t181 * t6 + t185 * t8;
t3 = -qJD(6) * t16 - t181 * t8 + t185 * t6;
t290 = t15 * t191 + t16 * t34 - t2 * t95 + t276 * t3;
t253 = pkin(9) - t189;
t141 = t253 * t233;
t155 = t253 * t187;
t142 = qJD(4) * t155;
t154 = t253 * t183;
t61 = t182 * t141 - t186 * t142 + t154 * t230 - t155 * t229;
t41 = -pkin(10) * t103 + t61;
t109 = -t186 * t154 - t182 * t155;
t62 = -qJD(5) * t109 + t186 * t141 + t142 * t182;
t42 = pkin(10) * t104 + t62;
t108 = t154 * t182 - t186 * t155;
t79 = pkin(10) * t200 + t108;
t80 = -pkin(10) * t201 + t109;
t46 = -t181 * t80 + t185 * t79;
t10 = qJD(6) * t46 + t181 * t42 + t185 * t41;
t47 = t181 * t79 + t185 * t80;
t11 = -qJD(6) * t47 - t181 * t41 + t185 * t42;
t289 = -t10 * t95 + t11 * t276 + t191 * t46 + t34 * t47;
t172 = pkin(4) * t186 + pkin(5);
t243 = t181 * t182;
t127 = -pkin(4) * t243 + t172 * t185;
t242 = t182 * t185;
t128 = pkin(4) * t242 + t172 * t181;
t227 = qJD(6) * t185;
t228 = qJD(6) * t181;
t89 = t172 * t227 + (-t182 * t228 + (t185 * t186 - t243) * qJD(5)) * pkin(4);
t90 = -t172 * t228 + (-t182 * t227 + (-t181 * t186 - t242) * qJD(5)) * pkin(4);
t288 = t127 * t191 + t128 * t34 + t276 * t90 - t89 * t95;
t287 = qJD(6) * (t181 * t276 + t185 * t95) - t181 * t34 - t185 * t191;
t281 = 0.2e1 * qJ(3);
t280 = m(4) * pkin(7);
t277 = t191 * t276;
t54 = -qJD(4) * t102 + t217;
t275 = t183 * t53 + t187 * t54;
t273 = Ifges(6,5) * t68 + Ifges(6,6) * t69 + Ifges(6,3) * t234;
t272 = qJD(5) * (-t182 * t200 - t186 * t201) - t103 * t182 + t104 * t186;
t224 = t183 * t235;
t270 = Ifges(5,5) * t224 + Ifges(5,6) * t196 + Ifges(5,3) * t234;
t266 = 2 * m(6);
t265 = 2 * m(7);
t264 = -0.2e1 * pkin(1);
t123 = -qJ(3) * t234 + t215;
t263 = 0.2e1 * t123;
t209 = -pkin(2) * t188 - t245;
t156 = -pkin(1) + t209;
t262 = -0.2e1 * t156;
t261 = t95 / 0.2e1;
t260 = t276 / 0.2e1;
t258 = -t201 / 0.2e1;
t257 = -t200 / 0.2e1;
t255 = -t188 / 0.2e1;
t179 = t188 * pkin(7);
t254 = t89 * mrSges(7,2);
t252 = Ifges(7,5) * t191 - Ifges(7,6) * t34;
t251 = Ifges(5,4) * t183;
t250 = Ifges(5,4) * t187;
t248 = t184 * Ifges(5,5);
t246 = -Ifges(6,5) * t104 - Ifges(6,6) * t103;
t212 = Ifges(5,1) * t183 + t250;
t119 = -t188 * t212 + t248;
t241 = t183 * t119;
t159 = Ifges(5,1) * t187 - t251;
t240 = t183 * t159;
t211 = Ifges(5,2) * t187 + t251;
t118 = t184 * Ifges(5,6) - t188 * t211;
t238 = t187 * t118;
t158 = -Ifges(5,2) * t183 + t250;
t237 = t187 * t158;
t169 = t183 * pkin(4) + qJ(3);
t161 = t188 * pkin(3) + t179;
t165 = pkin(4) * t232 + qJD(3);
t226 = 2 * mrSges(6,3);
t77 = t120 * t185 + t121 * t181;
t24 = qJD(6) * t77 + t181 * t69 + t185 * t68;
t78 = t120 * t181 - t121 * t185;
t25 = -qJD(6) * t78 - t181 * t68 + t185 * t69;
t225 = Ifges(7,5) * t24 + Ifges(7,6) * t25 + Ifges(7,3) * t234;
t129 = pkin(4) * t236 + t161;
t222 = t187 * t231;
t220 = mrSges(7,1) * t191 - t34 * mrSges(7,2);
t88 = t90 * mrSges(7,1);
t219 = t88 - t254;
t214 = t11 * mrSges(7,1) - t10 * mrSges(7,2) + t252;
t213 = mrSges(5,1) * t187 - mrSges(5,2) * t183;
t157 = mrSges(5,1) * t183 + mrSges(5,2) * t187;
t210 = -Ifges(5,5) * t183 - Ifges(5,6) * t187;
t101 = -t135 * t183 + t146;
t208 = t101 * t183 - t102 * t187;
t207 = -t103 * t201 - t104 * t200;
t150 = mrSges(5,3) * t183 * t188 + mrSges(5,1) * t184;
t151 = -mrSges(5,2) * t184 - mrSges(5,3) * t236;
t202 = -t183 * t150 + t187 * t151;
t199 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t225;
t198 = -t104 * mrSges(6,1) - t103 * mrSges(6,2) + t220;
t197 = (-mrSges(6,1) * t182 - mrSges(6,2) * t186) * qJD(5) * pkin(4);
t195 = -t222 + t224;
t194 = t103 * t52 - t104 * t51 + t13 * t201 - t14 * t200;
t193 = t62 * mrSges(6,1) - t61 * mrSges(6,2) + t214 + t246;
t192 = t103 * t109 - t104 * t108 - t200 * t62 + t201 * t61;
t110 = -pkin(4) * t223 + (-pkin(4) * t187 - t259) * t235;
t190 = t14 * mrSges(6,1) - t13 * mrSges(6,2) + t199 + t273;
t152 = t259 * t235;
t149 = t212 * qJD(4);
t148 = t211 * qJD(4);
t147 = t213 * qJD(4);
t136 = (-mrSges(7,1) * t181 - mrSges(7,2) * t185) * qJD(6) * pkin(5);
t134 = t213 * t188;
t117 = pkin(5) * t201 + t169;
t115 = mrSges(5,1) * t234 - mrSges(5,3) * t195;
t114 = -mrSges(5,2) * t234 + mrSges(5,3) * t196;
t113 = mrSges(6,1) * t184 + mrSges(6,3) * t121;
t112 = -mrSges(6,2) * t184 + mrSges(6,3) * t120;
t107 = -Ifges(6,1) * t200 - Ifges(6,4) * t201;
t106 = -Ifges(6,4) * t200 - Ifges(6,2) * t201;
t105 = mrSges(6,1) * t201 - mrSges(6,2) * t200;
t92 = -mrSges(5,1) * t196 + mrSges(5,2) * t195;
t91 = -pkin(5) * t120 + t129;
t86 = pkin(5) * t103 + t165;
t85 = -mrSges(6,1) * t120 - mrSges(6,2) * t121;
t84 = -t159 * t231 + (t188 * Ifges(5,5) + t184 * t212) * qJD(2);
t83 = -t158 * t231 + (t188 * Ifges(5,6) + t184 * t211) * qJD(2);
t76 = -Ifges(6,1) * t121 + Ifges(6,4) * t120 + Ifges(6,5) * t184;
t75 = -Ifges(6,4) * t121 + Ifges(6,2) * t120 + Ifges(6,6) * t184;
t71 = mrSges(7,1) * t184 - mrSges(7,3) * t78;
t70 = -mrSges(7,2) * t184 + mrSges(7,3) * t77;
t65 = -Ifges(6,1) * t104 - Ifges(6,4) * t103;
t64 = -Ifges(6,4) * t104 - Ifges(6,2) * t103;
t63 = mrSges(6,1) * t103 - mrSges(6,2) * t104;
t59 = -mrSges(6,2) * t234 + mrSges(6,3) * t69;
t58 = mrSges(6,1) * t234 - mrSges(6,3) * t68;
t57 = Ifges(7,1) * t276 + Ifges(7,4) * t95;
t56 = Ifges(7,4) * t276 + Ifges(7,2) * t95;
t55 = -mrSges(7,1) * t95 + mrSges(7,2) * t276;
t50 = -mrSges(7,1) * t77 + mrSges(7,2) * t78;
t49 = -pkin(5) * t69 + t110;
t44 = Ifges(7,1) * t78 + Ifges(7,4) * t77 + Ifges(7,5) * t184;
t43 = Ifges(7,4) * t78 + Ifges(7,2) * t77 + Ifges(7,6) * t184;
t28 = -mrSges(6,1) * t69 + mrSges(6,2) * t68;
t27 = Ifges(6,1) * t68 + Ifges(6,4) * t69 + Ifges(6,5) * t234;
t26 = Ifges(6,4) * t68 + Ifges(6,2) * t69 + Ifges(6,6) * t234;
t21 = -mrSges(7,2) * t234 + mrSges(7,3) * t25;
t20 = mrSges(7,1) * t234 - mrSges(7,3) * t24;
t19 = Ifges(7,1) * t191 - Ifges(7,4) * t34;
t18 = Ifges(7,4) * t191 - Ifges(7,2) * t34;
t17 = mrSges(7,1) * t34 + mrSges(7,2) * t191;
t7 = -mrSges(7,1) * t25 + mrSges(7,2) * t24;
t5 = Ifges(7,1) * t24 + Ifges(7,4) * t25 + Ifges(7,5) * t234;
t4 = Ifges(7,4) * t24 + Ifges(7,2) * t25 + Ifges(7,6) * t234;
t1 = [((mrSges(3,1) * t264 + mrSges(4,2) * t262 + t238 + t241 + 0.2e1 * (-Ifges(4,6) - Ifges(3,4)) * t184) * t184 + (mrSges(3,2) * t264 + mrSges(4,3) * t262 - Ifges(6,5) * t121 + Ifges(7,5) * t78 + Ifges(6,6) * t120 + Ifges(7,6) * t77 + (0.2e1 * Ifges(3,4) + 0.2e1 * Ifges(4,6) + t210) * t188 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + Ifges(5,3) + Ifges(6,3) + Ifges(7,3)) * t184) * t188) * qJD(2) + 0.2e1 * m(5) * (t101 * t54 + t102 * t53 - t152 * t161) + (mrSges(4,2) * t263 - t183 * t84 - t187 * t83 + (t118 * t183 + (-t119 - t248) * t187) * qJD(4)) * t188 + (-0.2e1 * t123 * mrSges(4,3) + t225 + t270 + t273) * t184 + m(4) * t156 * t263 + 0.2e1 * t161 * t92 + 0.2e1 * t54 * t150 + 0.2e1 * t53 * t151 - 0.2e1 * t152 * t134 + t120 * t26 - t121 * t27 + 0.2e1 * t129 * t28 + (t15 * t3 + t16 * t2 + t49 * t91) * t265 + (t110 * t129 + t13 * t52 + t14 * t51) * t266 + 0.2e1 * t15 * t20 + 0.2e1 * t16 * t21 + t25 * t43 + t24 * t44 + 0.2e1 * t49 * t50 + 0.2e1 * t51 * t58 + 0.2e1 * t52 * t59 + 0.2e1 * t2 * t70 + 0.2e1 * t3 * t71 + t69 * t75 + t68 * t76 + t77 * t4 + t78 * t5 + 0.2e1 * t91 * t7 + 0.2e1 * t110 * t85 + 0.2e1 * t13 * t112 + 0.2e1 * t14 * t113 + 0.2e1 * t102 * t114 + 0.2e1 * t101 * t115; -t290 * mrSges(7,3) + (-t148 * t255 + t189 * t115 - t54 * mrSges(5,3) + t84 / 0.2e1) * t187 + (-t149 * t255 + t189 * t114 - t53 * mrSges(5,3) - t83 / 0.2e1) * t183 + t191 * t44 / 0.2e1 - t194 * mrSges(6,3) + m(5) * (-qJ(3) * t152 + t189 * t275) + (-t238 / 0.2e1 - t241 / 0.2e1 + (-t187 * t159 / 0.2e1 + t183 * t158 / 0.2e1) * t188 + t208 * mrSges(5,3) + (-m(5) * t208 + t202) * t189) * qJD(4) - t34 * t43 / 0.2e1 + (t210 * qJD(4) + t246 + t252) * t184 / 0.2e1 + (m(4) * t179 + m(5) * t161 + t188 * mrSges(4,1) + t134) * qJD(3) + m(6) * (t108 * t14 + t109 * t13 + t110 * t169 + t129 * t165 + t51 * t62 + t52 * t61) + m(7) * (t10 * t16 + t11 * t15 + t117 * t49 + t2 * t47 + t3 * t46 + t86 * t91) + t161 * t147 + t165 * t85 + t169 * t28 - t152 * t157 + t120 * t64 / 0.2e1 - t121 * t65 / 0.2e1 + t129 * t63 + (t209 * t280 + (-pkin(2) * mrSges(4,1) + Ifges(5,5) * t187 / 0.2e1 - Ifges(5,6) * t183 / 0.2e1 + Ifges(6,5) * t257 + Ifges(6,6) * t258 + Ifges(7,5) * t260 + Ifges(7,6) * t261 + Ifges(3,5) - Ifges(4,4) + (-mrSges(3,1) + mrSges(4,2)) * pkin(7)) * t188 + (t237 / 0.2e1 + t240 / 0.2e1 - qJ(3) * mrSges(4,1) - Ifges(3,6) + Ifges(4,5) + (mrSges(3,2) - mrSges(4,3)) * pkin(7)) * t184) * qJD(2) + t27 * t257 + t26 * t258 + t5 * t260 + t4 * t261 + t46 * t20 + t47 * t21 + t49 * t55 + t25 * t56 / 0.2e1 + t24 * t57 / 0.2e1 + t10 * t70 + t11 * t71 + t77 * t18 / 0.2e1 + t78 * t19 / 0.2e1 + t86 * t50 + t91 * t17 + qJ(3) * t92 - t103 * t75 / 0.2e1 - t104 * t76 / 0.2e1 + t69 * t106 / 0.2e1 + t68 * t107 / 0.2e1 + t108 * t58 + t109 * t59 + t110 * t105 + t61 * t112 + t62 * t113 + t117 * t7; t147 * t281 - t103 * t106 - t104 * t107 + 0.2e1 * t165 * t105 + 0.2e1 * t117 * t17 - t201 * t64 - t200 * t65 + t183 * t148 - t187 * t149 + 0.2e1 * t169 * t63 + t95 * t18 + t276 * t19 + t191 * t57 - t34 * t56 + 0.2e1 * t86 * t55 + (-t237 - t240) * qJD(4) + (t108 * t62 + t109 * t61 + t165 * t169) * t266 + (t10 * t47 + t11 * t46 + t117 * t86) * t265 - t192 * t226 - 0.2e1 * t289 * mrSges(7,3) + (0.2e1 * mrSges(4,3) + 0.2e1 * t157 + (m(4) + m(5)) * t281) * qJD(3); t103 * t112 - t104 * t113 + t183 * t114 + t187 * t115 + t201 * t59 - t200 * t58 + t276 * t20 - t95 * t21 + t34 * t70 + t191 * t71 + t202 * qJD(4) + (mrSges(4,1) + t280) * t234 + m(7) * t290 + m(6) * t194 + m(5) * (-t208 * qJD(4) + t275); t207 * t226 + m(7) * t289 + m(6) * t192 + (-0.2e1 * t277 + 0.2e1 * t286) * mrSges(7,3); -0.2e1 * m(6) * t207 + 0.2e1 * m(7) * (t277 - t286); -Ifges(5,5) * t222 + t190 + m(7) * (t127 * t3 + t128 * t2 + t15 * t90 + t16 * t89) + (t112 * t229 - t113 * t230 + m(6) * (t13 * t182 + t14 * t186 + t229 * t52 - t230 * t51) + t186 * t58 + t182 * t59) * pkin(4) + t127 * t20 + t128 * t21 - t53 * mrSges(5,2) + t54 * mrSges(5,1) + t89 * t70 + t90 * t71 + t270; m(7) * (t10 * t128 + t11 * t127 + t46 * t90 + t47 * t89) + ((-mrSges(5,2) * t189 - Ifges(5,6)) * t187 + (-mrSges(5,1) * t189 - Ifges(5,5)) * t183) * qJD(4) - t288 * mrSges(7,3) + (m(6) * (t182 * t61 + t186 * t62 + (-t108 * t182 + t109 * t186) * qJD(5)) + t272 * mrSges(6,3)) * pkin(4) + t193; -m(6) * t272 * pkin(4) + m(7) * t288 - t157 * qJD(4) + t198; 0.2e1 * t88 + (t127 * t90 + t128 * t89) * t265 - 0.2e1 * t254 + 0.2e1 * t197; (m(7) * (-t15 * t228 + t16 * t227 + t181 * t2 + t185 * t3) + t70 * t227 + t181 * t21 - t71 * t228 + t185 * t20) * pkin(5) + t190; (m(7) * (t10 * t181 + t11 * t185 + (-t181 * t46 + t185 * t47) * qJD(6)) + t287 * mrSges(7,3)) * pkin(5) + t193; -m(7) * pkin(5) * t287 + t198; t197 + (-mrSges(7,1) * t228 + m(7) * (-t127 * t228 + t128 * t227 + t181 * t89 + t185 * t90) - mrSges(7,2) * t227) * pkin(5) + t219; 0.2e1 * t136; t199; t214; t220; t219; t136; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
