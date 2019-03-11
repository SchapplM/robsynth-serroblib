% Calculate time derivative of joint inertia matrix for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:04:28
% EndTime: 2019-03-09 11:04:42
% DurationCPUTime: 5.96s
% Computational Cost: add. (7314->570), mult. (19045->819), div. (0->0), fcn. (19001->10), ass. (0->229)
t300 = Ifges(6,4) - Ifges(5,5);
t299 = Ifges(6,5) - Ifges(5,6);
t197 = sin(pkin(11));
t199 = cos(pkin(11));
t202 = sin(qJ(4));
t253 = t202 * t199;
t276 = cos(qJ(4));
t167 = t197 * t276 + t253;
t161 = t167 * qJD(4);
t166 = t197 * t202 - t199 * t276;
t201 = sin(qJ(6));
t204 = cos(qJ(6));
t244 = qJD(6) * t204;
t212 = t201 * t161 + t166 * t244;
t245 = qJD(6) * t201;
t211 = -t204 * t161 + t166 * t245;
t298 = 2 * mrSges(6,1) + 2 * mrSges(5,3);
t297 = -t204 / 0.2e1;
t235 = qJD(4) * t276;
t246 = qJD(4) * t202;
t160 = t197 * t246 - t199 * t235;
t296 = t160 * t300 + t299 * t161;
t200 = cos(pkin(6));
t198 = sin(pkin(6));
t203 = sin(qJ(2));
t258 = t198 * t203;
t159 = t197 * t200 + t199 * t258;
t240 = t197 * t258;
t256 = t199 * t200;
t214 = t240 - t256;
t107 = t159 * t276 - t202 * t214;
t205 = cos(qJ(2));
t257 = t198 * t205;
t163 = t200 * t203 * pkin(1) + pkin(8) * t257;
t151 = t163 * qJD(2);
t248 = qJD(2) * t198;
t238 = t205 * t248;
t229 = t197 * t238;
t128 = pkin(3) * t229 + t151;
t209 = t276 * t214;
t228 = t199 * t238;
t76 = qJD(4) * t209 - t276 * t228 + (qJD(4) * t159 + t229) * t202;
t208 = qJ(5) * t76 - qJD(5) * t107 + t128;
t286 = pkin(4) + pkin(10);
t77 = qJD(4) * t107 + t167 * t238;
t14 = t286 * t77 + t208;
t255 = t199 * t205;
t134 = (-qJD(3) * t203 + (pkin(2) * t203 - qJ(3) * t205) * qJD(2)) * t198;
t247 = qJD(2) * t203;
t239 = t198 * t247;
t275 = pkin(1) * t205;
t243 = t200 * t275;
t150 = -pkin(8) * t239 + qJD(2) * t243;
t139 = qJD(3) * t200 + t150;
t87 = t199 * t134 - t139 * t197;
t63 = (pkin(3) * t203 - pkin(9) * t255) * t248 + t87;
t145 = qJ(3) * t200 + t163;
t146 = (-pkin(2) * t205 - qJ(3) * t203 - pkin(1)) * t198;
t100 = -t145 * t197 + t199 * t146;
t65 = -pkin(3) * t257 - pkin(9) * t159 + t100;
t88 = t197 * t134 + t199 * t139;
t81 = -pkin(9) * t229 + t88;
t101 = t199 * t145 + t197 * t146;
t82 = -pkin(9) * t214 + t101;
t216 = t202 * t81 + t82 * t235 + t65 * t246 - t276 * t63;
t3 = -t76 * pkin(5) - t239 * t286 + t216;
t33 = -t202 * t82 + t276 * t65;
t25 = pkin(4) * t257 - t33;
t21 = t107 * pkin(5) + pkin(10) * t257 + t25;
t106 = t159 * t202 + t209;
t188 = pkin(8) * t258;
t193 = -pkin(3) * t199 - pkin(2);
t108 = pkin(3) * t240 + t188 + (t193 - t275) * t200;
t207 = -t107 * qJ(5) + t108;
t23 = t106 * t286 + t207;
t5 = -t201 * t23 + t204 * t21;
t1 = qJD(6) * t5 + t14 * t204 + t201 * t3;
t6 = t201 * t21 + t204 * t23;
t2 = -qJD(6) * t6 - t14 * t201 + t204 * t3;
t295 = t1 * t201 + t2 * t204;
t213 = -t202 * t63 - t65 * t235 + t246 * t82 - t276 * t81;
t10 = -t198 * (qJ(5) * t247 - qJD(5) * t205) + t213;
t270 = pkin(9) + qJ(3);
t234 = t270 * t197;
t165 = t276 * t234;
t176 = t270 * t199;
t233 = t276 * qJD(3);
t102 = (qJD(3) * t197 + qJD(4) * t176) * t202 + qJD(4) * t165 - t199 * t233;
t294 = 2 * m(4);
t293 = 2 * m(5);
t292 = 2 * m(6);
t291 = 2 * m(7);
t290 = -0.2e1 * pkin(1);
t289 = -2 * mrSges(3,3);
t215 = -t106 * t201 + t204 * t257;
t39 = qJD(6) * t215 - t201 * t239 + t204 * t77;
t89 = t106 * t204 + t201 * t257;
t40 = qJD(6) * t89 + t201 * t77 + t204 * t239;
t9 = Ifges(7,1) * t40 + Ifges(7,4) * t39 - Ifges(7,5) * t76;
t288 = t9 / 0.2e1;
t49 = Ifges(7,1) * t212 - Ifges(7,4) * t211 - Ifges(7,5) * t160;
t287 = t49 / 0.2e1;
t223 = Ifges(7,5) * t201 + Ifges(7,6) * t204;
t285 = -t223 * qJD(6) / 0.2e1;
t266 = Ifges(7,4) * t201;
t224 = Ifges(7,2) * t204 + t266;
t171 = t224 * qJD(6);
t284 = -t171 / 0.2e1;
t265 = Ifges(7,4) * t204;
t226 = Ifges(7,1) * t201 + t265;
t172 = t226 * qJD(6);
t283 = -t172 / 0.2e1;
t282 = Ifges(7,5) * t297 + Ifges(7,6) * t201 / 0.2e1;
t179 = -Ifges(7,2) * t201 + t265;
t281 = t179 / 0.2e1;
t180 = Ifges(7,1) * t204 - t266;
t280 = t180 / 0.2e1;
t279 = -t201 / 0.2e1;
t278 = t204 / 0.2e1;
t271 = mrSges(6,2) - mrSges(5,1);
t34 = t202 * t65 + t276 * t82;
t269 = mrSges(7,3) * t166;
t268 = Ifges(4,4) * t197;
t267 = Ifges(4,4) * t199;
t264 = Ifges(4,6) * t197;
t263 = Ifges(4,6) * t199;
t262 = t128 * mrSges(5,1);
t261 = t128 * mrSges(5,2);
t129 = t176 * t202 + t165;
t104 = pkin(5) * t167 + t129;
t218 = -qJ(5) * t167 + t193;
t92 = t166 * t286 + t218;
t43 = t104 * t204 - t201 * t92;
t260 = qJD(6) * t43;
t44 = t104 * t201 + t204 * t92;
t259 = qJD(6) * t44;
t135 = mrSges(4,1) * t229 + mrSges(4,2) * t228;
t7 = Ifges(7,5) * t40 + Ifges(7,6) * t39 - Ifges(7,3) * t76;
t242 = Ifges(6,1) * t239 + Ifges(6,4) * t76 + Ifges(6,5) * t77;
t241 = -Ifges(5,5) * t76 - Ifges(5,6) * t77 + Ifges(5,3) * t239;
t56 = -t76 * mrSges(6,1) + mrSges(6,2) * t239;
t227 = t201 * t5 - t204 * t6;
t177 = mrSges(7,1) * t201 + mrSges(7,2) * t204;
t225 = -t197 * Ifges(4,2) + t267;
t45 = -mrSges(7,2) * t107 + mrSges(7,3) * t89;
t46 = mrSges(7,1) * t107 + mrSges(7,3) * t215;
t222 = -t201 * t46 + t204 * t45;
t103 = t176 * t235 + qJD(3) * t253 + (-t246 * t270 + t233) * t197;
t130 = t176 * t276 - t202 * t234;
t221 = -t102 * t130 + t103 * t129;
t220 = qJ(5) * t160 - qJD(5) * t167;
t24 = qJ(5) * t257 - t34;
t27 = -Ifges(7,4) * t215 + Ifges(7,2) * t89 + Ifges(7,6) * t107;
t28 = -Ifges(7,1) * t215 + Ifges(7,4) * t89 + Ifges(7,5) * t107;
t217 = t27 * t297 + t28 * t279;
t47 = t212 * Ifges(7,5) - t211 * Ifges(7,6) - Ifges(7,3) * t160;
t183 = Ifges(3,5) * t238;
t169 = -mrSges(7,1) * t244 + mrSges(7,2) * t245;
t162 = -t188 + t243;
t153 = t160 * mrSges(5,2);
t152 = t160 * mrSges(6,3);
t148 = t188 + (-pkin(2) - t275) * t200;
t142 = (mrSges(4,1) * t203 - mrSges(4,3) * t255) * t248;
t141 = (-mrSges(4,3) * t197 * t205 - mrSges(4,2) * t203) * t248;
t133 = -mrSges(4,1) * t257 - mrSges(4,3) * t159;
t132 = mrSges(4,2) * t257 - mrSges(4,3) * t214;
t125 = Ifges(5,1) * t167 - Ifges(5,4) * t166;
t124 = Ifges(5,4) * t167 - Ifges(5,2) * t166;
t123 = -Ifges(6,2) * t167 + Ifges(6,6) * t166;
t122 = -Ifges(6,6) * t167 + Ifges(6,3) * t166;
t121 = -mrSges(6,2) * t166 - mrSges(6,3) * t167;
t120 = (t203 * Ifges(4,5) + (t199 * Ifges(4,1) - t268) * t205) * t248;
t119 = (t203 * Ifges(4,6) + t205 * t225) * t248;
t118 = -mrSges(7,2) * t167 + t204 * t269;
t117 = mrSges(7,1) * t167 - t201 * t269;
t116 = (-mrSges(7,1) * t204 + mrSges(7,2) * t201) * t166;
t115 = pkin(4) * t166 + t218;
t114 = -Ifges(5,1) * t160 - Ifges(5,4) * t161;
t113 = -Ifges(5,4) * t160 - Ifges(5,2) * t161;
t112 = Ifges(6,2) * t160 + Ifges(6,6) * t161;
t111 = Ifges(6,6) * t160 + Ifges(6,3) * t161;
t110 = t161 * mrSges(5,1) - t153;
t109 = -t161 * mrSges(6,2) + t152;
t105 = -t166 * pkin(5) + t130;
t99 = -mrSges(5,1) * t257 - mrSges(5,3) * t107;
t98 = mrSges(5,2) * t257 - mrSges(5,3) * t106;
t97 = mrSges(6,1) * t107 - mrSges(6,2) * t257;
t96 = mrSges(6,1) * t106 + mrSges(6,3) * t257;
t95 = Ifges(7,5) * t167 + t166 * t226;
t94 = Ifges(7,6) * t167 + t166 * t224;
t93 = Ifges(7,3) * t167 + t166 * t223;
t91 = pkin(4) * t161 + t220;
t86 = -mrSges(7,1) * t160 - mrSges(7,3) * t212;
t85 = mrSges(7,2) * t160 - mrSges(7,3) * t211;
t84 = -t160 * pkin(5) + t103;
t83 = -pkin(5) * t161 - t102;
t78 = t161 * t286 + t220;
t69 = t76 * mrSges(5,2);
t68 = t76 * mrSges(6,3);
t60 = mrSges(7,1) * t211 + mrSges(7,2) * t212;
t58 = -mrSges(5,2) * t239 - mrSges(5,3) * t77;
t57 = mrSges(5,1) * t239 + mrSges(5,3) * t76;
t55 = mrSges(6,1) * t77 - mrSges(6,3) * t239;
t54 = -mrSges(6,2) * t106 - mrSges(6,3) * t107;
t53 = Ifges(5,1) * t107 - Ifges(5,4) * t106 - Ifges(5,5) * t257;
t52 = Ifges(5,4) * t107 - Ifges(5,2) * t106 - Ifges(5,6) * t257;
t51 = -Ifges(6,4) * t257 - Ifges(6,2) * t107 + Ifges(6,6) * t106;
t50 = -Ifges(6,5) * t257 - Ifges(6,6) * t107 + Ifges(6,3) * t106;
t48 = Ifges(7,4) * t212 - Ifges(7,2) * t211 - Ifges(7,6) * t160;
t42 = -mrSges(7,1) * t89 - mrSges(7,2) * t215;
t41 = t106 * pkin(4) + t207;
t36 = t77 * mrSges(5,1) - t69;
t35 = -t77 * mrSges(6,2) + t68;
t32 = -Ifges(5,1) * t76 - Ifges(5,4) * t77 + Ifges(5,5) * t239;
t31 = -Ifges(5,4) * t76 - Ifges(5,2) * t77 + Ifges(5,6) * t239;
t30 = Ifges(6,4) * t239 + Ifges(6,2) * t76 + Ifges(6,6) * t77;
t29 = Ifges(6,5) * t239 + Ifges(6,6) * t76 + Ifges(6,3) * t77;
t26 = -Ifges(7,5) * t215 + Ifges(7,6) * t89 + Ifges(7,3) * t107;
t22 = -pkin(5) * t106 - t24;
t20 = -mrSges(7,1) * t76 - mrSges(7,3) * t40;
t19 = mrSges(7,2) * t76 + mrSges(7,3) * t39;
t18 = pkin(4) * t77 + t208;
t17 = -t201 * t78 + t204 * t84 - t259;
t16 = t201 * t84 + t204 * t78 + t260;
t15 = -mrSges(7,1) * t39 + mrSges(7,2) * t40;
t11 = -pkin(4) * t239 + t216;
t8 = Ifges(7,4) * t40 + Ifges(7,2) * t39 - Ifges(7,6) * t76;
t4 = -pkin(5) * t77 - t10;
t12 = [((t163 * t289 + Ifges(4,5) * t159 + (-(2 * Ifges(3,6)) + t263) * t200 + (mrSges(3,1) * t290 + (-(2 * Ifges(3,4)) - t264) * t203) * t198 - t300 * t107 + t299 * t106) * t203 + (Ifges(3,5) * t200 + t162 * t289 - t197 * (Ifges(4,4) * t159 + Ifges(4,2) * t256) + t199 * (Ifges(4,1) * t159 + Ifges(4,4) * t256) + (mrSges(3,2) * t290 + 0.2e1 * (-Ifges(4,5) * t199 + Ifges(3,4) + t264) * t205) * t198 + (-t197 * t225 + (2 * Ifges(3,1)) - Ifges(6,1) - (2 * Ifges(3,2)) - (2 * Ifges(4,3)) - Ifges(5,3)) * t258) * t205) * t248 + t200 * t183 + (-t30 + t32 + t7 + 0.2e1 * t261) * t107 + (t29 - t31 + 0.2e1 * t262) * t106 - 0.2e1 * t151 * (mrSges(3,1) * t200 - mrSges(3,3) * t258) - t242 * t257 - t241 * t257 + 0.2e1 * t150 * (-mrSges(3,2) * t200 + mrSges(3,3) * t257) - t214 * t119 + 0.2e1 * t151 * (mrSges(4,1) * t214 + t159 * mrSges(4,2)) + (t100 * t87 + t101 * t88 + t148 * t151) * t294 + (t1 * t6 + t2 * t5 + t22 * t4) * t291 + (t10 * t24 + t11 * t25 + t18 * t41) * t292 - 0.2e1 * t213 * t98 + (t108 * t128 - t213 * t34 - t216 * t33) * t293 - 0.2e1 * t216 * t99 - t215 * t9 + (t51 - t26 - t53) * t76 + (t50 - t52) * t77 + 0.2e1 * t6 * t19 + 0.2e1 * t5 * t20 + 0.2e1 * t22 * t15 + t39 * t27 + t40 * t28 + 0.2e1 * t41 * t35 + 0.2e1 * t4 * t42 + 0.2e1 * t1 * t45 + 0.2e1 * t2 * t46 + 0.2e1 * t18 * t54 + 0.2e1 * t24 * t55 + 0.2e1 * t25 * t56 + 0.2e1 * t33 * t57 + 0.2e1 * t34 * t58 + t89 * t8 + 0.2e1 * t10 * t96 + 0.2e1 * t11 * t97 + 0.2e1 * t108 * t36 + 0.2e1 * t88 * t132 + 0.2e1 * t87 * t133 + 0.2e1 * m(3) * (t150 * t163 - t151 * t162) + 0.2e1 * t101 * t141 + 0.2e1 * t100 * t142 + 0.2e1 * t148 * t135 + t159 * t120; m(4) * (-pkin(2) * t151 + (-t100 * t197 + t101 * t199) * qJD(3) + (-t197 * t87 + t199 * t88) * qJ(3)) + t183 + (t33 * mrSges(5,3) - t25 * mrSges(6,1) + t51 / 0.2e1 - t26 / 0.2e1 - t53 / 0.2e1) * t160 + (qJD(3) * t132 + qJ(3) * t141 + t88 * mrSges(4,3) + t119 / 0.2e1 - t151 * mrSges(4,1)) * t199 + (t24 * mrSges(6,1) - t34 * mrSges(5,3) + t50 / 0.2e1 - t52 / 0.2e1 - t217) * t161 + (-qJD(3) * t133 - qJ(3) * t142 - t87 * mrSges(4,3) + t120 / 0.2e1 + t151 * mrSges(4,2)) * t197 + (t58 - t55) * t130 + (t56 - t57) * t129 + (t97 - t99) * t103 + (t96 - t98) * t102 + (t47 / 0.2e1 - t112 / 0.2e1 + t114 / 0.2e1) * t107 + (((-Ifges(3,6) + Ifges(4,5) * t197 / 0.2e1 + t263 / 0.2e1 + (Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1) * t167 + (-Ifges(5,6) / 0.2e1 + Ifges(6,5) / 0.2e1) * t166) * t203 + (t199 * (Ifges(4,1) * t197 + t267) / 0.2e1 - t197 * (Ifges(4,2) * t199 + t268) / 0.2e1) * t205) * qJD(2) - t296 * t205 / 0.2e1) * t198 + (t111 / 0.2e1 - t113 / 0.2e1) * t106 + m(6) * (-t10 * t130 + t102 * t24 + t103 * t25 + t11 * t129 + t115 * t18 + t41 * t91) + m(7) * (t1 * t44 + t105 * t4 + t16 * t6 + t17 * t5 + t2 * t43 + t22 * t83) + (t8 * t278 + t201 * t288 + t10 * mrSges(6,1) + t213 * mrSges(5,3) + t29 / 0.2e1 - t31 / 0.2e1 + t262 + (t27 * t279 + t278 * t28) * qJD(6)) * t166 + (t216 * mrSges(5,3) + t11 * mrSges(6,1) + t7 / 0.2e1 + t32 / 0.2e1 + t261 - t30 / 0.2e1) * t167 + m(5) * (-t102 * t34 - t103 * t33 + t128 * t193 + t129 * t216 - t130 * t213) - t215 * t287 + (t122 / 0.2e1 - t124 / 0.2e1) * t77 + (-t93 / 0.2e1 + t123 / 0.2e1 - t125 / 0.2e1) * t76 + t43 * t20 + t44 * t19 + t16 * t45 + t17 * t46 + t22 * t60 + t83 * t42 + t6 * t85 + t5 * t86 + t89 * t48 / 0.2e1 + t91 * t54 + t39 * t94 / 0.2e1 + t40 * t95 / 0.2e1 + t105 * t15 + t41 * t109 + t108 * t110 + t115 * t35 + t4 * t116 + t2 * t117 + t1 * t118 + t18 * t121 - pkin(2) * t135 - t150 * mrSges(3,2) - t151 * mrSges(3,1) + t193 * t36; 0.2e1 * t105 * t60 + 0.2e1 * t115 * t109 + 0.2e1 * t193 * t110 + 0.2e1 * t83 * t116 + 0.2e1 * t17 * t117 + 0.2e1 * t16 * t118 + 0.2e1 * t91 * t121 + 0.2e1 * t43 * t86 + 0.2e1 * t44 * t85 + (t103 * t298 - t112 + t114 + t47) * t167 + (-t129 * t298 + t123 - t125 - t93) * t160 + t221 * t293 + (t115 * t91 + t221) * t292 + (t105 * t83 + t16 * t44 + t17 * t43) * t291 + (-t130 * t298 + t201 * t95 + t204 * t94 + t122 - t124) * t161 + (t201 * t49 + t204 * t48 + t111 - t113 + t102 * t298 + (-t201 * t94 + t204 * t95) * qJD(6)) * t166 + (qJ(3) * t294 + 0.2e1 * mrSges(4,3)) * qJD(3) * (t197 ^ 2 + t199 ^ 2); t204 * t19 - t201 * t20 + t68 - t69 - t271 * t77 + (-t201 * t45 - t204 * t46) * qJD(6) + m(7) * (t1 * t204 - t2 * t201 + (-t201 * t6 - t204 * t5) * qJD(6)) + m(6) * t18 + m(5) * t128 + m(4) * t151 + t135; -t201 * t86 + t204 * t85 + t152 - t153 - t271 * t161 + (-t117 * t204 - t118 * t201) * qJD(6) + m(7) * (t16 * t204 - t17 * t201 + (-t201 * t44 - t204 * t43) * qJD(6)) + m(6) * t91; 0; (-t2 * mrSges(7,3) - t20 * t286 + t288) * t204 + (-t286 * t19 - t1 * mrSges(7,3) - t8 / 0.2e1) * t201 + (t42 - t96) * qJD(5) + (t15 - t55) * qJ(5) + m(7) * (qJ(5) * t4 + qJD(5) * t22 - t286 * t295) + (t227 * mrSges(7,3) - (-m(7) * t227 + t222) * t286 + t217) * qJD(6) + m(6) * (-pkin(4) * t11 - qJ(5) * t10 - qJD(5) * t24) - t10 * mrSges(6,3) + t11 * mrSges(6,2) + t213 * mrSges(5,2) - t216 * mrSges(5,1) - pkin(4) * t56 + t241 + t242 - t22 * t169 + t107 * t285 + t89 * t284 - t215 * t283 + t4 * t177 + t76 * t282 + t39 * t281 + t40 * t280; qJ(5) * t60 + qJD(5) * t116 - t105 * t169 + t167 * t285 + t83 * t177 + t160 * t282 + t271 * t103 + (-mrSges(6,3) + mrSges(5,2)) * t102 + m(7) * (qJ(5) * t83 + qJD(5) * t105) + m(6) * (-pkin(4) * t103 - qJ(5) * t102 + qJD(5) * t130) + (pkin(4) * t160 - qJ(5) * t161 - qJD(5) * t166) * mrSges(6,1) + (t166 * t284 + t161 * t281 - t17 * mrSges(7,3) + t287 + (t166 * t280 - t44 * mrSges(7,3) - t94 / 0.2e1) * qJD(6) - (m(7) * (t17 + t259) + qJD(6) * t118 + t86) * t286) * t204 + (t166 * t283 + t161 * t280 - t16 * mrSges(7,3) - t48 / 0.2e1 + (-t166 * t179 / 0.2e1 + t43 * mrSges(7,3) - t95 / 0.2e1) * qJD(6) - (m(7) * (t16 - t260) + t85 - qJD(6) * t117) * t286) * t201 + t296; 0; -0.2e1 * qJ(5) * t169 + t171 * t201 - t172 * t204 + (-t179 * t204 - t180 * t201) * qJD(6) + 0.2e1 * (mrSges(6,3) + t177 + (m(6) + m(7)) * qJ(5)) * qJD(5); t201 * t19 + t204 * t20 + t222 * qJD(6) + m(7) * (-qJD(6) * t227 + t295) + m(6) * t11 + t56; -t160 * mrSges(6,1) + t201 * t85 + t204 * t86 + (-t117 * t201 + t118 * t204) * qJD(6) + m(7) * (t16 * t201 + t17 * t204 + (-t201 * t43 + t204 * t44) * qJD(6)) + m(6) * t103; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t17 - mrSges(7,2) * t16 + t47; t169; ((mrSges(7,2) * t286 - Ifges(7,6)) * t204 + (mrSges(7,1) * t286 - Ifges(7,5)) * t201) * qJD(6); -t177 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
