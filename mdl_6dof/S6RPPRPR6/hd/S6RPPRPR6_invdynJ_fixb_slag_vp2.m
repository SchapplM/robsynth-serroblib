% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:43
% EndTime: 2019-03-09 01:50:55
% DurationCPUTime: 8.24s
% Computational Cost: add. (2513->482), mult. (4429->619), div. (0->0), fcn. (2049->6), ass. (0->215)
t303 = mrSges(6,2) - mrSges(5,1);
t297 = Ifges(6,4) - Ifges(5,5);
t296 = Ifges(6,5) - Ifges(5,6);
t119 = sin(qJ(6));
t122 = cos(qJ(6));
t120 = sin(qJ(4));
t205 = t120 * qJD(1);
t213 = qJD(4) * t119;
t67 = t122 * t205 - t213;
t123 = cos(qJ(4));
t202 = qJD(1) * qJD(4);
t75 = qJDD(1) * t120 + t123 * t202;
t24 = qJD(6) * t67 + qJDD(4) * t122 + t119 * t75;
t74 = -t123 * qJDD(1) + t120 * t202;
t64 = qJDD(6) - t74;
t11 = mrSges(7,1) * t64 - mrSges(7,3) * t24;
t211 = qJD(4) * t122;
t68 = t119 * t205 + t211;
t25 = -qJD(6) * t68 - qJDD(4) * t119 + t122 * t75;
t12 = -mrSges(7,2) * t64 + mrSges(7,3) * t25;
t204 = t123 * qJD(1);
t92 = qJD(6) + t204;
t35 = -mrSges(7,2) * t92 + mrSges(7,3) * t67;
t36 = mrSges(7,1) * t92 - mrSges(7,3) * t68;
t280 = -t119 * t36 + t122 * t35;
t302 = -qJD(6) * t280 - t122 * t11 - t119 * t12;
t118 = pkin(1) + qJ(3);
t106 = t123 * qJ(5);
t262 = pkin(8) * t120;
t151 = -t106 + t262;
t219 = pkin(4) * t205 - qJD(2);
t32 = (t151 + t118) * qJD(1) + t219;
t265 = pkin(4) + pkin(8);
t102 = qJD(1) * qJ(2) + qJD(3);
t86 = -qJD(1) * pkin(7) + t102;
t43 = (pkin(5) * qJD(1) - t86) * t123;
t300 = t43 + qJD(5);
t33 = -qJD(4) * t265 + t300;
t10 = t119 * t33 + t122 * t32;
t212 = qJD(4) * t120;
t112 = qJD(1) * qJD(2);
t287 = qJDD(1) * qJ(2) + t112;
t85 = qJDD(3) + t287;
t69 = -qJDD(1) * pkin(7) + t85;
t30 = t123 * t69 - t86 * t212;
t170 = qJDD(5) - t30;
t13 = -pkin(5) * t74 - qJDD(4) * t265 + t170;
t143 = qJDD(1) * t118 - qJDD(2);
t175 = -qJD(5) * t123 + qJD(3);
t127 = qJ(5) * t74 + qJD(1) * t175 + t143;
t8 = t265 * t75 + t127;
t9 = -t119 * t32 + t122 * t33;
t1 = qJD(6) * t9 + t119 * t13 + t122 * t8;
t2 = -qJD(6) * t10 - t119 * t8 + t122 * t13;
t171 = t1 * t119 + t122 * t2;
t207 = qJD(6) * t122;
t209 = qJD(6) * t119;
t301 = -t10 * t207 + t9 * t209 - t171;
t121 = sin(qJ(1));
t124 = cos(qJ(1));
t172 = g(1) * t124 + g(2) * t121;
t266 = -m(6) - m(7);
t267 = -m(5) - m(4);
t196 = t266 + t267;
t272 = t24 / 0.2e1;
t271 = t25 / 0.2e1;
t270 = t64 / 0.2e1;
t299 = mrSges(3,2) - mrSges(4,3);
t298 = -mrSges(4,2) - mrSges(3,3);
t47 = mrSges(6,1) * t75 - qJDD(4) * mrSges(6,3);
t295 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t75 - t47;
t251 = Ifges(5,4) * t120;
t162 = t123 * Ifges(5,1) - t251;
t294 = Ifges(5,5) * qJD(4) + t68 * Ifges(7,5) + t67 * Ifges(7,6) + t92 * Ifges(7,3) + qJD(1) * t162;
t194 = mrSges(6,1) * t204;
t293 = -mrSges(5,3) * t204 - t303 * qJD(4) - t194;
t29 = -mrSges(7,1) * t67 + mrSges(7,2) * t68;
t195 = mrSges(6,1) * t205;
t82 = -qJD(4) * mrSges(6,3) + t195;
t292 = -t82 + t29;
t80 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t205;
t291 = -t82 + t80;
t210 = qJD(4) * t123;
t290 = pkin(4) * t210 + qJ(5) * t212;
t26 = -qJDD(4) * pkin(4) + t170;
t214 = qJD(4) * qJ(5);
t76 = t120 * t86;
t51 = -t76 - t214;
t289 = -qJD(4) * t51 - t26;
t288 = t123 * t172;
t246 = Ifges(6,6) * t123;
t285 = t120 * (-Ifges(5,2) * t123 - t251) + t123 * (Ifges(6,2) * t120 + t246);
t284 = t297 * t120 + t296 * t123;
t31 = t120 * t69 + t86 * t210;
t282 = t120 * t31 + t123 * t30;
t23 = -qJDD(4) * qJ(5) - qJD(4) * qJD(5) - t31;
t281 = -t120 * t23 - t123 * t26;
t152 = t120 * Ifges(6,3) - t246;
t263 = Ifges(7,4) * t68;
t19 = Ifges(7,2) * t67 + Ifges(7,6) * t92 + t263;
t63 = Ifges(7,4) * t67;
t20 = Ifges(7,1) * t68 + Ifges(7,5) * t92 + t63;
t278 = Ifges(6,5) * qJD(4) + qJD(1) * t152 + t119 * t20 + t122 * t19;
t167 = t120 * mrSges(5,1) + t123 * mrSges(5,2);
t277 = -mrSges(2,1) - t167 + t299;
t220 = t123 * t124;
t222 = t121 * t123;
t223 = t120 * t124;
t276 = -g(1) * (pkin(4) * t220 + qJ(5) * t223) - g(2) * (t120 * t121 * qJ(5) + pkin(4) * t222);
t48 = -t74 * mrSges(6,1) + qJDD(4) * mrSges(6,2);
t275 = t302 - t48;
t274 = m(7) * pkin(5) + mrSges(6,1) + mrSges(2,2) + mrSges(5,3) + t298;
t126 = qJD(1) ^ 2;
t273 = Ifges(7,1) * t272 + Ifges(7,4) * t271 + Ifges(7,5) * t270;
t268 = t68 / 0.2e1;
t109 = t120 * pkin(4);
t117 = qJ(2) - pkin(7);
t255 = pkin(5) - t117;
t250 = Ifges(5,4) * t123;
t249 = Ifges(7,4) * t119;
t248 = Ifges(7,4) * t122;
t247 = Ifges(6,6) * t120;
t242 = t120 * mrSges(6,2);
t241 = t120 * mrSges(7,3);
t233 = t123 * t86;
t71 = pkin(4) * t204 + qJ(5) * t205;
t42 = -pkin(5) * t205 + t76;
t37 = t42 + t214;
t231 = qJD(4) * t37;
t229 = qJDD(1) * pkin(1);
t226 = t119 * t120;
t225 = t119 * t123;
t224 = t120 * t122;
t221 = t122 * t123;
t218 = t109 - t106;
t217 = t124 * pkin(1) + t121 * qJ(2);
t216 = qJD(2) * t120;
t215 = qJD(2) * t123;
t208 = qJD(6) * t120;
t203 = qJD(1) * qJD(3);
t199 = Ifges(7,5) * t24 + Ifges(7,6) * t25 + Ifges(7,3) * t64;
t198 = t80 + t292;
t192 = t124 * qJ(3) + t217;
t186 = -t209 / 0.2e1;
t108 = t124 * qJ(2);
t181 = -pkin(7) * t124 + t108;
t180 = qJD(4) * t255;
t179 = (t120 ^ 2 + t123 ^ 2) * t86;
t178 = -t202 / 0.2e1;
t173 = -t218 - t262;
t169 = t10 * t122 - t9 * t119;
t168 = mrSges(5,1) * t123 - mrSges(5,2) * t120;
t166 = mrSges(7,1) * t122 - mrSges(7,2) * t119;
t165 = mrSges(7,1) * t119 + mrSges(7,2) * t122;
t164 = -t123 * mrSges(6,2) + mrSges(6,3) * t120;
t163 = t123 * mrSges(6,3) + t242;
t161 = Ifges(7,1) * t122 - t249;
t160 = Ifges(7,1) * t119 + t248;
t159 = -t120 * Ifges(5,2) + t250;
t157 = -Ifges(7,2) * t119 + t248;
t156 = Ifges(7,2) * t122 + t249;
t154 = Ifges(7,5) * t122 - Ifges(7,6) * t119;
t153 = Ifges(7,5) * t119 + Ifges(7,6) * t122;
t148 = t119 * t35 + t122 * t36;
t50 = -t173 + t118;
t79 = t255 * t123;
t28 = t119 * t79 + t122 * t50;
t27 = -t119 * t50 + t122 * t79;
t49 = -qJD(4) * pkin(4) + qJD(5) - t233;
t147 = -t120 * t51 - t123 * t49;
t146 = t49 * t120 - t51 * t123;
t70 = t143 + t203;
t87 = qJD(1) * t118 - qJD(2);
t145 = qJD(3) * t87 + t118 * t70;
t144 = -pkin(7) * t121 + t192;
t72 = t163 * qJD(1);
t142 = -t280 + t72;
t38 = (-t106 + t118) * qJD(1) + t219;
t140 = t38 * t164;
t139 = t87 * t168;
t137 = t123 * (-Ifges(5,1) * t120 - t250);
t135 = t148 - t293;
t134 = -t119 * t208 + t122 * t210;
t133 = t119 * t210 + t120 * t207;
t132 = -Ifges(7,5) * t120 + t123 * t160;
t131 = -Ifges(7,6) * t120 + t123 * t156;
t130 = -Ifges(7,3) * t120 + t123 * t153;
t129 = qJD(6) * t169 + t171;
t103 = qJDD(2) - t229;
t98 = Ifges(6,6) * t205;
t78 = t255 * t120;
t77 = t218 + t118;
t73 = t167 * qJD(1);
t62 = t166 * t120;
t60 = t119 * t222 - t122 * t124;
t59 = t119 * t124 + t121 * t221;
t58 = t119 * t220 + t121 * t122;
t57 = t119 * t121 - t122 * t220;
t56 = Ifges(6,4) * qJD(4) - Ifges(6,2) * t204 + t98;
t53 = Ifges(5,6) * qJD(4) + qJD(1) * t159;
t45 = qJDD(4) * mrSges(5,1) + mrSges(5,3) * t74;
t44 = pkin(8) * t204 + t71;
t41 = t175 + t290;
t40 = -t123 * t180 + t216;
t39 = -t120 * t180 - t215;
t34 = qJD(3) + (qJD(4) * pkin(8) - qJD(5)) * t123 + t290;
t17 = t119 * t42 + t122 * t44;
t16 = -t119 * t44 + t122 * t42;
t15 = pkin(4) * t75 + t127;
t14 = -pkin(5) * t75 - t23;
t7 = -mrSges(7,1) * t25 + mrSges(7,2) * t24;
t6 = -qJD(6) * t28 - t119 * t34 + t122 * t39;
t5 = qJD(6) * t27 + t119 * t39 + t122 * t34;
t3 = t24 * Ifges(7,4) + t25 * Ifges(7,2) + t64 * Ifges(7,6);
t4 = [(qJD(4) * t132 + t161 * t208) * t268 + 0.2e1 * t287 * mrSges(3,3) + (-m(5) * t181 - t60 * mrSges(7,1) - t59 * mrSges(7,2) + t266 * (t121 * t106 + t181) + (-m(3) - m(4)) * t108 + t274 * t124 + (m(3) * pkin(1) + m(6) * t109 - t163 - (-m(7) * t265 - mrSges(7,3)) * t120 - t196 * t118 - t277) * t121) * g(1) + t70 * t167 - t75 * t159 / 0.2e1 - t74 * t162 / 0.2e1 - t15 * t163 + t75 * t152 / 0.2e1 + m(4) * (qJ(2) * t85 + qJD(2) * t102 + t145) + (t103 - t229) * mrSges(3,2) + (m(6) * (qJD(4) * t146 + t281) + m(5) * t282 + (-t48 + t45) * t123 + t291 * t210 - t293 * t212 + t295 * t120) * t117 + t37 * (-mrSges(7,1) * t134 + mrSges(7,2) * t133) + qJD(4) * t139 + qJD(4) * t140 + t284 * qJD(4) ^ 2 / 0.2e1 + t285 * t178 + (t85 + t287) * mrSges(4,2) - t282 * mrSges(5,3) + (-t49 * mrSges(6,1) - t294 / 0.2e1 + t56 / 0.2e1 - t9 * mrSges(7,1) + t10 * mrSges(7,2)) * t212 + t291 * t216 + t293 * t215 + (t203 + t70) * mrSges(4,3) + (mrSges(4,3) * t118 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) + t1 * (-mrSges(7,2) * t123 + mrSges(7,3) * t224) + m(5) * (qJD(2) * t179 + t145) + m(6) * (t147 * qJD(2) + t15 * t77 + t38 * t41) + t67 * (qJD(4) * t131 + t157 * t208) / 0.2e1 + t92 * (qJD(4) * t130 + t154 * t208) / 0.2e1 + (t51 * mrSges(6,1) + t278 / 0.2e1 - t53 / 0.2e1) * t210 + (Ifges(7,3) * t123 + t120 * t153) * t270 + (Ifges(7,6) * t123 + t120 * t156) * t271 + (Ifges(7,5) * t123 + t120 * t160) * t272 + t226 * t273 + (-Ifges(5,1) * t74 - Ifges(5,4) * t75 + Ifges(5,5) * qJDD(4) + t199) * t123 / 0.2e1 + m(7) * (t1 * t28 + t10 * t5 - t14 * t78 + t2 * t27 + t37 * t40 + t6 * t9) + (t120 * (Ifges(6,3) * t123 + t247) + t137) * t202 / 0.2e1 + (qJD(6) * t20 + t3) * t224 / 0.2e1 - t74 * t123 * Ifges(6,2) + m(3) * (-pkin(1) * t103 + (t287 + t112) * qJ(2)) - t123 * (Ifges(6,4) * qJDD(4) + Ifges(6,6) * t75) / 0.2e1 + t74 * t247 / 0.2e1 - t281 * mrSges(6,1) + t120 * (Ifges(6,5) * qJDD(4) + Ifges(6,6) * t74 + Ifges(6,3) * t75) / 0.2e1 - t120 * (-Ifges(5,4) * t74 - Ifges(5,2) * t75 + Ifges(5,6) * qJDD(4)) / 0.2e1 + t118 * (t75 * mrSges(5,1) - t74 * mrSges(5,2)) + t77 * (-t75 * mrSges(6,2) + t74 * mrSges(6,3)) - t78 * t7 - t41 * t72 + qJD(3) * t73 - t14 * t62 + t40 * t29 + t5 * t35 + t6 * t36 + t27 * t11 + t28 * t12 + t120 * t19 * t186 + t10 * mrSges(7,3) * t134 + (t296 * t120 - t297 * t123) * qJDD(4) / 0.2e1 + (-m(3) * t217 - m(4) * t192 - m(5) * t144 + t58 * mrSges(7,1) - t57 * mrSges(7,2) + t266 * (pkin(4) * t223 + t144) + t274 * t121 + (-m(7) * t151 - t241 + t242 - (-m(6) * qJ(5) - mrSges(6,3)) * t123 + t277) * t124) * g(2) - t9 * mrSges(7,3) * t133 + t2 * (mrSges(7,1) * t123 - mrSges(7,3) * t226); t119 * t11 - t122 * t12 + t303 * t75 + (-mrSges(6,3) + mrSges(5,2)) * t74 + t299 * qJDD(1) + t148 * qJD(6) + (-m(3) * qJ(2) + t298) * t126 + m(7) * (-t1 * t122 + t2 * t119 + (t10 * t119 + t9 * t122) * qJD(6)) + m(3) * t103 - m(6) * t15 - m(5) * t70 + (-t198 * t120 + t135 * t123 - m(7) * (-t10 * t225 + t120 * t37 - t221 * t9) - m(6) * t147 - m(5) * t179) * qJD(1) + (-g(1) * t121 + g(2) * t124) * (m(3) - t196) + (-qJD(1) * t102 - t70) * m(4); m(4) * t85 - t126 * mrSges(4,3) + qJDD(1) * mrSges(4,2) + (t7 + t135 * qJD(4) + m(7) * (t10 * t213 + t211 * t9 + t14) + m(6) * (qJD(4) * t49 - t23) + m(5) * t31 + t295) * t120 + (t45 + t198 * qJD(4) + m(7) * (t231 + t301) + m(6) * t289 + m(5) * t30 + t275) * t123 + (-m(6) * t38 - m(7) * t169 + t267 * t87 + t142 - t73) * qJD(1) + t172 * t196; (qJ(5) * t14 - t10 * t17 - t16 * t9 + t300 * t37 + t276) * m(7) + t301 * mrSges(7,3) + (m(6) * t218 - m(7) * t173 - t123 * t165 - t163 + t167 + t241) * g(3) + t14 * t165 + t53 * t204 / 0.2e1 + (-(m(7) * pkin(8) + mrSges(7,3)) * t123 - t165 * t120 - t164 - t168) * t172 + (-m(7) * t129 + t302) * t265 - t51 * t194 + t284 * t178 + (t285 / 0.2e1 - t137 / 0.2e1) * t126 + t92 * t37 * t166 - t291 * t233 + t292 * qJD(5) + t293 * t76 + t294 * t205 / 0.2e1 + (-pkin(4) * t26 - qJ(5) * t23 - qJD(5) * t51 - t146 * t86 - t38 * t71 + t276) * m(6) - t278 * t204 / 0.2e1 - t19 * t207 / 0.2e1 + t154 * t270 + t157 * t271 + t161 * t272 + t122 * t273 - (Ifges(6,3) * t204 + t56 + t98) * t205 / 0.2e1 - (t153 * t92 + t156 * t67 + t160 * t68) * qJD(6) / 0.2e1 - (t130 * t92 + t131 * t67 + t132 * t68) * qJD(1) / 0.2e1 - t119 * t3 / 0.2e1 + t71 * t72 - pkin(4) * t48 + t43 * t29 + t30 * mrSges(5,1) - t31 * mrSges(5,2) - t17 * t35 - t16 * t36 - t23 * mrSges(6,3) + t26 * mrSges(6,2) + t49 * t195 + (-t139 - t140 - t10 * (mrSges(7,2) * t120 + mrSges(7,3) * t221) - t9 * (-mrSges(7,1) * t120 - mrSges(7,3) * t225)) * qJD(1) + t296 * t75 + t297 * t74 + (Ifges(6,1) + Ifges(5,3)) * qJDD(4) + (-t47 + t7) * qJ(5) + t20 * t186; t266 * t120 * g(3) - t142 * t204 - t292 * qJD(4) + (t169 * t204 + t129 - t231 + t288) * m(7) + (t204 * t38 + t288 - t289) * m(6) - t275; -t1 * mrSges(7,2) + t2 * mrSges(7,1) - t37 * (mrSges(7,1) * t68 + mrSges(7,2) * t67) - t68 * (Ifges(7,1) * t67 - t263) / 0.2e1 + t19 * t268 - t92 * (Ifges(7,5) * t67 - Ifges(7,6) * t68) / 0.2e1 - t9 * t35 + t10 * t36 - g(1) * (mrSges(7,1) * t57 + mrSges(7,2) * t58) - g(2) * (-mrSges(7,1) * t59 + mrSges(7,2) * t60) - g(3) * t62 + (t10 * t68 + t67 * t9) * mrSges(7,3) + t199 - (-Ifges(7,2) * t68 + t20 + t63) * t67 / 0.2e1;];
tau  = t4;
