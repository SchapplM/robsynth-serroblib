% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP7
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:56
% EndTime: 2019-12-31 20:00:11
% DurationCPUTime: 6.71s
% Computational Cost: add. (3595->428), mult. (9492->563), div. (0->0), fcn. (6352->6), ass. (0->207)
t287 = Ifges(5,1) + Ifges(6,1);
t286 = Ifges(6,4) + Ifges(5,5);
t288 = Ifges(5,6) - Ifges(6,6);
t135 = sin(qJ(4));
t137 = cos(qJ(4));
t270 = -t288 * t135 + t286 * t137;
t220 = Ifges(6,5) * t135;
t222 = Ifges(5,4) * t135;
t269 = t287 * t137 + t220 - t222;
t136 = sin(qJ(2));
t138 = cos(qJ(2));
t214 = sin(pkin(8));
t215 = cos(pkin(8));
t141 = -t136 * t214 + t138 * t215;
t106 = t141 * qJD(1);
t282 = qJD(4) - t106;
t247 = -t106 / 0.2e1;
t285 = Ifges(4,2) * t247 - Ifges(4,6) * qJD(2) / 0.2e1;
t284 = t135 * t286 + t137 * t288;
t219 = Ifges(6,5) * t137;
t221 = Ifges(5,4) * t137;
t283 = t135 * t287 - t219 + t221;
t132 = -pkin(2) * t138 - pkin(1);
t205 = qJD(1) * t132;
t122 = qJD(3) + t205;
t212 = Ifges(4,5) * qJD(2);
t281 = -t212 / 0.2e1 - t122 * mrSges(4,2);
t118 = t215 * t136 + t214 * t138;
t107 = t118 * qJD(1);
t150 = t137 * qJD(2) - t107 * t135;
t109 = t141 * qJD(2);
t99 = qJD(1) * t109;
t55 = qJD(4) * t150 + t137 * t99;
t89 = qJD(2) * t135 + t107 * t137;
t56 = qJD(4) * t89 + t135 * t99;
t108 = t118 * qJD(2);
t98 = qJD(1) * t108;
t278 = t286 * t98 + (-Ifges(5,4) + Ifges(6,5)) * t56 + t287 * t55;
t238 = Ifges(6,5) * t150;
t87 = Ifges(5,4) * t150;
t277 = t286 * t282 + t287 * t89 - t238 + t87;
t35 = mrSges(5,1) * t98 - mrSges(5,3) * t55;
t36 = -t98 * mrSges(6,1) + t55 * mrSges(6,2);
t276 = -t35 + t36;
t37 = -mrSges(5,2) * t98 - mrSges(5,3) * t56;
t38 = -mrSges(6,2) * t56 + mrSges(6,3) * t98;
t275 = t37 + t38;
t216 = t107 * mrSges(4,3);
t229 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t150 + mrSges(5,2) * t89 + t216;
t58 = mrSges(6,2) * t150 + mrSges(6,3) * t282;
t241 = mrSges(5,3) * t150;
t59 = -mrSges(5,2) * t282 + t241;
t228 = t58 + t59;
t240 = mrSges(5,3) * t89;
t60 = mrSges(5,1) * t282 - t240;
t61 = -mrSges(6,1) * t282 + t89 * mrSges(6,2);
t227 = t60 - t61;
t158 = pkin(4) * t135 - qJ(5) * t137;
t230 = -qJ(3) - pkin(6);
t125 = t230 * t136;
t120 = qJD(1) * t125;
t126 = t230 * t138;
t121 = qJD(1) * t126;
t181 = t215 * t121;
t81 = t120 * t214 - t181;
t273 = -qJD(5) * t135 + t158 * t282 - t81;
t62 = -t106 * pkin(3) - t107 * pkin(7) + t122;
t115 = qJD(2) * pkin(2) + t120;
t80 = t214 * t115 - t181;
t75 = qJD(2) * pkin(7) + t80;
t21 = -t135 * t75 + t137 * t62;
t272 = qJD(5) - t21;
t172 = mrSges(6,1) * t135 - mrSges(6,3) * t137;
t174 = mrSges(5,1) * t135 + mrSges(5,2) * t137;
t110 = t214 * t121;
t79 = t115 * t215 + t110;
t74 = -qJD(2) * pkin(3) - t79;
t23 = -pkin(4) * t150 - t89 * qJ(5) + t74;
t271 = -t23 * t172 - t74 * t174;
t200 = qJD(4) * t137;
t208 = t106 * t137;
t268 = t200 - t208;
t201 = qJD(4) * t135;
t209 = t106 * t135;
t267 = -t201 + t209;
t182 = qJD(2) * t230;
t105 = -t136 * qJD(3) + t138 * t182;
t139 = qJD(1) * t105;
t104 = qJD(3) * t138 + t136 * t182;
t93 = t104 * qJD(1);
t54 = t139 * t214 + t215 * t93;
t199 = qJD(1) * qJD(2);
t186 = t136 * t199;
t178 = pkin(2) * t186;
t57 = pkin(3) * t98 - pkin(7) * t99 + t178;
t3 = t135 * t57 + t137 * t54 + t62 * t200 - t201 * t75;
t22 = t135 * t62 + t137 * t75;
t4 = -qJD(4) * t22 - t135 * t54 + t137 * t57;
t266 = -t4 * t135 + t3 * t137;
t1 = qJ(5) * t98 + qJD(5) * t282 + t3;
t2 = -pkin(4) * t98 - t4;
t265 = t1 * t137 + t2 * t135;
t203 = qJD(1) * t138;
t211 = Ifges(3,6) * qJD(2);
t224 = Ifges(3,4) * t136;
t263 = t211 / 0.2e1 + (t138 * Ifges(3,2) + t224) * qJD(1) / 0.2e1 + pkin(6) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t203);
t78 = -pkin(3) * t141 - t118 * pkin(7) + t132;
t85 = t125 * t214 - t126 * t215;
t226 = t135 * t78 + t137 * t85;
t67 = t104 * t215 + t105 * t214;
t202 = qJD(2) * t136;
t69 = pkin(2) * t202 + pkin(3) * t108 - pkin(7) * t109;
t9 = -qJD(4) * t226 - t135 * t67 + t137 * t69;
t223 = Ifges(4,4) * t107;
t255 = -t223 / 0.2e1 + t285;
t161 = Ifges(6,3) * t135 + t219;
t167 = -Ifges(5,2) * t135 + t221;
t244 = t135 / 0.2e1;
t245 = -t135 / 0.2e1;
t251 = t89 / 0.2e1;
t253 = -t150 / 0.2e1;
t254 = t150 / 0.2e1;
t86 = Ifges(6,5) * t89;
t29 = Ifges(6,6) * t282 - Ifges(6,3) * t150 + t86;
t239 = Ifges(5,4) * t89;
t32 = Ifges(5,2) * t150 + Ifges(5,6) * t282 + t239;
t262 = t161 * t253 + t167 * t254 + t29 * t244 + t32 * t245 - t271 + t269 * t251 + t270 * t282 / 0.2e1;
t15 = -pkin(4) * t282 + t272;
t16 = qJ(5) * t282 + t22;
t261 = t122 * mrSges(4,1) + t21 * mrSges(5,1) - t15 * mrSges(6,1) - t22 * mrSges(5,2) + t16 * mrSges(6,3) + t255;
t260 = -0.2e1 * pkin(1);
t258 = t55 / 0.2e1;
t257 = -t56 / 0.2e1;
t256 = t56 / 0.2e1;
t252 = -t89 / 0.2e1;
t250 = t98 / 0.2e1;
t249 = -t282 / 0.2e1;
t246 = -t107 / 0.2e1;
t243 = -t137 / 0.2e1;
t242 = t137 / 0.2e1;
t204 = qJD(1) * t136;
t237 = pkin(6) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t204);
t53 = -t215 * t139 + t214 * t93;
t84 = -t215 * t125 - t126 * t214;
t231 = t53 * t84;
t198 = pkin(2) * t204;
t68 = pkin(3) * t107 - pkin(7) * t106 + t198;
t82 = t120 * t215 + t110;
t27 = t135 * t68 + t137 * t82;
t225 = Ifges(4,1) * t107;
t101 = Ifges(4,4) * t106;
t217 = t106 * mrSges(4,3);
t213 = Ifges(3,5) * qJD(2);
t197 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t196 = -Ifges(5,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t195 = Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t192 = t215 * pkin(2);
t191 = t214 * pkin(2);
t187 = t98 * mrSges(4,1) + t99 * mrSges(4,2);
t185 = t138 * t199;
t131 = -t192 - pkin(3);
t177 = -t1 * t135 + t137 * t2;
t176 = -t3 * t135 - t4 * t137;
t66 = t104 * t214 - t215 * t105;
t175 = mrSges(5,1) * t137 - mrSges(5,2) * t135;
t173 = mrSges(6,1) * t137 + mrSges(6,3) * t135;
t166 = Ifges(5,2) * t137 + t222;
t160 = -Ifges(6,3) * t137 + t220;
t159 = t137 * pkin(4) + t135 * qJ(5);
t157 = -t16 * t135 + t15 * t137;
t156 = t135 * t15 + t137 * t16;
t155 = -t22 * t135 - t21 * t137;
t154 = t135 * t21 - t137 * t22;
t26 = -t135 * t82 + t137 * t68;
t39 = -t135 * t85 + t137 * t78;
t8 = t135 * t69 + t137 * t67 + t78 * t200 - t201 * t85;
t140 = t4 * mrSges(5,1) - t2 * mrSges(6,1) - t3 * mrSges(5,2) + t1 * mrSges(6,3);
t133 = Ifges(3,4) * t203;
t116 = -t159 + t131;
t114 = Ifges(3,1) * t204 + t133 + t213;
t97 = Ifges(6,2) * t98;
t96 = Ifges(5,3) * t98;
t91 = -qJD(2) * mrSges(4,2) + t217;
t76 = -mrSges(4,1) * t106 + mrSges(4,2) * t107;
t72 = t101 + t212 + t225;
t51 = Ifges(6,4) * t55;
t50 = Ifges(5,5) * t55;
t49 = Ifges(5,6) * t56;
t48 = Ifges(6,6) * t56;
t43 = -mrSges(6,1) * t150 - mrSges(6,3) * t89;
t42 = pkin(4) * t89 - qJ(5) * t150;
t41 = t118 * t158 + t84;
t31 = Ifges(6,4) * t89 + Ifges(6,2) * t282 - Ifges(6,6) * t150;
t30 = Ifges(5,5) * t89 + Ifges(5,6) * t150 + Ifges(5,3) * t282;
t25 = pkin(4) * t141 - t39;
t24 = -qJ(5) * t141 + t226;
t20 = -pkin(4) * t107 - t26;
t19 = qJ(5) * t107 + t27;
t18 = mrSges(5,1) * t56 + mrSges(5,2) * t55;
t17 = mrSges(6,1) * t56 - mrSges(6,3) * t55;
t12 = t55 * Ifges(5,4) - t56 * Ifges(5,2) + t98 * Ifges(5,6);
t11 = t55 * Ifges(6,5) + t98 * Ifges(6,6) + t56 * Ifges(6,3);
t10 = t158 * t109 + (qJD(4) * t159 - qJD(5) * t137) * t118 + t66;
t7 = t56 * pkin(4) - t55 * qJ(5) - t89 * qJD(5) + t53;
t6 = -pkin(4) * t108 - t9;
t5 = qJ(5) * t108 - qJD(5) * t141 + t8;
t13 = [t84 * t18 + t5 * t58 + t8 * t59 + t9 * t60 + t6 * t61 + t10 * t43 + t25 * t36 + t24 * t38 + t39 * t35 + t226 * t37 + t41 * t17 + t132 * t187 + t67 * t91 + t229 * t66 + (t84 * t99 - t85 * t98) * mrSges(4,3) + m(4) * (t54 * t85 - t66 * t79 + t67 * t80 + t231) + m(5) * (t21 * t9 + t22 * t8 + t226 * t3 + t39 * t4 + t66 * t74 + t231) + m(6) * (t1 * t24 + t10 * t23 + t15 * t6 + t16 * t5 + t2 * t25 + t41 * t7) + (t114 / 0.2e1 - t237 + t213 / 0.2e1 + (mrSges(3,2) * t260 + 0.3e1 / 0.2e1 * Ifges(3,4) * t138) * qJD(1)) * t138 * qJD(2) + (-t80 * mrSges(4,3) + t31 / 0.2e1 + t30 / 0.2e1 + t197 * t89 - t196 * t150 + t195 * t282 + t261 + t255) * t108 + (t225 / 0.2e1 - t79 * mrSges(4,3) + t262 + t277 * t242 + t155 * mrSges(5,3) + t157 * mrSges(6,2) + t72 / 0.2e1 + t101 / 0.2e1 - t281) * t109 - (-t54 * mrSges(4,3) - Ifges(4,4) * t99 + t50 / 0.2e1 - t49 / 0.2e1 + t96 / 0.2e1 + t51 / 0.2e1 + t97 / 0.2e1 + t48 / 0.2e1 + t196 * t56 + t197 * t55 + (Ifges(4,2) + t195) * t98 + t140) * t141 + (t11 * t244 + t7 * t172 + t161 * t256 + t167 * t257 - Ifges(4,4) * t98 + Ifges(4,1) * t99 + (mrSges(4,3) + t174) * t53 + t176 * mrSges(5,3) + t177 * mrSges(6,2) + (-mrSges(6,2) * t156 + mrSges(5,3) * t154 + t160 * t254 + t166 * t253 + t173 * t23 + t175 * t74 + t243 * t32 + t249 * t284 + t252 * t283) * qJD(4) + t269 * t258 + t270 * t250 + (qJD(4) * t277 + t12) * t245 + (qJD(4) * t29 + t278) * t242) * t118 + (-t211 / 0.2e1 + (mrSges(3,1) * t260 - 0.3e1 / 0.2e1 * t224 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t138) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t141 + mrSges(4,2) * t118) + t76 + m(4) * (t122 + t205)) * pkin(2) - t263) * t202; (t72 + t101) * t247 - t229 * t81 - t54 * mrSges(4,2) - t19 * t58 - t27 * t59 - t26 * t60 - t20 * t61 - (-Ifges(5,6) * t253 - Ifges(6,6) * t254 - t286 * t252 + (-Ifges(5,3) - Ifges(6,2)) * t249 + t261 + t285) * t107 + (-mrSges(3,1) * t185 + mrSges(3,2) * t186) * pkin(6) + t273 * t43 + (t15 * t268 + t16 * t267 + t265) * mrSges(6,2) + (-t21 * t268 + t22 * t267 + t266) * mrSges(5,3) + (t32 / 0.2e1 - t29 / 0.2e1) * t209 + t263 * t204 + t283 * t258 + t284 * t250 + (-t228 * t201 + t275 * t137 + t276 * t135 - t227 * t200 + (m(5) * t155 + m(6) * t157) * qJD(4) + t265 * m(6) + t266 * m(5)) * (t191 + pkin(7)) + (-mrSges(4,1) - t175) * t53 + (Ifges(4,1) * t246 + t161 * t254 + t167 * t253 + t270 * t249 + t269 * t252 + t271 + t281) * t106 + (-t122 * t198 + t79 * t81 - t80 * t82 + (t214 * t54 - t215 * t53) * pkin(2)) * m(4) - (-Ifges(3,2) * t204 + t114 + t133) * t203 / 0.2e1 + (-t136 * (Ifges(3,1) * t138 - t224) / 0.2e1 + pkin(1) * (mrSges(3,1) * t136 + mrSges(3,2) * t138)) * qJD(1) ^ 2 - t7 * t173 + t278 * t244 - t76 * t198 - (Ifges(3,5) * t138 - Ifges(3,6) * t136) * t199 / 0.2e1 + t262 * qJD(4) + (t116 * t7 - t15 * t20 - t16 * t19 + t23 * t273) * m(6) + (t131 * t53 - t21 * t26 - t22 * t27 - t74 * t81) * m(5) - Ifges(3,6) * t186 + (-t223 + t31 + t30) * t246 + (-t191 * t98 - t192 * t99) * mrSges(4,3) + t80 * t216 + t160 * t256 + t166 * t257 + t12 * t242 + t11 * t243 + t203 * t237 + t79 * t217 + (t200 / 0.2e1 - t208 / 0.2e1) * t277 + t131 * t18 + t116 * t17 - Ifges(4,6) * t98 + Ifges(4,5) * t99 - t82 * t91 + Ifges(3,5) * t185; -t106 * t91 - (t43 + t229) * t107 + (t228 * t282 - t276) * t137 + (-t227 * t282 + t275) * t135 + t187 + (-t107 * t23 + t156 * t282 - t177) * m(6) + (-t74 * t107 - t154 * t282 - t176) * m(5) + (-t106 * t80 + t107 * t79 + t178) * m(4); t51 + t50 - t49 + t48 + qJD(5) * t58 - t42 * t43 - pkin(4) * t36 + qJ(5) * t38 + t140 + (-t15 * t150 + t16 * t89) * mrSges(6,2) + t32 * t251 + (Ifges(6,3) * t89 + t238) * t254 + (t227 + t240) * t22 + (-t228 + t241) * t21 + t96 + t97 - t23 * (mrSges(6,1) * t89 - mrSges(6,3) * t150) - t74 * (mrSges(5,1) * t89 + mrSges(5,2) * t150) + (t150 * t286 - t288 * t89) * t249 + (-pkin(4) * t2 + qJ(5) * t1 - t15 * t22 + t16 * t272 - t23 * t42) * m(6) + (-Ifges(5,2) * t89 + t277 + t87) * t253 + (t150 * t287 - t239 + t29 + t86) * t252; -t282 * t58 + t89 * t43 + 0.2e1 * (t2 / 0.2e1 + t16 * t249 + t23 * t251) * m(6) + t36;];
tauc = t13(:);
