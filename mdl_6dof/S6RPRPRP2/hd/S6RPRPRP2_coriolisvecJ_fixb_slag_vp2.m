% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2018-11-23 15:57
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:57:26
% EndTime: 2018-11-23 15:57:33
% DurationCPUTime: 6.63s
% Computational Cost: add. (5052->480), mult. (12401->624), div. (0->0), fcn. (8194->8), ass. (0->229)
t314 = Ifges(6,1) + Ifges(7,1);
t313 = Ifges(7,4) + Ifges(6,5);
t315 = Ifges(6,6) - Ifges(7,6);
t146 = sin(qJ(5));
t148 = cos(qJ(5));
t296 = -t315 * t146 + t313 * t148;
t254 = Ifges(7,5) * t146;
t257 = Ifges(6,4) * t146;
t295 = t314 * t148 + t254 - t257;
t147 = sin(qJ(3));
t149 = cos(qJ(3));
t242 = sin(pkin(10));
t243 = cos(pkin(10));
t152 = -t147 * t242 + t149 * t243;
t118 = t152 * qJD(1);
t307 = qJD(5) - t118;
t278 = -t118 / 0.2e1;
t312 = Ifges(5,2) * t278 - Ifges(5,6) * qJD(3) / 0.2e1;
t311 = t313 * t146 + t315 * t148;
t253 = Ifges(7,5) * t148;
t256 = Ifges(6,4) * t148;
t310 = t314 * t146 - t253 + t256;
t224 = qJD(5) * t148;
t225 = qJD(5) * t146;
t222 = t147 * qJD(4);
t138 = sin(pkin(9)) * pkin(1) + pkin(7);
t130 = t138 * qJD(1);
t192 = qJ(4) * qJD(1) + t130;
t223 = t147 * qJD(2);
t99 = t149 * t192 + t223;
t150 = -qJD(1) * t222 - qJD(3) * t99;
t143 = t149 * qJD(2);
t227 = qJD(3) * t147;
t104 = qJD(3) * t143 - t130 * t227;
t226 = qJD(4) * t149;
t87 = (-qJ(4) * t227 + t226) * qJD(1) + t104;
t33 = t150 * t242 + t243 * t87;
t200 = t243 * t99;
t98 = -t147 * t192 + t143;
t92 = qJD(3) * pkin(3) + t98;
t49 = t242 * t92 + t200;
t46 = qJD(3) * pkin(8) + t49;
t207 = -cos(pkin(9)) * pkin(1) - pkin(2);
t129 = -pkin(3) * t149 + t207;
t230 = qJD(1) * t129;
t117 = qJD(4) + t230;
t128 = t243 * t147 + t242 * t149;
t119 = t128 * qJD(1);
t59 = -t118 * pkin(4) - t119 * pkin(8) + t117;
t120 = t128 * qJD(3);
t109 = qJD(1) * t120;
t121 = t152 * qJD(3);
t110 = qJD(1) * t121;
t221 = qJD(1) * qJD(3);
t199 = t147 * t221;
t191 = pkin(3) * t199;
t68 = pkin(4) * t109 - pkin(8) * t110 + t191;
t3 = t146 * t68 + t148 * t33 + t59 * t224 - t225 * t46;
t20 = t146 * t59 + t148 * t46;
t4 = -qJD(5) * t20 - t146 * t33 + t148 * t68;
t188 = -t146 * t4 + t148 * t3;
t19 = -t146 * t46 + t148 * t59;
t309 = -t19 * t224 - t20 * t225 + t188;
t298 = qJD(6) - t19;
t11 = -pkin(5) * t307 + t298;
t12 = qJ(6) * t307 + t20;
t1 = qJ(6) * t109 + qJD(6) * t307 + t3;
t2 = -pkin(5) * t109 - t4;
t190 = t1 * t148 + t146 * t2;
t308 = t11 * t224 - t12 * t225 + t190;
t240 = Ifges(5,5) * qJD(3);
t306 = -t240 / 0.2e1 - t117 * mrSges(5,2);
t97 = qJD(3) * t146 + t119 * t148;
t67 = qJD(5) * t97 + t110 * t146;
t42 = -mrSges(6,2) * t109 - mrSges(6,3) * t67;
t43 = -mrSges(7,2) * t67 + mrSges(7,3) * t109;
t262 = t42 + t43;
t161 = t148 * qJD(3) - t119 * t146;
t66 = qJD(5) * t161 + t110 * t148;
t40 = mrSges(6,1) * t109 - mrSges(6,3) * t66;
t41 = -t109 * mrSges(7,1) + t66 * mrSges(7,2);
t263 = -t40 + t41;
t305 = t263 * t146 + t262 * t148;
t302 = (-Ifges(6,4) + Ifges(7,5)) * t67 + t314 * t66 + t313 * t109;
t269 = Ifges(7,5) * t161;
t94 = Ifges(6,4) * t161;
t301 = t313 * t307 + t314 * t97 - t269 + t94;
t69 = mrSges(7,2) * t161 + mrSges(7,3) * t307;
t272 = mrSges(6,3) * t161;
t70 = -mrSges(6,2) * t307 + t272;
t261 = t69 + t70;
t271 = mrSges(6,3) * t97;
t71 = mrSges(6,1) * t307 - t271;
t72 = -mrSges(7,1) * t307 + mrSges(7,2) * t97;
t260 = t71 - t72;
t169 = pkin(5) * t146 - qJ(6) * t148;
t50 = t242 * t98 + t200;
t299 = -qJD(6) * t146 + t169 * t307 - t50;
t248 = t119 * mrSges(5,3);
t244 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t161 + mrSges(6,2) * t97 + t248;
t183 = mrSges(7,1) * t146 - mrSges(7,3) * t148;
t185 = mrSges(6,1) * t146 + mrSges(6,2) * t148;
t89 = t242 * t99;
t48 = t243 * t92 - t89;
t45 = -qJD(3) * pkin(4) - t48;
t21 = -pkin(5) * t161 - t97 * qJ(6) + t45;
t297 = -t21 * t183 - t45 * t185;
t229 = qJD(1) * t147;
t80 = -pkin(4) * t152 - t128 * pkin(8) + t129;
t231 = qJ(4) + t138;
t125 = t231 * t147;
t126 = t231 * t149;
t85 = -t125 * t242 + t126 * t243;
t259 = t146 * t80 + t148 * t85;
t193 = qJD(3) * t231;
t100 = -t147 * t193 + t226;
t101 = -t149 * t193 - t222;
t58 = t100 * t243 + t101 * t242;
t78 = pkin(3) * t227 + pkin(4) * t120 - pkin(8) * t121;
t9 = -qJD(5) * t259 - t146 * t58 + t148 * t78;
t246 = t119 * Ifges(5,4);
t286 = -t246 / 0.2e1 + t312;
t172 = Ifges(7,3) * t146 + t253;
t178 = -Ifges(6,2) * t146 + t256;
t275 = t146 / 0.2e1;
t276 = -t146 / 0.2e1;
t282 = t97 / 0.2e1;
t284 = -t161 / 0.2e1;
t285 = t161 / 0.2e1;
t93 = Ifges(7,5) * t97;
t34 = Ifges(7,6) * t307 - Ifges(7,3) * t161 + t93;
t270 = Ifges(6,4) * t97;
t37 = Ifges(6,2) * t161 + Ifges(6,6) * t307 + t270;
t291 = t172 * t284 + t178 * t285 + t34 * t275 + t37 * t276 - t297 + t295 * t282 + t296 * t307 / 0.2e1;
t290 = t117 * mrSges(5,1) + t19 * mrSges(6,1) - t11 * mrSges(7,1) - t20 * mrSges(6,2) + t12 * mrSges(7,3) + t286;
t289 = t66 / 0.2e1;
t288 = -t67 / 0.2e1;
t287 = t67 / 0.2e1;
t283 = -t97 / 0.2e1;
t281 = t109 / 0.2e1;
t280 = -t307 / 0.2e1;
t277 = -t119 / 0.2e1;
t274 = -t148 / 0.2e1;
t273 = t148 / 0.2e1;
t32 = -t243 * t150 + t242 * t87;
t84 = t243 * t125 + t126 * t242;
t264 = t32 * t84;
t51 = t243 * t98 - t89;
t219 = pkin(3) * t229;
t77 = pkin(4) * t119 - pkin(8) * t118 + t219;
t23 = t146 * t77 + t148 * t51;
t258 = Ifges(4,4) * t147;
t114 = Ifges(5,4) * t118;
t252 = t109 * mrSges(5,3);
t251 = t110 * mrSges(5,3);
t250 = t118 * mrSges(5,3);
t247 = t119 * Ifges(5,1);
t245 = t152 * t32;
t241 = Ifges(4,5) * qJD(3);
t239 = Ifges(4,6) * qJD(3);
t237 = t118 * t146;
t236 = t121 * t146;
t235 = t121 * t148;
t232 = t148 * t118;
t228 = qJD(1) * t149;
t55 = -mrSges(7,1) * t161 - mrSges(7,3) * t97;
t220 = t55 + t244;
t218 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t217 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t216 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t215 = mrSges(4,3) * t229;
t214 = mrSges(4,3) * t228;
t141 = Ifges(4,4) * t228;
t206 = t243 * pkin(3);
t205 = t242 * pkin(3);
t201 = m(4) * t138 + mrSges(4,3);
t194 = t109 * mrSges(5,1) + t110 * mrSges(5,2);
t139 = -t206 - pkin(4);
t189 = -t1 * t146 + t148 * t2;
t187 = -t146 * t3 - t148 * t4;
t132 = t207 * qJD(1);
t57 = t100 * t242 - t243 * t101;
t186 = mrSges(6,1) * t148 - mrSges(6,2) * t146;
t184 = mrSges(7,1) * t148 + mrSges(7,3) * t146;
t177 = Ifges(6,2) * t148 + t257;
t171 = -Ifges(7,3) * t148 + t254;
t170 = t148 * pkin(5) + t146 * qJ(6);
t168 = t11 * t148 - t12 * t146;
t167 = t11 * t146 + t12 * t148;
t166 = -t20 * t146 - t19 * t148;
t165 = t146 * t19 - t148 * t20;
t22 = -t146 * t51 + t148 * t77;
t29 = -t146 * t85 + t148 * t80;
t112 = t130 * t149 + t223;
t8 = t146 * t78 + t148 * t58 + t80 * t224 - t225 * t85;
t151 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t133 = -qJD(3) * mrSges(4,2) + t214;
t131 = qJD(3) * mrSges(4,1) - t215;
t124 = -t170 + t139;
t123 = Ifges(4,1) * t229 + t141 + t241;
t122 = t239 + (t149 * Ifges(4,2) + t258) * qJD(1);
t111 = -t130 * t147 + t143;
t108 = Ifges(7,2) * t109;
t107 = Ifges(6,3) * t109;
t105 = t112 * qJD(3);
t102 = -qJD(3) * mrSges(5,2) + t250;
t86 = -mrSges(5,1) * t118 + mrSges(5,2) * t119;
t82 = t114 + t240 + t247;
t65 = Ifges(7,4) * t66;
t64 = Ifges(6,5) * t66;
t63 = Ifges(6,6) * t67;
t62 = Ifges(7,6) * t67;
t54 = pkin(5) * t97 - qJ(6) * t161;
t44 = t128 * t169 + t84;
t36 = t97 * Ifges(7,4) + Ifges(7,2) * t307 - Ifges(7,6) * t161;
t35 = t97 * Ifges(6,5) + Ifges(6,6) * t161 + Ifges(6,3) * t307;
t27 = pkin(5) * t152 - t29;
t26 = -qJ(6) * t152 + t259;
t25 = mrSges(6,1) * t67 + mrSges(6,2) * t66;
t24 = mrSges(7,1) * t67 - mrSges(7,3) * t66;
t18 = -pkin(5) * t119 - t22;
t17 = qJ(6) * t119 + t23;
t14 = t66 * Ifges(6,4) - t67 * Ifges(6,2) + t109 * Ifges(6,6);
t13 = t66 * Ifges(7,5) + t109 * Ifges(7,6) + t67 * Ifges(7,3);
t10 = t169 * t121 + (qJD(5) * t170 - qJD(6) * t148) * t128 + t57;
t7 = -pkin(5) * t120 - t9;
t6 = t67 * pkin(5) - t66 * qJ(6) - t97 * qJD(6) + t32;
t5 = qJ(6) * t120 - qJD(6) * t152 + t8;
t15 = [t29 * t40 + t27 * t41 + t259 * t42 + t26 * t43 + t44 * t24 + t10 * t55 + t5 * t69 + t8 * t70 + t9 * t71 + t7 * t72 + t84 * t25 + t58 * t102 + t129 * t194 + t244 * t57 + (-t109 * t85 + t110 * t84) * mrSges(5,3) + m(5) * (t33 * t85 - t48 * t57 + t49 * t58 + t264) + m(7) * (t1 * t26 + t10 * t21 + t11 * t7 + t12 * t5 + t2 * t27 + t44 * t6) + m(6) * (t19 * t9 + t20 * t8 + t259 * t3 + t29 * t4 + t45 * t57 + t264) + (-t49 * mrSges(5,3) + t35 / 0.2e1 + t36 / 0.2e1 + t218 * t97 - t216 * t161 + t217 * t307 + t290 + t286) * t120 + (t291 - t48 * mrSges(5,3) + t247 / 0.2e1 + t166 * mrSges(6,3) + t168 * mrSges(7,2) + t114 / 0.2e1 + t82 / 0.2e1 + t301 * t273 - t306) * t121 + (t201 * t104 + (t123 / 0.2e1 - t138 * t131 + 0.3e1 / 0.2e1 * t141 + t241 / 0.2e1 - t201 * t111 + 0.2e1 * t132 * mrSges(4,2)) * qJD(3)) * t149 - (-t33 * mrSges(5,3) + t64 / 0.2e1 - t63 / 0.2e1 + t107 / 0.2e1 + t65 / 0.2e1 + t108 / 0.2e1 + t62 / 0.2e1 - Ifges(5,4) * t110 + t216 * t67 + t218 * t66 + (Ifges(5,2) + t217) * t109 + t151) * t152 + (t6 * t183 + t172 * t287 + t178 * t288 + t13 * t275 + Ifges(5,1) * t110 - Ifges(5,4) * t109 + (mrSges(5,3) + t185) * t32 + t187 * mrSges(6,3) + t189 * mrSges(7,2) + (-mrSges(7,2) * t167 + mrSges(6,3) * t165 + t171 * t285 + t177 * t284 + t184 * t21 + t186 * t45 + t274 * t37 + t280 * t311 + t283 * t310) * qJD(5) + t295 * t289 + t296 * t281 + (qJD(5) * t301 + t14) * t276 + (qJD(5) * t34 + t302) * t273) * t128 + (t201 * t105 + (-t122 / 0.2e1 + t132 * mrSges(4,1) - t138 * t133 - t239 / 0.2e1 - t201 * t112 + (t207 * mrSges(4,1) - 0.3e1 / 0.2e1 * t258 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t149) * qJD(1) + (qJD(1) * (-mrSges(5,1) * t152 + mrSges(5,2) * t128) + m(5) * (t117 + t230) + t86) * pkin(3)) * qJD(3)) * t147; -(t24 + t25 + t251) * t152 + t220 * t120 + (-t147 * t131 + t149 * t133 + (-t147 ^ 2 - t149 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + (-t146 * t260 + t148 * t261 + t102) * t121 + m(5) * (-t120 * t48 + t121 * t49 - t245) + m(7) * (t11 * t236 + t12 * t235 + t120 * t21 - t152 * t6) + m(6) * (t120 * t45 - t19 * t236 + t20 * t235 - t245) + m(4) * (t104 * t147 - t105 * t149 + (-t111 * t147 + t112 * t149) * qJD(3)) + (-t252 + (-t146 * t261 - t148 * t260) * qJD(5) + m(5) * t33 + m(7) * t308 + m(6) * t309 + t305) * t128; (-mrSges(5,1) - t186) * t32 - (-Ifges(4,2) * t229 + t123 + t141) * t228 / 0.2e1 + t221 * Ifges(4,5) * t149 / 0.2e1 + (t224 / 0.2e1 - t232 / 0.2e1) * t301 + (t36 + t35 - t246) * t277 + t171 * t287 + t177 * t288 + t48 * t250 + t14 * t273 + t13 * t274 - t244 * t50 + t299 * t55 - (-Ifges(7,6) * t285 - Ifges(6,6) * t284 - t313 * t283 + (-Ifges(6,3) - Ifges(7,2)) * t280 + t290 + t312) * t119 + (t37 / 0.2e1 - t34 / 0.2e1) * t237 + (-t132 * (mrSges(4,1) * t147 + mrSges(4,2) * t149) - (Ifges(4,1) * t149 - t258) * t229 / 0.2e1) * qJD(1) + ((t242 * t33 - t243 * t32) * pkin(3) - t117 * t219 + t48 * t50 - t49 * t51) * m(5) + (t82 + t114) * t278 + t302 * t275 + (t215 + t131) * t112 + (-t11 * t18 - t12 * t17 + t124 * t6 + t21 * t299) * m(7) + (t139 * t32 - t19 * t22 - t20 * t23 - t45 * t50) * m(6) + t291 * qJD(5) + (t214 - t133) * t111 - t33 * mrSges(5,2) - t17 * t69 - t23 * t70 - t22 * t71 - t18 * t72 - t6 * t184 - t51 * t102 - t104 * mrSges(4,2) - t105 * mrSges(4,1) - Ifges(5,6) * t109 + Ifges(5,5) * t110 + t124 * t24 + t139 * t25 - Ifges(4,6) * t199 / 0.2e1 - t86 * t219 + (t190 * m(7) + t188 * m(6) + (m(6) * t166 + m(7) * t168) * qJD(5) - t260 * t224 - t261 * t225 + t305) * (t205 + pkin(8)) + t122 * t229 / 0.2e1 + (Ifges(5,1) * t277 + t172 * t285 + t178 * t284 + t296 * t280 + t295 * t283 + t297 + t306) * t118 + t49 * t248 + (-t11 * t232 + t12 * t237 + t308) * mrSges(7,2) + (t19 * t232 + t20 * t237 + t309) * mrSges(6,3) + t310 * t289 - t206 * t251 + t311 * t281 - t205 * t252; -t118 * t102 - t220 * t119 + (t261 * t307 - t263) * t148 + (-t260 * t307 + t262) * t146 + t194 + (-t119 * t21 + t167 * t307 - t189) * m(7) + (-t119 * t45 - t165 * t307 - t187) * m(6) + (-t118 * t49 + t119 * t48 + t191) * m(5); (Ifges(7,3) * t97 + t269) * t285 + t37 * t282 + (-t11 * t161 + t12 * t97) * mrSges(7,2) + t107 + t108 + t151 + t65 + t64 - t63 + t62 - pkin(5) * t41 + qJ(6) * t43 - t54 * t55 + qJD(6) * t69 - t21 * (mrSges(7,1) * t97 - mrSges(7,3) * t161) - t45 * (mrSges(6,1) * t97 + mrSges(6,2) * t161) + (t260 + t271) * t20 + (-t261 + t272) * t19 + (t161 * t313 - t315 * t97) * t280 + (-pkin(5) * t2 + qJ(6) * t1 - t11 * t20 + t12 * t298 - t21 * t54) * m(7) + (-Ifges(6,2) * t97 + t301 + t94) * t284 + (t161 * t314 - t270 + t34 + t93) * t283; -t307 * t69 + t97 * t55 + 0.2e1 * (t2 / 0.2e1 + t12 * t280 + t21 * t282) * m(7) + t41;];
tauc  = t15(:);
