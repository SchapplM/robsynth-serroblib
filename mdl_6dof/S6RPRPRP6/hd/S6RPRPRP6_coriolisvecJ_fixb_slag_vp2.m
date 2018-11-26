% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2018-11-23 16:00
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:59:48
% EndTime: 2018-11-23 15:59:53
% DurationCPUTime: 5.70s
% Computational Cost: add. (5035->442), mult. (13113->547), div. (0->0), fcn. (9190->6), ass. (0->208)
t304 = Ifges(6,4) + Ifges(7,4);
t305 = Ifges(6,1) + Ifges(7,1);
t292 = Ifges(6,5) + Ifges(7,5);
t303 = Ifges(6,2) + Ifges(7,2);
t289 = Ifges(6,6) + Ifges(7,6);
t161 = cos(qJ(5));
t159 = sin(qJ(5));
t158 = cos(pkin(9));
t250 = cos(qJ(3));
t208 = t250 * t158;
t203 = qJD(1) * t208;
t157 = sin(pkin(9));
t160 = sin(qJ(3));
t230 = t157 * t160;
t130 = qJD(1) * t230 - t203;
t137 = t157 * t250 + t160 * t158;
t131 = t137 * qJD(1);
t210 = -pkin(2) * t158 - pkin(1);
t140 = qJD(1) * t210 + qJD(2);
t167 = -qJ(4) * t131 + t140;
t261 = pkin(3) + pkin(8);
t52 = t130 * t261 + t167;
t245 = pkin(7) + qJ(2);
t143 = t245 * t157;
t138 = qJD(1) * t143;
t144 = t245 * t158;
t139 = qJD(1) * t144;
t103 = t250 * t138 + t139 * t160;
t175 = pkin(4) * t131 + t103;
t59 = -qJD(3) * t261 + qJD(4) + t175;
t15 = -t159 * t52 + t161 * t59;
t16 = t159 * t59 + t161 * t52;
t181 = t15 * t159 - t16 * t161;
t194 = mrSges(7,1) * t161 - mrSges(7,2) * t159;
t196 = mrSges(6,1) * t161 - mrSges(6,2) * t159;
t109 = -qJD(3) * t159 + t130 * t161;
t11 = qJ(6) * t109 + t16;
t110 = qJD(3) * t161 + t130 * t159;
t10 = -qJ(6) * t110 + t15;
t126 = qJD(5) + t131;
t9 = pkin(5) * t126 + t10;
t198 = t11 * t161 - t9 * t159;
t254 = -t159 / 0.2e1;
t301 = t304 * t109;
t274 = t305 * t110 + t292 * t126 + t301;
t300 = t304 * t110;
t275 = t303 * t109 + t289 * t126 + t300;
t156 = qJD(3) * qJ(4);
t104 = -t160 * t138 + t250 * t139;
t79 = -pkin(4) * t130 + t104;
t65 = t156 + t79;
t31 = -pkin(5) * t109 + qJD(6) + t65;
t302 = t31 * t194 + t65 * t196 + t274 * t254 - t275 * t161 / 0.2e1 + t181 * mrSges(6,3) - t198 * mrSges(7,3);
t294 = mrSges(4,3) + mrSges(5,1);
t293 = Ifges(4,4) + Ifges(5,6);
t299 = t304 * t161;
t298 = t304 * t159;
t214 = Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1;
t255 = t126 / 0.2e1;
t257 = t110 / 0.2e1;
t259 = t109 / 0.2e1;
t270 = t305 * t159 + t299;
t271 = t303 * t161 + t298;
t272 = t292 * t159 + t289 * t161;
t276 = qJD(3) / 0.2e1;
t277 = -qJD(3) / 0.2e1;
t291 = Ifges(4,2) + Ifges(5,3);
t295 = t131 / 0.2e1;
t296 = -t130 / 0.2e1;
t83 = pkin(3) * t130 + t167;
t99 = -t156 - t104;
t297 = -(Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1) * qJD(3) + t214 * t131 - t140 * mrSges(4,1) - t99 * mrSges(5,1) + t83 * mrSges(5,2) + t104 * mrSges(4,3) - Ifges(5,5) * t276 - Ifges(4,6) * t277 - t255 * t272 - t257 * t270 - t259 * t271 + t291 * t296 + t293 * t295 + t302;
t281 = -qJD(4) - t103;
t247 = mrSges(5,2) - mrSges(4,1);
t290 = Ifges(5,2) + Ifges(4,1);
t207 = qJD(3) * t230;
t119 = qJD(1) * t207 - qJD(3) * t203;
t133 = t137 * qJD(3);
t120 = qJD(1) * t133;
t80 = qJD(5) * t109 + t120 * t159;
t81 = -qJD(5) * t110 + t120 * t161;
t288 = -t289 * t119 + t303 * t81 + t304 * t80;
t287 = -t292 * t119 + t304 * t81 + t305 * t80;
t221 = qJD(6) * t161;
t223 = qJD(5) * t159;
t229 = qJ(6) + t261;
t67 = t161 * t79;
t232 = qJ(4) * t130;
t68 = t131 * t261 + t232;
t286 = pkin(5) * t130 - t67 - (-qJ(6) * t131 - t68) * t159 + t223 * t229 - t221;
t142 = t229 * t161;
t231 = qJ(6) * t161;
t24 = t159 * t79 + t161 * t68;
t285 = -qJD(5) * t142 - qJD(6) * t159 - t131 * t231 - t24;
t148 = qJD(2) * t208;
t206 = qJD(3) * t250;
t220 = qJD(1) * qJD(2);
t63 = -t160 * (qJD(3) * t139 + t157 * t220) + qJD(1) * t148 - t138 * t206;
t60 = -qJD(3) * qJD(4) - t63;
t209 = -pkin(5) * t161 - pkin(4);
t222 = qJD(5) * t161;
t284 = pkin(5) * t222 - t131 * t209 - t281;
t283 = -t303 * t159 + t299;
t282 = t305 * t161 - t298;
t179 = qJ(4) * t119 - qJD(4) * t131;
t32 = t120 * t261 + t179;
t169 = t137 * qJD(2);
t64 = qJD(1) * t169 + qJD(3) * t104;
t36 = -t119 * pkin(4) + t64;
t5 = t159 * t36 + t161 * t32 + t59 * t222 - t223 * t52;
t6 = -qJD(5) * t16 - t159 * t32 + t161 * t36;
t269 = t159 * t5 + t161 * t6;
t268 = t293 * t296;
t267 = (m(3) * qJ(2) + mrSges(3,3)) * (t157 ^ 2 + t158 ^ 2);
t88 = (qJD(2) * t157 + qJD(3) * t144) * t160 + t143 * t206 - t148;
t211 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t213 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t216 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t98 = -qJD(3) * pkin(3) - t281;
t264 = -(Ifges(5,4) / 0.2e1 - Ifges(4,5) / 0.2e1) * qJD(3) + t213 * t109 + t216 * t110 + t211 * t126 + t98 * mrSges(5,1) + t15 * mrSges(6,1) + t9 * mrSges(7,1) + t140 * mrSges(4,2) - t16 * mrSges(6,2) - t11 * mrSges(7,2) + t103 * mrSges(4,3) - t83 * mrSges(5,3) + Ifges(5,4) * t277 + Ifges(4,5) * t276 + t268 + t290 * t295 + t289 * t259 + t292 * t257 + (Ifges(7,3) + Ifges(6,3)) * t255;
t263 = t80 / 0.2e1;
t262 = t81 / 0.2e1;
t260 = -t109 / 0.2e1;
t258 = -t110 / 0.2e1;
t256 = -t126 / 0.2e1;
t251 = t161 / 0.2e1;
t43 = -mrSges(7,1) * t119 - mrSges(7,3) * t80;
t44 = -mrSges(6,1) * t119 - mrSges(6,3) * t80;
t244 = -t43 - t44;
t45 = mrSges(7,2) * t119 + mrSges(7,3) * t81;
t46 = mrSges(6,2) * t119 + mrSges(6,3) * t81;
t243 = t45 + t46;
t136 = -t208 + t230;
t176 = -qJ(4) * t137 + t210;
t77 = t136 * t261 + t176;
t105 = t250 * t143 + t144 * t160;
t92 = pkin(4) * t137 + t105;
t29 = t159 * t92 + t161 * t77;
t84 = -mrSges(7,2) * t126 + mrSges(7,3) * t109;
t85 = -mrSges(6,2) * t126 + mrSges(6,3) * t109;
t242 = t84 + t85;
t86 = mrSges(7,1) * t126 - mrSges(7,3) * t110;
t87 = mrSges(6,1) * t126 - mrSges(6,3) * t110;
t241 = -t86 - t87;
t236 = t105 * t64;
t113 = t130 * mrSges(5,1) - qJD(3) * mrSges(5,3);
t70 = -mrSges(6,1) * t109 + mrSges(6,2) * t110;
t233 = t70 - t113;
t111 = -qJD(3) * mrSges(4,2) - t130 * mrSges(4,3);
t228 = t111 - t113;
t227 = -qJD(3) * t247 - t131 * t294;
t69 = -mrSges(7,1) * t109 + mrSges(7,2) * t110;
t219 = -t69 - t233;
t215 = -Ifges(5,2) / 0.2e1 - Ifges(4,1) / 0.2e1;
t212 = Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1;
t25 = -t81 * mrSges(7,1) + t80 * mrSges(7,2);
t205 = -qJ(6) * t136 - t77;
t1 = -pkin(5) * t119 - qJ(6) * t80 - qJD(6) * t110 + t6;
t2 = qJ(6) * t81 + qJD(6) * t109 + t5;
t202 = -t1 * t161 - t2 * t159;
t201 = -t1 * t159 + t161 * t2;
t200 = -t159 * t6 + t161 * t5;
t197 = t11 * t159 + t161 * t9;
t195 = mrSges(6,1) * t159 + mrSges(6,2) * t161;
t193 = mrSges(7,1) * t159 + mrSges(7,2) * t161;
t182 = t15 * t161 + t159 * t16;
t132 = -t158 * t206 + t207;
t178 = qJ(4) * t132 - qJD(4) * t137;
t49 = t133 * t261 + t178;
t106 = -t160 * t143 + t144 * t250;
t89 = qJD(3) * t106 + t169;
t54 = -t132 * pkin(4) + t89;
t7 = t159 * t54 + t161 * t49 + t92 * t222 - t223 * t77;
t172 = -t161 * t133 + t136 * t223;
t170 = -t159 * t242 + t161 * t241;
t168 = t6 * mrSges(6,1) + t1 * mrSges(7,1) - t5 * mrSges(6,2) - t2 * mrSges(7,2);
t53 = -pkin(4) * t133 - t88;
t33 = -pkin(4) * t120 - t60;
t152 = pkin(5) * t159 + qJ(4);
t141 = t229 * t159;
t118 = Ifges(6,3) * t119;
t117 = Ifges(7,3) * t119;
t116 = t119 * mrSges(4,2);
t115 = t119 * mrSges(5,3);
t102 = pkin(3) * t136 + t176;
t101 = -mrSges(5,2) * t130 - mrSges(5,3) * t131;
t100 = pkin(3) * t131 + t232;
t93 = -t136 * pkin(4) + t106;
t91 = t161 * t92;
t76 = Ifges(6,5) * t80;
t75 = Ifges(7,5) * t80;
t74 = Ifges(6,6) * t81;
t73 = Ifges(7,6) * t81;
t71 = pkin(3) * t133 + t178;
t58 = t136 * t209 + t106;
t57 = pkin(3) * t120 + t179;
t51 = t161 * t54;
t28 = -t159 * t77 + t91;
t27 = pkin(5) * t172 + t53;
t26 = -mrSges(6,1) * t81 + mrSges(6,2) * t80;
t23 = -t159 * t68 + t67;
t22 = t136 * t231 + t29;
t17 = -pkin(5) * t81 + t33;
t13 = pkin(5) * t137 + t159 * t205 + t91;
t8 = -qJD(5) * t29 - t159 * t49 + t51;
t4 = -qJ(6) * t172 + t136 * t221 + t7;
t3 = -pkin(5) * t132 + t51 + t205 * t222 + (-qJ(6) * t133 - qJD(5) * t92 - qJD(6) * t136 - t49) * t159;
t12 = [(t130 * t214 + t131 * t215 - t264) * t132 + (t130 * t212 - t297) * t133 + (-t17 * t194 - t33 * t196 + t60 * mrSges(5,1) - t63 * mrSges(4,3) - t57 * mrSges(5,2) + t291 * t120 + (-t159 * t216 - t161 * t213 + t293) * t119 + t201 * mrSges(7,3) + t200 * mrSges(6,3) + (-t182 * mrSges(6,3) - t197 * mrSges(7,3) + t31 * t193 + t65 * t195 + t283 * t259 + t282 * t257 + (-t159 * t289 + t161 * t292) * t255 + t275 * t254) * qJD(5) + t270 * t263 + t271 * t262 + t287 * t159 / 0.2e1 + (t274 * qJD(5) + t288) * t251) * t136 + (-t57 * mrSges(5,3) - t117 / 0.2e1 - t118 / 0.2e1 + t76 / 0.2e1 + t75 / 0.2e1 + t74 / 0.2e1 + t73 / 0.2e1 + t213 * t81 + t216 * t80 + t294 * t64 - t293 * t120 + (-t211 - t290) * t119 + t168) * t137 + 0.2e1 * t267 * t220 + t294 * (-t105 * t119 - t106 * t120) - t227 * t89 - t228 * t88 + m(4) * (t103 * t89 - t104 * t88 + t106 * t63 + t236) + m(5) * (t102 * t57 - t106 * t60 + t71 * t83 + t88 * t99 + t89 * t98 + t236) + t210 * (t120 * mrSges(4,1) - t116) + m(7) * (t1 * t13 + t11 * t4 + t17 * t58 + t2 * t22 + t27 * t31 + t3 * t9) + m(6) * (t15 * t8 + t16 * t7 + t28 * t6 + t29 * t5 + t33 * t93 + t53 * t65) + t13 * t43 + t28 * t44 + t22 * t45 + t29 * t46 + t58 * t25 + t27 * t69 + t53 * t70 + t4 * t84 + t7 * t85 + t3 * t86 + t8 * t87 + t93 * t26 + t71 * t101 + t102 * (-t120 * mrSges(5,2) + t115); t115 - t116 + t243 * t161 + t244 * t159 - t247 * t120 + (t111 - t219) * t130 + t170 * qJD(5) + (t170 + t227) * t131 - m(4) * (t103 * t131 - t104 * t130) - t267 * qJD(1) ^ 2 + (-t126 * t197 + t130 * t31 + t201) * m(7) + (-t126 * t182 + t130 * t65 + t200) * m(6) + (-t130 * t99 - t131 * t98 + t57) * m(5); -t269 * mrSges(6,3) + ((-t212 - t215) * t130 + t297) * t131 + (-(-m(6) * t181 - t159 * t87 + t161 * t85) * t261 + t271 * t260 + t270 * t258 + t272 * t256 + t302) * qJD(5) + t284 * t69 + t285 * t84 + (-t1 * t142 + t11 * t285 - t141 * t2 + t152 * t17 + t284 * t31 + t286 * t9) * m(7) + t286 * t86 - m(6) * (t15 * t23 + t16 * t24 - t175 * t65) + t175 * t70 + m(6) * (qJ(4) * t33 + qJD(4) * t65 - t261 * t269) - (t159 * t46 + t161 * t44) * t261 + t227 * t104 + t228 * t103 + (pkin(3) * mrSges(5,1) + t159 * t213 - t161 * t216 + Ifges(5,4) - Ifges(4,5)) * t119 + t202 * mrSges(7,3) + (-pkin(3) * t64 - qJ(4) * t60 - t100 * t83 - t104 * t98 + t281 * t99) * m(5) + t282 * t263 + t283 * t262 + (t264 + t268) * t130 + t247 * t64 + t287 * t251 + t288 * t254 + (-qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6)) * t120 + t233 * qJD(4) + t17 * t193 + t33 * t195 + qJ(4) * t26 - t60 * mrSges(5,3) - t63 * mrSges(4,2) - t24 * t85 - t23 * t87 - t100 * t101 - t141 * t45 - t142 * t43 + t152 * t25; -t119 * mrSges(5,1) + t131 * t101 + t219 * qJD(3) + (t126 * t242 - t244) * t161 + (t126 * t241 + t243) * t159 + (-qJD(3) * t31 + t126 * t198 - t202) * m(7) + (-qJD(3) * t65 - t126 * t181 + t269) * m(6) + (qJD(3) * t99 + t131 * t83 + t64) * m(5); t168 + (-(t10 - t9) * t11 + (-t110 * t31 + t1) * pkin(5)) * m(7) + (-t110 * t69 + t43) * pkin(5) - t117 - t118 + t76 + t75 + t74 + t73 + (t109 * t9 + t11 * t110) * mrSges(7,3) + (t109 * t15 + t110 * t16) * mrSges(6,3) - t10 * t84 - t15 * t85 + t11 * t86 + t16 * t87 - t31 * (mrSges(7,1) * t110 + mrSges(7,2) * t109) - t65 * (mrSges(6,1) * t110 + mrSges(6,2) * t109) + (t305 * t109 - t300) * t258 + t275 * t257 + (t109 * t292 - t110 * t289) * t256 + (-t303 * t110 + t274 + t301) * t260; -t109 * t84 + t110 * t86 + 0.2e1 * (t17 / 0.2e1 + t11 * t260 + t9 * t257) * m(7) + t25;];
tauc  = t12(:);
