% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2018-11-23 15:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRPPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:10:39
% EndTime: 2018-11-23 15:10:44
% DurationCPUTime: 5.01s
% Computational Cost: add. (4016->512), mult. (9987->731), div. (0->0), fcn. (6584->10), ass. (0->251)
t297 = pkin(4) + pkin(8);
t312 = mrSges(5,1) + mrSges(4,3);
t322 = mrSges(5,2) - mrSges(4,1);
t189 = cos(pkin(6));
t195 = cos(qJ(3));
t264 = t189 * t195;
t167 = qJD(1) * t264;
t192 = sin(qJ(3));
t193 = sin(qJ(2));
t187 = sin(pkin(6));
t255 = qJD(1) * t187;
t235 = t193 * t255;
t152 = qJD(2) * pkin(8) + t235;
t226 = pkin(4) * qJD(2) + t152;
t214 = t226 * t192;
t96 = t167 - t214;
t321 = qJD(4) - t96;
t186 = sin(pkin(11));
t188 = cos(pkin(11));
t254 = qJD(1) * t192;
t166 = t189 * t254;
t196 = cos(qJ(2));
t253 = qJD(2) * t187;
t230 = qJD(1) * t253;
t224 = t196 * t230;
t52 = t192 * t224 + (t195 * t226 - qJD(5) + t166) * qJD(3);
t213 = -qJ(4) * t195 + qJ(5) * t192;
t247 = qJD(4) * t192;
t200 = qJD(3) * t213 - qJD(5) * t195 - t247;
t244 = qJD(2) * qJD(3);
t229 = t192 * t244;
t256 = pkin(3) * t229 + t193 * t230;
t60 = qJD(2) * t200 + t256;
t16 = -t186 * t60 + t188 * t52;
t269 = t186 * t192;
t209 = pkin(5) * t195 - pkin(9) * t269;
t204 = t209 * qJD(3);
t14 = qJD(2) * t204 + t16;
t17 = t186 * t52 + t188 * t60;
t223 = t188 * t229;
t15 = pkin(9) * t223 + t17;
t191 = sin(qJ(6));
t194 = cos(qJ(6));
t245 = t188 * qJD(3);
t246 = t186 * qJD(2);
t138 = t195 * t246 - t245;
t190 = -pkin(3) - qJ(5);
t74 = qJD(3) * t190 + t321;
t227 = -qJ(4) * t192 - pkin(2);
t136 = t190 * t195 + t227;
t234 = t196 * t255;
t95 = qJD(2) * t136 - t234;
t25 = -t186 * t95 + t188 * t74;
t252 = qJD(2) * t192;
t18 = pkin(5) * t252 + pkin(9) * t138 + t25;
t250 = qJD(2) * t195;
t139 = -t186 * qJD(3) - t188 * t250;
t26 = t186 * t74 + t188 * t95;
t21 = pkin(9) * t139 + t26;
t5 = t18 * t194 - t191 * t21;
t1 = qJD(6) * t5 + t14 * t191 + t15 * t194;
t6 = t18 * t191 + t194 * t21;
t2 = -qJD(6) * t6 + t14 * t194 - t15 * t191;
t211 = t194 * t186 + t191 * t188;
t205 = t211 * t192;
t112 = qJD(2) * t205;
t125 = t211 * qJD(6);
t261 = -t125 - t112;
t231 = t188 * t252;
t111 = -t191 * t192 * t246 + t194 * t231;
t309 = -t186 * t191 + t188 * t194;
t124 = t309 * qJD(6);
t262 = t124 + t111;
t320 = -t1 * t211 - t2 * t309 - t261 * t5 - t262 * t6;
t178 = pkin(3) * t252;
t118 = qJD(2) * t213 + t178;
t106 = t195 * t152 + t166;
t97 = pkin(4) * t250 + t106;
t44 = -t118 * t186 + t188 * t97;
t30 = qJD(2) * t209 + t44;
t45 = t188 * t118 + t186 * t97;
t33 = pkin(9) * t231 + t45;
t284 = -pkin(9) + t190;
t149 = t284 * t186;
t150 = t284 * t188;
t89 = t149 * t194 + t150 * t191;
t318 = -qJD(5) * t309 - qJD(6) * t89 + t191 * t33 - t194 * t30;
t315 = -qJD(2) / 0.2e1;
t314 = -qJD(3) / 0.2e1;
t313 = qJD(3) / 0.2e1;
t88 = -t149 * t191 + t150 * t194;
t311 = -qJD(5) * t211 + qJD(6) * t88 - t191 * t30 - t194 * t33;
t225 = t138 * t191 + t194 * t139;
t76 = t138 * t194 - t139 * t191;
t27 = -mrSges(7,1) * t225 - mrSges(7,2) * t76;
t87 = -mrSges(6,1) * t139 - mrSges(6,2) * t138;
t310 = t27 + t87;
t105 = t152 * t192 - t167;
t308 = -qJD(4) - t105;
t202 = qJD(3) * t205;
t38 = qJD(2) * t202 + qJD(6) * t225;
t249 = qJD(3) * t192;
t201 = t309 * t249;
t39 = qJD(2) * t201 + qJD(6) * t76;
t307 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t38 + Ifges(7,6) * t39;
t185 = qJD(3) * qJ(4);
t102 = -t185 - t106;
t154 = -pkin(3) * t195 + t227;
t107 = qJD(2) * t154 - t234;
t277 = qJD(2) * pkin(2);
t153 = -t234 - t277;
t215 = t25 * t186 - t26 * t188;
t281 = Ifges(6,4) * t186;
t217 = Ifges(6,2) * t188 + t281;
t280 = Ifges(6,4) * t188;
t282 = Ifges(6,1) * t186;
t218 = t280 + t282;
t283 = mrSges(6,2) * t186;
t219 = mrSges(6,1) * t188 - t283;
t288 = t188 / 0.2e1;
t289 = t186 / 0.2e1;
t80 = qJD(5) + t185 + t97;
t306 = -t215 * mrSges(6,3) + t102 * mrSges(5,1) + t153 * mrSges(4,1) - t80 * t219 + Ifges(4,6) * t314 + (Ifges(4,4) * t192 + t195 * Ifges(4,2)) * t315 + Ifges(5,5) * t313 + (-Ifges(5,6) * t192 - t195 * Ifges(5,3)) * qJD(2) / 0.2e1 - t106 * mrSges(4,3) - t107 * mrSges(5,2) - t138 * t218 / 0.2e1 + t139 * t217 / 0.2e1 + (-t138 * Ifges(6,1) + t139 * Ifges(6,4) + Ifges(6,5) * t252) * t289 + (-t138 * Ifges(6,4) + t139 * Ifges(6,2) + Ifges(6,6) * t252) * t288;
t305 = -m(5) / 0.2e1;
t304 = t38 / 0.2e1;
t303 = t39 / 0.2e1;
t302 = -t225 / 0.2e1;
t301 = t225 / 0.2e1;
t300 = t76 / 0.2e1;
t299 = -t76 / 0.2e1;
t298 = -(Ifges(6,6) * t195 + t192 * t217) * t244 / 0.2e1;
t115 = t309 * t195;
t296 = -t115 / 0.2e1;
t116 = t211 * t195;
t295 = -t116 / 0.2e1;
t294 = -t211 / 0.2e1;
t293 = t309 / 0.2e1;
t172 = qJD(6) + t252;
t292 = -t172 / 0.2e1;
t291 = t172 / 0.2e1;
t290 = -t186 / 0.2e1;
t287 = Ifges(7,4) * t76;
t279 = Ifges(6,5) * t186;
t278 = Ifges(6,6) * t188;
t268 = t187 * t193;
t238 = t192 * t268;
t126 = t238 - t264;
t232 = t196 * t253;
t248 = qJD(3) * t195;
t69 = t152 * t248 + (qJD(3) * t189 + t232) * t254;
t276 = t126 * t69;
t113 = -mrSges(6,1) * t223 + t229 * t283;
t13 = -t39 * mrSges(7,1) + t38 * mrSges(7,2);
t272 = t13 + t113;
t180 = pkin(3) * t249;
t100 = t180 + t200;
t148 = t297 * t248;
t54 = t188 * t100 + t186 * t148;
t267 = t187 * t196;
t265 = t188 * t195;
t263 = t192 * t196;
t162 = t297 * t192;
t85 = t188 * t136 + t186 * t162;
t144 = (mrSges(5,2) * t195 - mrSges(5,3) * t192) * qJD(2);
t260 = t144 + (-mrSges(4,1) * t195 + mrSges(4,2) * t192) * qJD(2);
t259 = qJD(3) * t167 + t195 * t224;
t258 = -qJD(3) * t322 - t252 * t312;
t159 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t250;
t160 = -mrSges(5,1) * t250 - qJD(3) * mrSges(5,3);
t257 = t159 - t160;
t163 = t297 * t195;
t251 = qJD(2) * t193;
t243 = pkin(8) * t69 / 0.2e1;
t242 = -t160 + t310;
t241 = 0.3e1 / 0.2e1 * Ifges(4,4) + 0.3e1 / 0.2e1 * Ifges(5,6);
t240 = -Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t239 = -Ifges(4,6) / 0.2e1 + Ifges(5,5) / 0.2e1;
t184 = qJD(3) * qJD(4);
t237 = t184 + t259;
t236 = -pkin(5) * t188 - pkin(4);
t233 = t187 * t251;
t228 = t195 * t244;
t53 = -t100 * t186 + t188 * t148;
t222 = -t234 / 0.2e1;
t216 = t16 * t188 + t17 * t186;
t143 = t188 * t162;
t61 = pkin(5) * t192 + t143 + (pkin(9) * t195 - t136) * t186;
t67 = -pkin(9) * t265 + t85;
t19 = -t191 * t67 + t194 * t61;
t20 = t191 * t61 + t194 * t67;
t91 = t126 * t188 + t186 * t267;
t92 = t126 * t186 - t188 * t267;
t31 = -t191 * t92 + t194 * t91;
t32 = t191 * t91 + t194 * t92;
t68 = -t152 * t249 + t259;
t208 = -t279 / 0.2e1 - t278 / 0.2e1;
t127 = t189 * t192 + t195 * t268;
t207 = -qJ(4) * t248 - t247;
t206 = -m(4) * t106 + m(5) * t102 - t257;
t203 = (qJD(2) * t236 - t152) * t192;
t101 = -qJD(3) * pkin(3) - t308;
t177 = Ifges(4,4) * t250;
t199 = t101 * mrSges(5,1) + t105 * mrSges(4,3) + t153 * mrSges(4,2) + t25 * mrSges(6,1) + t5 * mrSges(7,1) + Ifges(4,5) * t313 + t177 / 0.2e1 + Ifges(5,4) * t314 + (-Ifges(5,2) * t192 - Ifges(5,6) * t195) * t315 + t172 * Ifges(7,3) - t76 * Ifges(7,5) + t225 * Ifges(7,6) + t139 * Ifges(6,6) - t138 * Ifges(6,5) - t107 * mrSges(5,3) - t26 * mrSges(6,2) - t6 * mrSges(7,2) + (Ifges(4,1) + Ifges(6,3)) * t252 / 0.2e1;
t197 = qJD(2) ^ 2;
t173 = pkin(5) * t186 + qJ(4);
t170 = Ifges(7,3) * t228;
t147 = t297 * t249;
t146 = -qJ(4) * t250 + t178;
t135 = (mrSges(4,1) * t192 + mrSges(4,2) * t195) * t244;
t134 = (-mrSges(5,2) * t192 - mrSges(5,3) * t195) * t244;
t123 = pkin(5) * t265 + t163;
t121 = t180 + t207;
t120 = (mrSges(6,3) * t188 * t192 - mrSges(6,2) * t195) * t244;
t119 = (mrSges(6,1) * t195 - mrSges(6,3) * t269) * t244;
t114 = (-pkin(8) + t236) * t249;
t110 = mrSges(6,1) * t252 + mrSges(6,3) * t138;
t109 = -mrSges(6,2) * t252 + mrSges(6,3) * t139;
t104 = (t186 * t263 + t188 * t193) * t255;
t103 = (-t186 * t193 + t188 * t263) * t255;
t99 = (Ifges(6,5) * t195 + t192 * t218) * t244;
t94 = qJD(3) * t127 + t192 * t232;
t93 = qJD(3) * t238 - t189 * t248 - t195 * t232;
t90 = qJD(2) * t207 + t256;
t84 = -t136 * t186 + t143;
t75 = Ifges(7,4) * t225;
t70 = t167 + t203;
t64 = -t184 - t68;
t63 = t125 * t195 + t201;
t62 = -qJD(6) * t115 + t202;
t59 = t186 * t94 + t188 * t233;
t58 = -t186 * t233 + t188 * t94;
t57 = mrSges(7,1) * t172 + mrSges(7,3) * t76;
t56 = -mrSges(7,2) * t172 + mrSges(7,3) * t225;
t51 = -qJD(3) * t214 + t237;
t50 = -pkin(5) * t139 + t80;
t43 = pkin(9) * t192 * t245 + t54;
t42 = t103 * t191 + t104 * t194;
t41 = t103 * t194 - t104 * t191;
t40 = qJD(3) * t203 + t237;
t37 = t204 + t53;
t29 = -mrSges(7,2) * t228 + mrSges(7,3) * t39;
t28 = mrSges(7,1) * t228 - mrSges(7,3) * t38;
t24 = -Ifges(7,1) * t76 + Ifges(7,5) * t172 + t75;
t23 = Ifges(7,2) * t225 + Ifges(7,6) * t172 - t287;
t12 = t38 * Ifges(7,1) + t39 * Ifges(7,4) + Ifges(7,5) * t228;
t11 = t38 * Ifges(7,4) + t39 * Ifges(7,2) + Ifges(7,6) * t228;
t8 = -qJD(6) * t32 - t191 * t59 + t194 * t58;
t7 = qJD(6) * t31 + t191 * t58 + t194 * t59;
t4 = -qJD(6) * t20 - t191 * t43 + t194 * t37;
t3 = qJD(6) * t19 + t191 * t37 + t194 * t43;
t9 = [t59 * t109 + t58 * t110 + t91 * t119 + t92 * t120 + t31 * t28 + t32 * t29 + t7 * t56 + t8 * t57 - t258 * t94 + t272 * t127 + (-t159 - t242) * t93 + (-t197 * t193 * mrSges(3,1) + (-mrSges(3,2) * t197 - t134 - t135) * t196) * t187 + (t260 * t268 + t312 * qJD(3) * (t126 * t195 - t127 * t192)) * qJD(2) + m(6) * (t127 * t51 + t16 * t91 + t17 * t92 + t25 * t58 + t26 * t59 - t80 * t93) + m(7) * (t1 * t32 + t127 * t40 + t2 * t31 + t5 * t8 - t50 * t93 + t6 * t7) + m(5) * (t101 * t94 + t102 * t93 + t276 - t127 * t64 + (t107 * t251 - t196 * t90) * t187) + m(4) * (t105 * t94 - t106 * t93 + t276 + t127 * t68 + (t153 - t234) * t233); (-t1 * t115 + t116 * t2 - t5 * t62 + t6 * t63) * mrSges(7,3) + t40 * (mrSges(7,1) * t115 - mrSges(7,2) * t116) + (-Ifges(7,4) * t116 - Ifges(7,2) * t115) * t303 + (-Ifges(7,1) * t116 - Ifges(7,4) * t115) * t304 - t260 * t235 + m(7) * (t1 * t20 + t114 * t50 + t123 * t40 + t19 * t2 + t3 * t6 + t4 * t5) + m(6) * (-t147 * t80 + t16 * t84 + t163 * t51 + t17 * t85 + t25 * t53 + t26 * t54) + m(5) * (t107 * t121 + t154 * t90) - m(6) * (t103 * t25 + t104 * t26) - m(7) * (t41 * t5 + t42 * t6) + 0.2e1 * (t107 * t305 + (-t153 / 0.2e1 - t277 / 0.2e1) * m(4)) * t235 + (t54 - t104) * t109 + (t53 - t103) * t110 + t19 * t28 + t20 * t29 + (t68 * mrSges(4,3) - t64 * mrSges(5,1) + t99 * t290 + t188 * t298 + t51 * t219 + t90 * mrSges(5,2) + (t16 * t186 - t17 * t188) * mrSges(6,3) + (-mrSges(4,1) * t251 + (-m(6) * t80 - m(7) * t50 + t206 - t310) * t196) * t255 + (t240 * qJD(3) + ((t208 + t241) * t195 + Ifges(7,5) * t295 + Ifges(7,6) * t296) * qJD(2) + t199) * qJD(3) + (0.2e1 * t64 * t305 + m(4) * t68 + (m(4) * t105 + m(5) * t101 - t258) * qJD(3)) * pkin(8)) * t195 + (-t90 * mrSges(5,3) + t16 * mrSges(6,1) - t17 * mrSges(6,2) + t170 / 0.2e1 + t312 * t69 + (mrSges(4,2) * t251 + t196 * t258) * t255 + 0.2e1 * (t101 * t222 + t243) * m(5) + 0.2e1 * (t105 * t222 + t243) * m(4) + (t206 * pkin(8) + t239 * qJD(3) + ((0.3e1 / 0.2e1 * t279 + 0.3e1 / 0.2e1 * t278 - t241) * t192 + (0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(6,3) - 0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2) - t188 ^ 2 * Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1 + (-t280 - t282 / 0.2e1) * t186) * t195) * qJD(2) + t306) * qJD(3) + t307) * t192 + t62 * t24 / 0.2e1 + t63 * t23 / 0.2e1 + t50 * (-mrSges(7,1) * t63 + mrSges(7,2) * t62) + (t3 - t42) * t56 + (Ifges(7,5) * t62 + Ifges(7,6) * t63) * t291 + t114 * t27 + t84 * t119 + t85 * t120 + t123 * t13 - pkin(2) * t135 + t121 * t144 - t147 * t87 + t154 * t134 + t163 * t113 + (t4 - t41) * t57 + t12 * t295 + t11 * t296 + (Ifges(7,1) * t62 + Ifges(7,4) * t63) * t299 + (Ifges(7,4) * t62 + Ifges(7,2) * t63) * t301; (qJ(4) * t51 + t216 * t190 + (-t186 * t26 - t188 * t25) * qJD(5) - t25 * t44 - t26 * t45 + t321 * t80) * m(6) + (Ifges(7,1) * t309 - Ifges(7,4) * t211) * t304 + t40 * (mrSges(7,1) * t211 + mrSges(7,2) * t309) + (Ifges(7,4) * t309 - Ifges(7,2) * t211) * t303 + (-t111 / 0.2e1 - t124 / 0.2e1) * t23 + (mrSges(7,1) * t262 + mrSges(7,2) * t261) * t50 + t257 * t105 + t258 * t106 + t322 * t69 + (-qJD(5) * t110 - t16 * mrSges(6,3) + t190 * t119 + t99 / 0.2e1 + t51 * mrSges(6,2)) * t188 + t242 * qJD(4) + t320 * mrSges(7,3) + (-t112 / 0.2e1 - t125 / 0.2e1) * t24 + (-Ifges(7,5) * t125 - Ifges(7,6) * t124) * t291 + (-Ifges(7,1) * t125 - Ifges(7,4) * t124) * t299 + (-Ifges(7,4) * t125 - Ifges(7,2) * t124) * t301 + ((-Ifges(5,6) * t250 / 0.2e1 - t199 + (-pkin(3) * mrSges(5,1) + Ifges(6,5) * t288 + Ifges(7,5) * t293 + Ifges(6,6) * t290 + Ifges(7,6) * t294 + t240) * qJD(3) - t177 / 0.2e1) * t195 + ((Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1 + t208) * t252 + ((Ifges(6,1) * t188 - t281) * t289 + (-Ifges(6,2) * t186 + t280) * t288 - qJ(4) * mrSges(5,1) + t239) * qJD(3) + (-Ifges(5,2) / 0.2e1 + Ifges(5,3) / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t250 - t306) * t192) * qJD(2) + (t51 * mrSges(6,1) - t17 * mrSges(6,3) - qJD(5) * t109 + t190 * t120 + t298) * t186 + (-pkin(3) * t69 - qJ(4) * t64 - t101 * t106 + t102 * t308 - t107 * t146) * m(5) + t311 * t56 - t64 * mrSges(5,3) - t68 * mrSges(4,2) - t70 * t27 + t318 * t57 + (t1 * t89 + t173 * t40 + t2 * t88 + t311 * t6 + (qJD(4) - t70) * t50 + t318 * t5) * m(7) + (Ifges(7,5) * t112 + Ifges(7,6) * t111) * t292 + t88 * t28 + t89 * t29 - t96 * t87 - t45 * t109 - t44 * t110 + qJ(4) * t113 - t146 * t144 + t173 * t13 + t12 * t293 + t11 * t294 + (Ifges(7,1) * t112 + Ifges(7,4) * t111) * t300 + (Ifges(7,4) * t112 + Ifges(7,2) * t111) * t302; t188 * t119 + t186 * t120 + t211 * t29 + t309 * t28 + t261 * t57 + t262 * t56 - t242 * qJD(3) + (mrSges(5,1) * t248 + (t109 * t188 - t110 * t186 + t144) * t192) * qJD(2) + (-qJD(3) * t50 - t320) * m(7) + (-qJD(3) * t80 - t215 * t252 + t216) * m(6) + (qJD(3) * t102 + t107 * t252 + t69) * m(5); -t139 * t109 - t138 * t110 - t225 * t56 - t76 * t57 + t272 + (-t225 * t6 - t5 * t76 + t40) * m(7) + (-t138 * t25 - t139 * t26 + t51) * m(6); t170 - t50 * (-mrSges(7,1) * t76 + mrSges(7,2) * t225) + (Ifges(7,1) * t225 + t287) * t300 + t23 * t299 + (Ifges(7,5) * t225 + Ifges(7,6) * t76) * t292 - t5 * t56 + t6 * t57 + (t225 * t5 - t6 * t76) * mrSges(7,3) + (Ifges(7,2) * t76 + t24 + t75) * t302 + t307;];
tauc  = t9(:);
