% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2018-11-23 15:13
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:13:12
% EndTime: 2018-11-23 15:13:18
% DurationCPUTime: 6.43s
% Computational Cost: add. (3208->482), mult. (7947->638), div. (0->0), fcn. (4802->8), ass. (0->227)
t312 = Ifges(6,4) + Ifges(7,4);
t287 = -qJD(2) / 0.2e1;
t313 = Ifges(6,1) + Ifges(7,1);
t300 = Ifges(6,5) + Ifges(7,5);
t311 = Ifges(6,2) + Ifges(7,2);
t299 = Ifges(6,6) + Ifges(7,6);
t156 = sin(qJ(5));
t159 = cos(qJ(5));
t162 = -pkin(3) - pkin(9);
t155 = cos(pkin(6));
t160 = cos(qJ(3));
t244 = t155 * t160;
t137 = qJD(1) * t244;
t157 = sin(qJ(3));
t158 = sin(qJ(2));
t154 = sin(pkin(6));
t234 = qJD(1) * t154;
t208 = t158 * t234;
t121 = qJD(2) * pkin(8) + t208;
t200 = pkin(4) * qJD(2) + t121;
t176 = t200 * t157;
t65 = t137 - t176;
t48 = qJD(3) * t162 + qJD(4) - t65;
t201 = -qJ(4) * t157 - pkin(2);
t107 = t160 * t162 + t201;
t161 = cos(qJ(2));
t207 = t161 * t234;
t67 = qJD(2) * t107 - t207;
t15 = -t156 * t67 + t159 * t48;
t16 = t156 * t48 + t159 * t67;
t178 = t15 * t156 - t159 * t16;
t191 = mrSges(7,1) * t159 - mrSges(7,2) * t156;
t193 = mrSges(6,1) * t159 - mrSges(6,2) * t156;
t231 = qJD(2) * t157;
t143 = qJD(5) + t231;
t227 = qJD(3) * t159;
t229 = qJD(2) * t160;
t112 = -t156 * t229 + t227;
t8 = -qJ(6) * t112 + t15;
t7 = pkin(5) * t143 + t8;
t111 = -qJD(3) * t156 - t159 * t229;
t9 = qJ(6) * t111 + t16;
t194 = t156 * t7 - t159 * t9;
t262 = -t159 / 0.2e1;
t264 = -t156 / 0.2e1;
t308 = t312 * t111;
t282 = t112 * t313 + t143 * t300 + t308;
t307 = t312 * t112;
t283 = t111 * t311 + t143 * t299 + t307;
t153 = qJD(3) * qJ(4);
t233 = qJD(1) * t157;
t136 = t155 * t233;
t86 = t121 * t160 + t136;
t66 = pkin(4) * t229 + t86;
t53 = t153 + t66;
t30 = -pkin(5) * t111 + qJD(6) + t53;
t310 = mrSges(6,3) * t178 + mrSges(7,3) * t194 + t30 * t191 + t53 * t193 + t262 * t283 + t264 * t282;
t309 = -t229 / 0.2e1;
t306 = t312 * t159;
t305 = t312 * t156;
t249 = qJD(2) * pkin(2);
t122 = -t207 - t249;
t265 = t143 / 0.2e1;
t267 = t112 / 0.2e1;
t278 = t156 * t313 + t306;
t279 = t159 * t311 + t305;
t280 = t156 * t300 + t159 * t299;
t285 = qJD(3) / 0.2e1;
t286 = -qJD(3) / 0.2e1;
t288 = t111 / 0.2e1;
t70 = -t153 - t86;
t125 = -pkin(3) * t160 + t201;
t87 = qJD(2) * t125 - t207;
t304 = -t122 * mrSges(4,1) - t70 * mrSges(5,1) + t87 * mrSges(5,2) + t86 * mrSges(4,3) - Ifges(5,5) * t285 - Ifges(4,6) * t286 - t280 * t265 - t278 * t267 - t279 * t288 + t310 + ((-Ifges(4,2) - Ifges(5,3)) * t160 + (-Ifges(4,4) - Ifges(5,6)) * t157) * t287;
t303 = -Ifges(4,1) / 0.2e1;
t271 = pkin(4) + pkin(8);
t302 = Ifges(4,4) * t309;
t284 = mrSges(5,1) + mrSges(4,3);
t301 = mrSges(5,2) - mrSges(4,1);
t220 = qJD(2) * qJD(3);
t202 = t160 * t220;
t219 = qJD(3) * qJD(5);
t222 = qJD(5) * t160;
t228 = qJD(3) * t157;
t77 = -t156 * t219 + (t156 * t228 - t159 * t222) * qJD(2);
t204 = t156 * t222;
t170 = t157 * t227 + t204;
t78 = qJD(2) * t170 - t159 * t219;
t298 = t202 * t299 + t311 * t78 + t312 * t77;
t297 = t202 * t300 + t312 * t78 + t313 * t77;
t224 = qJD(5) * t156;
t147 = pkin(3) * t231;
t179 = pkin(9) * t157 - qJ(4) * t160;
t92 = qJD(2) * t179 + t147;
t23 = -t156 * t92 + t159 * t66;
t240 = qJ(6) - t162;
t296 = -qJD(6) * t159 + t224 * t240 - (-qJ(6) * t156 * t157 + pkin(5) * t160) * qJD(2) - t23;
t124 = t240 * t159;
t24 = t156 * t66 + t159 * t92;
t295 = -qJ(6) * t159 * t231 - qJD(5) * t124 - qJD(6) * t156 - t24;
t209 = -pkin(5) * t159 - pkin(4);
t223 = qJD(5) * t159;
t294 = pkin(5) * t223 + qJD(4) - t137 - (qJD(2) * t209 - t121) * t157;
t85 = t121 * t157 - t137;
t293 = -qJD(4) - t85;
t292 = -t156 * t311 + t306;
t291 = t159 * t313 - t305;
t232 = qJD(2) * t154;
t203 = qJD(1) * t232;
t198 = t161 * t203;
t31 = t157 * t198 + (t160 * t200 + t136) * qJD(3);
t225 = qJD(4) * t157;
t167 = qJD(3) * t179 - t225;
t235 = qJD(3) * t147 + t158 * t203;
t46 = qJD(2) * t167 + t235;
t3 = t156 * t31 + t159 * t46 + t223 * t48 - t224 * t67;
t4 = -qJD(5) * t16 - t156 * t46 + t159 * t31;
t277 = t156 * t3 + t159 * t4;
t211 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t212 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t214 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t69 = -qJD(3) * pkin(3) - t293;
t275 = -t212 * t111 - t214 * t112 - t211 * t143 - t122 * mrSges(4,2) - t15 * mrSges(6,1) - t69 * mrSges(5,1) - t7 * mrSges(7,1) - t85 * mrSges(4,3) - Ifges(5,4) * t286 - (-t157 * Ifges(5,2) - Ifges(5,6) * t160) * t287 + t231 * t303 - Ifges(4,5) * t285 + t302 + t16 * mrSges(6,2) + t87 * mrSges(5,3) + t9 * mrSges(7,2) - t299 * t288 - t300 * t267 - (Ifges(7,3) + Ifges(6,3)) * t265;
t274 = -m(5) / 0.2e1;
t270 = -t111 / 0.2e1;
t268 = -t112 / 0.2e1;
t266 = -t143 / 0.2e1;
t205 = t161 * t232;
t226 = qJD(3) * t160;
t45 = t121 * t226 + (qJD(3) * t155 + t205) * t233;
t245 = t154 * t158;
t210 = t157 * t245;
t96 = t210 - t244;
t258 = t45 * t96;
t49 = mrSges(7,1) * t202 - mrSges(7,3) * t77;
t50 = mrSges(6,1) * t202 - mrSges(6,3) * t77;
t257 = t49 + t50;
t51 = -mrSges(7,2) * t202 + mrSges(7,3) * t78;
t52 = -mrSges(6,2) * t202 + mrSges(6,3) * t78;
t256 = t51 + t52;
t81 = -mrSges(7,2) * t143 + mrSges(7,3) * t111;
t82 = -mrSges(6,2) * t143 + mrSges(6,3) * t111;
t255 = t81 + t82;
t83 = mrSges(7,1) * t143 - mrSges(7,3) * t112;
t84 = mrSges(6,1) * t143 - mrSges(6,3) * t112;
t254 = t83 + t84;
t129 = -mrSges(5,1) * t229 - qJD(3) * mrSges(5,3);
t60 = -mrSges(6,1) * t111 + mrSges(6,2) * t112;
t246 = -t129 + t60;
t131 = t271 * t157;
t57 = t107 * t159 + t131 * t156;
t243 = t156 * t161;
t242 = t159 * t160;
t241 = t159 * t161;
t115 = (mrSges(5,2) * t160 - mrSges(5,3) * t157) * qJD(2);
t239 = t115 + (-mrSges(4,1) * t160 + mrSges(4,2) * t157) * qJD(2);
t238 = qJD(3) * t137 + t160 * t198;
t237 = -qJD(3) * t301 - t231 * t284;
t128 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t229;
t236 = t129 - t128;
t132 = t271 * t160;
t230 = qJD(2) * t158;
t221 = qJD(6) * t160;
t218 = pkin(8) * t45 / 0.2e1;
t59 = -mrSges(7,1) * t111 + mrSges(7,2) * t112;
t217 = -t59 - t246;
t216 = -Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t215 = Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1;
t213 = -0.3e1 / 0.2e1 * Ifges(5,6) - 0.3e1 / 0.2e1 * Ifges(4,4);
t206 = t154 * t230;
t25 = -mrSges(7,1) * t78 + mrSges(7,2) * t77;
t199 = qJ(6) * t160 - t107;
t197 = -t207 / 0.2e1;
t1 = pkin(5) * t202 - qJ(6) * t77 - qJD(6) * t112 + t4;
t2 = qJ(6) * t78 + qJD(6) * t111 + t3;
t196 = -t1 * t159 - t156 * t2;
t192 = mrSges(6,1) * t156 + mrSges(6,2) * t159;
t190 = mrSges(7,1) * t156 + mrSges(7,2) * t159;
t44 = -t121 * t228 + t238;
t174 = t154 * t241 - t156 * t96;
t63 = t154 * t243 + t159 * t96;
t97 = t155 * t157 + t160 * t245;
t172 = -m(4) * t86 + m(5) * t70 + t236;
t171 = -qJ(4) * t226 - t225;
t119 = t271 * t226;
t149 = pkin(3) * t228;
t88 = t149 + t167;
t10 = -t107 * t224 + t119 * t156 + t131 * t223 + t159 * t88;
t169 = -t156 * t254 + t159 * t255;
t168 = t4 * mrSges(6,1) + t1 * mrSges(7,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t152 = qJD(3) * qJD(4);
t27 = -qJD(3) * t176 + t152 + t238;
t163 = qJD(2) ^ 2;
t145 = pkin(5) * t156 + qJ(4);
t141 = Ifges(6,3) * t202;
t140 = Ifges(7,3) * t202;
t123 = t240 * t156;
t118 = t271 * t228;
t117 = -qJ(4) * t229 + t147;
t114 = t159 * t131;
t105 = (mrSges(4,1) * t157 + mrSges(4,2) * t160) * t220;
t104 = (-mrSges(5,2) * t157 - mrSges(5,3) * t160) * t220;
t103 = t159 * t119;
t95 = pkin(5) * t242 + t132;
t94 = t149 + t171;
t80 = (t157 * t243 + t158 * t159) * t234;
t79 = (-t156 * t158 + t157 * t241) * t234;
t75 = Ifges(6,5) * t77;
t74 = Ifges(7,5) * t77;
t73 = Ifges(6,6) * t78;
t72 = Ifges(7,6) * t78;
t68 = -pkin(5) * t204 + (-pkin(8) + t209) * t228;
t62 = qJD(3) * t97 + t157 * t205;
t61 = qJD(3) * t210 - t155 * t226 - t160 * t205;
t58 = qJD(2) * t171 + t235;
t56 = -t107 * t156 + t114;
t43 = -qJ(6) * t242 + t57;
t33 = -t152 - t44;
t32 = pkin(5) * t157 + t156 * t199 + t114;
t26 = -mrSges(6,1) * t78 + mrSges(6,2) * t77;
t14 = -pkin(5) * t78 + t27;
t13 = qJD(5) * t174 - t156 * t206 + t159 * t62;
t12 = qJD(5) * t63 + t156 * t62 + t159 * t206;
t11 = -qJD(5) * t57 - t156 * t88 + t103;
t6 = qJ(6) * t170 - t159 * t221 + t10;
t5 = pkin(5) * t226 + t103 + t199 * t223 + (-qJ(6) * t228 - qJD(5) * t131 + t221 - t88) * t156;
t17 = [(t25 + t26) * t97 - t256 * t174 + t257 * t63 - t237 * t62 + t254 * t13 + t255 * t12 + (-t128 + t217) * t61 + (-t163 * t158 * mrSges(3,1) + (-mrSges(3,2) * t163 - t104 - t105) * t161) * t154 + (t239 * t245 + t284 * qJD(3) * (-t157 * t97 + t160 * t96)) * qJD(2) + m(4) * (t44 * t97 + t258 - t61 * t86 + t62 * t85 + (t122 - t207) * t206) + m(5) * (-t33 * t97 + t258 + t61 * t70 + t62 * t69 + (-t161 * t58 + t230 * t87) * t154) + m(6) * (t12 * t16 + t13 * t15 - t174 * t3 + t27 * t97 + t4 * t63 - t53 * t61) + m(7) * (t1 * t63 + t12 * t9 + t13 * t7 + t14 * t97 - t174 * t2 - t30 * t61); m(6) * (t10 * t16 + t11 * t15 - t118 * t53 + t132 * t27 + t3 * t57 + t4 * t56) + m(5) * (t125 * t58 + t87 * t94) + t125 * t104 + t132 * t26 + t94 * t115 - t118 * t60 + t95 * t25 - pkin(2) * t105 + t6 * t81 + t10 * t82 + t5 * t83 + t11 * t84 + t68 * t59 + t56 * t50 + t57 * t52 + t32 * t49 + t43 * t51 + (t58 * mrSges(5,2) + t27 * t193 + t14 * t191 + t44 * mrSges(4,3) - t33 * mrSges(5,1) + (t1 * t156 - t159 * t2) * mrSges(7,3) + (t156 * t4 - t159 * t3) * mrSges(6,3) + (-mrSges(4,1) * t230 + (-m(6) * t53 - m(7) * t30 + t172 - t59 - t60) * t161) * t234 + (-t53 * t192 - t30 * t190 + (t156 * t9 + t159 * t7) * mrSges(7,3) + (t15 * t159 + t156 * t16) * mrSges(6,3) + t292 * t270 + t291 * t268 + (t156 * t299 - t159 * t300) * t265 + t283 * t156 / 0.2e1) * qJD(5) + (t216 * qJD(3) + ((-t156 * t214 - t159 * t212 - t213) * t160 + (-0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(4,1) + t211) * t157) * qJD(2) - t275) * qJD(3) - t278 * t77 / 0.2e1 - t279 * t78 / 0.2e1 + t297 * t264 + (qJD(5) * t282 + t298) * t262 + (m(4) * t44 + 0.2e1 * t33 * t274 + (m(4) * t85 + m(5) * t69 - t237) * qJD(3)) * pkin(8)) * t160 + m(7) * (t1 * t32 + t14 * t95 + t2 * t43 + t30 * t68 + t5 * t7 + t6 * t9) + 0.2e1 * (t87 * t274 + (-t249 / 0.2e1 - t122 / 0.2e1) * m(4)) * t208 - t254 * t79 - t255 * t80 + (t141 / 0.2e1 + t140 / 0.2e1 + t74 / 0.2e1 + t75 / 0.2e1 + t72 / 0.2e1 + t73 / 0.2e1 - t58 * mrSges(5,3) + t212 * t78 + t214 * t77 + t284 * t45 + (mrSges(4,2) * t230 + t161 * t237) * t234 + 0.2e1 * (t197 * t69 + t218) * m(5) + 0.2e1 * (t197 * t85 + t218) * m(4) + (pkin(8) * t172 + t215 * qJD(3) + t213 * t231 - t304) * qJD(3) + t168) * t157 - t239 * t208 - m(6) * (t15 * t79 + t16 * t80) - m(7) * (t7 * t79 + t80 * t9); t292 * t78 / 0.2e1 + (-pkin(3) * t45 - qJ(4) * t33 - t117 * t87 + t293 * t70 - t69 * t86) * m(5) + t294 * t59 + t145 * t25 + t196 * mrSges(7,3) + t27 * t192 - m(6) * (t15 * t23 + t16 * t24 + t53 * t65) - t123 * t51 - t124 * t49 - t117 * t115 - t24 * t82 - t23 * t84 - t65 * t60 - t33 * mrSges(5,3) - t44 * mrSges(4,2) + qJ(4) * t26 + m(6) * (qJ(4) * t27 + qJD(4) * t53 + t162 * t277) - t277 * mrSges(6,3) + t291 * t77 / 0.2e1 + t295 * t81 + (-t1 * t124 - t123 * t2 + t14 * t145 + t294 * t30 + t295 * t9 + t296 * t7) * m(7) + t296 * t83 + t297 * t159 / 0.2e1 + t298 * t264 + t301 * t45 + (t156 * t52 + t159 * t50) * t162 + t14 * t190 + ((-m(6) * t178 - t156 * t84 + t159 * t82) * t162 + t279 * t270 + t278 * t268 + t280 * t266 + t310) * qJD(5) + ((t302 + (-pkin(3) * mrSges(5,1) - t156 * t212 + t159 * t214 + t216) * qJD(3) + Ifges(5,6) * t309 + t275) * t160 + ((Ifges(4,2) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t303) * t229 + (-mrSges(5,1) * qJ(4) + t215) * qJD(3) + (Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1) * t231 + t304) * t157) * qJD(2) - t236 * t85 + t237 * t86 + t246 * qJD(4); t257 * t159 + t256 * t156 + t217 * qJD(3) + t169 * qJD(5) + (mrSges(5,1) * t226 + (t115 + t169) * t157) * qJD(2) + (-qJD(3) * t30 - t143 * t194 - t196) * m(7) + (-qJD(3) * t53 - t143 * t178 + t277) * m(6) + (qJD(3) * t70 + t231 * t87 + t45) * m(5); t141 + t140 + t74 + t75 + t72 + t73 - t53 * (mrSges(6,1) * t112 + mrSges(6,2) * t111) - t30 * (mrSges(7,1) * t112 + mrSges(7,2) * t111) - t8 * t81 - t15 * t82 + t9 * t83 + t16 * t84 + t168 + (-t112 * t59 + t49) * pkin(5) + (t111 * t7 + t112 * t9) * mrSges(7,3) + (t111 * t15 + t112 * t16) * mrSges(6,3) + (-(-t7 + t8) * t9 + (-t112 * t30 + t1) * pkin(5)) * m(7) + (t111 * t313 - t307) * t268 + t283 * t267 + (t111 * t300 - t112 * t299) * t266 + (-t112 * t311 + t282 + t308) * t270; -t111 * t81 + t112 * t83 + 0.2e1 * (t14 / 0.2e1 + t9 * t270 + t7 * t267) * m(7) + t25;];
tauc  = t17(:);
