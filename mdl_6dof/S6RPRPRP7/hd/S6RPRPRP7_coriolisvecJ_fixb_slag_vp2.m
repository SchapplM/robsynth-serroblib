% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2018-11-23 16:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:00:41
% EndTime: 2018-11-23 16:00:50
% DurationCPUTime: 9.70s
% Computational Cost: add. (5028->471), mult. (11425->640), div. (0->0), fcn. (7410->6), ass. (0->215)
t318 = Ifges(6,4) + Ifges(7,4);
t319 = Ifges(6,1) + Ifges(7,1);
t317 = Ifges(6,5) + Ifges(7,5);
t316 = Ifges(6,2) + Ifges(7,2);
t315 = Ifges(6,6) + Ifges(7,6);
t314 = Ifges(6,3) + Ifges(7,3);
t154 = sin(pkin(9));
t246 = cos(pkin(9));
t274 = cos(qJ(3));
t187 = t246 * t274;
t173 = qJD(1) * t187;
t156 = sin(qJ(3));
t234 = qJD(1) * t156;
t118 = t154 * t234 - t173;
t155 = sin(qJ(5));
t157 = cos(qJ(5));
t95 = qJD(3) * t157 + t118 * t155;
t329 = t318 * t95;
t96 = qJD(3) * t155 - t118 * t157;
t328 = t318 * t96;
t162 = -t154 * t274 - t156 * t246;
t119 = t162 * qJD(1);
t324 = qJD(5) - t119;
t308 = t315 * t324 + t316 * t95 + t328;
t307 = t317 * t324 + t319 * t96 + t329;
t231 = qJD(5) * t157;
t232 = qJD(5) * t155;
t204 = t274 * qJD(4);
t158 = -pkin(1) - pkin(7);
t139 = qJD(1) * t158 + qJD(2);
t233 = qJD(3) * t156;
t216 = t139 * t233;
t160 = -t216 + (qJ(4) * t233 - t204) * qJD(1);
t205 = qJD(3) * t274;
t192 = t139 * t205;
t230 = t156 * qJD(4);
t92 = t192 + (-qJ(4) * t205 - t230) * qJD(1);
t44 = t154 * t160 + t246 * t92;
t229 = qJD(1) * qJD(3);
t203 = t156 * t229;
t110 = -qJD(3) * t173 + t154 * t203;
t111 = qJD(3) * t119;
t152 = qJD(1) * qJD(2);
t189 = qJD(1) * t205;
t131 = pkin(3) * t189 + t152;
t49 = -pkin(4) * t110 - pkin(8) * t111 + t131;
t130 = t274 * t139;
t206 = qJD(1) * t274;
t115 = -qJ(4) * t206 + t130;
t106 = qJD(3) * pkin(3) + t115;
t114 = -qJ(4) * t234 + t139 * t156;
t197 = t246 * t114;
t65 = t154 * t106 + t197;
t57 = qJD(3) * pkin(8) + t65;
t132 = pkin(3) * t234 + qJD(1) * qJ(2) + qJD(4);
t70 = -pkin(4) * t119 + pkin(8) * t118 + t132;
t5 = t155 * t49 + t157 * t44 + t70 * t231 - t232 * t57;
t23 = t155 * t70 + t157 * t57;
t6 = -qJD(5) * t23 - t155 * t44 + t157 * t49;
t186 = -t6 * t155 + t5 * t157;
t22 = -t155 * t57 + t157 * t70;
t327 = -t22 * t231 - t23 * t232 + t186;
t326 = t318 * t157;
t325 = t318 * t155;
t323 = qJ(2) * (m(3) + m(4)) + mrSges(3,3);
t120 = -qJD(3) * t187 + t154 * t233;
t321 = t120 / 0.2e1;
t121 = t162 * qJD(3);
t320 = t121 / 0.2e1;
t273 = pkin(3) * t154;
t145 = pkin(8) + t273;
t235 = qJ(6) + t145;
t195 = qJD(5) * t235;
t103 = t154 * t114;
t78 = t115 * t246 - t103;
t193 = pkin(3) * t206;
t80 = -t118 * pkin(4) - t119 * pkin(8) + t193;
t24 = -t155 * t78 + t157 * t80;
t244 = t119 * t157;
t313 = pkin(5) * t118 + qJ(6) * t244 - qJD(6) * t155 - t157 * t195 - t24;
t58 = qJD(5) * t95 + t111 * t157;
t59 = -qJD(5) * t96 - t111 * t155;
t312 = -t110 * t315 + t316 * t59 + t318 * t58;
t311 = -t317 * t110 + t318 * t59 + t319 * t58;
t245 = t119 * t155;
t25 = t155 * t80 + t157 * t78;
t310 = qJ(6) * t245 + qJD(6) * t157 - t155 * t195 - t25;
t309 = t314 * t324 + t315 * t95 + t317 * t96;
t77 = t115 * t154 + t197;
t306 = -t77 + (t232 - t245) * pkin(5);
t21 = -mrSges(6,1) * t59 + mrSges(6,2) * t58;
t254 = t111 * mrSges(5,3);
t305 = t254 + t21;
t253 = t118 * mrSges(5,3);
t304 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t95 + mrSges(6,2) * t96 - t253;
t303 = Ifges(5,5) * qJD(3);
t302 = Ifges(5,6) * qJD(3);
t236 = qJ(4) - t158;
t183 = mrSges(7,1) * t155 + mrSges(7,2) * t157;
t184 = mrSges(6,1) * t155 + mrSges(6,2) * t157;
t64 = t106 * t246 - t103;
t56 = -qJD(3) * pkin(4) - t64;
t28 = -t95 * pkin(5) + qJD(6) + t56;
t301 = t28 * t183 + t56 * t184;
t300 = -t155 * t315 + t157 * t317;
t299 = -t155 * t316 + t326;
t298 = t157 * t319 - t325;
t296 = -t110 * t314 + t315 * t59 + t317 * t58;
t127 = t154 * t156 - t187;
t43 = t154 * t92 - t246 * t160;
t247 = t127 * t43;
t294 = -t120 * t65 + t121 * t64 + t247;
t1 = -pkin(5) * t110 - qJ(6) * t58 - qJD(6) * t96 + t6;
t12 = qJ(6) * t95 + t23;
t2 = qJ(6) * t59 + qJD(6) * t95 + t5;
t11 = -qJ(6) * t96 + t22;
t9 = pkin(5) * t324 + t11;
t293 = -t1 * t155 - t12 * t232 + t2 * t157 - t9 * t231;
t292 = -m(6) * t56 - t304;
t291 = qJD(1) ^ 2;
t290 = t58 / 0.2e1;
t289 = t59 / 0.2e1;
t288 = -t95 / 0.2e1;
t286 = -t96 / 0.2e1;
t285 = t96 / 0.2e1;
t284 = -t110 / 0.2e1;
t283 = -t324 / 0.2e1;
t281 = -t118 / 0.2e1;
t280 = t118 / 0.2e1;
t279 = -t119 / 0.2e1;
t133 = t236 * t156;
t134 = t236 * t274;
t90 = -t133 * t154 + t246 * t134;
t270 = t43 * t90;
t37 = -mrSges(7,1) * t110 - mrSges(7,3) * t58;
t38 = -mrSges(6,1) * t110 - mrSges(6,3) * t58;
t264 = -t37 - t38;
t39 = mrSges(7,2) * t110 + mrSges(7,3) * t59;
t40 = mrSges(6,2) * t110 + mrSges(6,3) * t59;
t263 = t39 + t40;
t66 = -mrSges(7,2) * t324 + t95 * mrSges(7,3);
t67 = -mrSges(6,2) * t324 + mrSges(6,3) * t95;
t262 = t66 + t67;
t68 = mrSges(7,1) * t324 - mrSges(7,3) * t96;
t69 = mrSges(6,1) * t324 - mrSges(6,3) * t96;
t261 = t68 + t69;
t147 = t156 * pkin(3) + qJ(2);
t87 = -pkin(4) * t162 + pkin(8) * t127 + t147;
t91 = -t133 * t246 - t154 * t134;
t88 = t157 * t91;
t30 = t155 * t87 + t88;
t260 = Ifges(4,4) * t156;
t255 = t110 * mrSges(5,3);
t252 = t118 * Ifges(5,4);
t251 = t119 * mrSges(5,3);
t18 = -pkin(5) * t59 + t43;
t248 = t127 * t18;
t243 = t120 * t155;
t242 = t120 * t157;
t241 = t121 * t155;
t240 = t121 * t157;
t239 = t127 * t155;
t238 = t127 * t157;
t137 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t206;
t237 = t156 * t137;
t140 = pkin(3) * t205 + qJD(2);
t112 = t233 * t236 - t204;
t113 = -qJD(3) * t134 - t230;
t74 = t154 * t112 + t113 * t246;
t79 = -pkin(4) * t120 - pkin(8) * t121 + t140;
t228 = t155 * t79 + t157 * t74 + t87 * t231;
t47 = -mrSges(7,1) * t95 + mrSges(7,2) * t96;
t227 = t47 + t304;
t221 = Ifges(4,4) * t274;
t217 = t246 * pkin(3);
t215 = t127 * t231;
t20 = -t59 * mrSges(7,1) + t58 * mrSges(7,2);
t199 = t231 / 0.2e1;
t198 = -t155 * t74 + t157 * t79;
t29 = -t155 * t91 + t157 * t87;
t196 = -t110 * mrSges(5,1) + t111 * mrSges(5,2);
t73 = -t246 * t112 + t113 * t154;
t146 = -t217 - pkin(4);
t176 = t23 * t155 + t22 * t157;
t174 = -qJ(6) * t121 + qJD(6) * t127;
t172 = -Ifges(4,5) * t156 - Ifges(4,6) * t274;
t169 = t215 - t241;
t168 = t127 * t232 + t240;
t167 = qJ(2) * (mrSges(4,1) * t274 - mrSges(4,2) * t156);
t166 = t156 * (-Ifges(4,2) * t274 - t260);
t165 = (Ifges(4,1) * t274 - t260) * qJD(1);
t164 = (-Ifges(4,2) * t156 + t221) * qJD(1);
t129 = (t156 * mrSges(4,1) + mrSges(4,2) * t274) * qJD(1);
t163 = -t155 * t262 - t157 * t261;
t161 = (-Ifges(4,1) * t156 - t221) * t274;
t136 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t234;
t135 = -t157 * pkin(5) + t146;
t125 = t235 * t157;
t124 = t235 * t155;
t123 = Ifges(4,5) * qJD(3) + t165;
t122 = Ifges(4,6) * qJD(3) + t164;
t116 = Ifges(5,4) * t119;
t101 = -qJD(3) * mrSges(5,2) + t251;
t86 = -mrSges(5,1) * t119 - mrSges(5,2) * t118;
t83 = -t118 * Ifges(5,1) + t116 + t303;
t82 = t119 * Ifges(5,2) - t252 + t302;
t62 = -pkin(5) * t239 + t90;
t27 = -pkin(5) * t169 + t73;
t26 = qJ(6) * t239 + t30;
t19 = -pkin(5) * t162 + qJ(6) * t238 + t29;
t8 = -qJD(5) * t30 + t198;
t7 = -t232 * t91 + t228;
t4 = qJ(6) * t215 + (-qJD(5) * t91 + t174) * t155 + t228;
t3 = -pkin(5) * t120 + t174 * t157 + (-t88 + (-qJ(6) * t127 - t87) * t155) * qJD(5) + t198;
t10 = [-(t164 + t122) * t205 / 0.2e1 + m(5) * (t131 * t147 + t132 * t140 + t44 * t91 + t65 * t74 + t270) + m(6) * (t22 * t8 + t23 * t7 + t29 * t6 + t30 * t5 + t270) + (qJD(5) * t307 + t312) * t239 / 0.2e1 + (-t127 * t110 + t111 * t162 + t119 * t320 + t120 * t281) * Ifges(5,4) + (t110 * t162 + t119 * t321) * Ifges(5,2) + (-t127 * t300 - t162 * t314) * t284 + (-t127 * t299 - t162 * t315) * t289 + (-t127 * t298 - t162 * t317) * t290 - t296 * t162 / 0.2e1 + (t162 * t44 - t294) * mrSges(5,3) + t1 * (-mrSges(7,1) * t162 + mrSges(7,3) * t238) + t6 * (-mrSges(6,1) * t162 + mrSges(6,3) * t238) + t2 * (mrSges(7,2) * t162 + mrSges(7,3) * t239) + t5 * (mrSges(6,2) * t162 + mrSges(6,3) * t239) + t131 * (-mrSges(5,1) * t162 - mrSges(5,2) * t127) + (0.2e1 * t167 - t166 + t161) * t229 + (t136 * t205 - t137 * t233) * t158 + (t127 * t199 - t241 / 0.2e1) * t308 + 0.2e1 * t323 * t152 - t311 * t238 / 0.2e1 + (-t120 * t315 + t168 * t318 + t169 * t316) * t95 / 0.2e1 + (-t120 * t317 + t168 * t319 + t169 * t318) * t285 + t147 * t196 + (-m(5) * t64 - t292) * t73 + t305 * t90 - (t123 + t165) * t233 / 0.2e1 + (-t111 * t127 + t121 * t281) * Ifges(5,1) + qJD(3) * (Ifges(5,5) * t121 + Ifges(5,6) * t120) / 0.2e1 + t132 * (-mrSges(5,1) * t120 + mrSges(5,2) * t121) + (-t120 * t314 + t168 * t317 + t169 * t315) * t324 / 0.2e1 + m(7) * (t1 * t19 + t12 * t4 + t18 * t62 + t2 * t26 + t27 * t28 + t3 * t9) + t83 * t320 + t82 * t321 + t28 * (-mrSges(7,1) * t169 + mrSges(7,2) * t168) + t9 * (-mrSges(7,1) * t120 - mrSges(7,3) * t168) + t22 * (-mrSges(6,1) * t120 - mrSges(6,3) * t168) + t23 * (mrSges(6,2) * t120 + mrSges(6,3) * t169) + t12 * (mrSges(7,2) * t120 + mrSges(7,3) * t169) + t56 * (-mrSges(6,1) * t169 + mrSges(6,2) * t168) + t307 * t240 / 0.2e1 - t309 * t120 / 0.2e1 + qJD(3) ^ 2 * t172 / 0.2e1 - t184 * t247 - t183 * t248 + t19 * t37 + t29 * t38 + t26 * t39 + t30 * t40 + t27 * t47 + t91 * t255 + t62 * t20 + t4 * t66 + t7 * t67 + t3 * t68 + t8 * t69 + t74 * t101 + 0.2e1 * qJD(2) * t129 + t140 * t86; (t136 * t274 - t237) * qJD(3) - t323 * t291 + (t20 + t305) * t127 - t227 * t121 + (t155 * t261 - t157 * t262 - t101) * t120 + m(5) * t294 + m(6) * (-t121 * t56 + t22 * t243 - t23 * t242 + t247) + m(7) * (-t12 * t242 - t121 * t28 + t243 * t9 + t248) - (m(5) * t44 + m(6) * t327 + m(7) * t293 + t163 * qJD(5) + t264 * t155 + t263 * t157 + t255) * t162 + (-m(5) * t132 - t129 - t86 - m(6) * t176 - m(7) * (t12 * t155 + t9 * t157) + t163) * qJD(1); (m(6) * (-qJD(5) * t176 + t186) - t69 * t231 - t67 * t232 - t155 * t38 + t157 * t40) * t145 + (-t232 / 0.2e1 + t245 / 0.2e1) * t308 + (t199 - t244 / 0.2e1) * t307 + (t157 * t316 + t325) * t289 + (t155 * t319 + t326) * t290 + (t22 * t244 + t23 * t245 + t327) * mrSges(6,3) + t311 * t155 / 0.2e1 + t312 * t157 / 0.2e1 + (-t1 * t124 + t12 * t310 + t125 * t2 + t135 * t18 + t28 * t306 + t313 * t9) * m(7) + t313 * t68 + (t155 * t317 + t157 * t315) * t284 + (t9 * mrSges(7,1) + t22 * mrSges(6,1) - t12 * mrSges(7,2) - t23 * mrSges(6,2) + Ifges(5,2) * t279 - t302 / 0.2e1 + t132 * mrSges(5,1) - t315 * t288 - t317 * t286 - t314 * t283) * t118 + (t12 * t245 + t244 * t9 + t293) * mrSges(7,3) + t292 * t77 - t86 * t193 - Ifges(4,6) * t189 + (-t132 * t193 + t64 * t77 - t65 * t78 + (t154 * t44 - t246 * t43) * pkin(3)) * m(5) + (Ifges(5,1) * t280 - t303 / 0.2e1 - t132 * mrSges(5,2) + t299 * t288 + t298 * t286 + t300 * t283 - t301) * t119 + t306 * t47 + (m(6) * t146 - mrSges(6,1) * t157 + mrSges(6,2) * t155 - mrSges(5,1)) * t43 + (t116 + t83) * t279 + t301 * qJD(5) - t65 * t253 - t217 * t254 + (t298 * t96 + t299 * t95 + t300 * t324) * qJD(5) / 0.2e1 + (-t161 / 0.2e1 + t166 / 0.2e1 - t167) * t291 + t123 * t234 / 0.2e1 - t172 * t229 / 0.2e1 - Ifges(4,5) * t203 - t136 * t130 - mrSges(4,1) * t216 + t122 * t206 / 0.2e1 + (t252 + t309) * t280 + t310 * t66 - mrSges(4,2) * t192 - m(6) * (t22 * t24 + t23 * t25) - t44 * mrSges(5,2) + t139 * t237 + t64 * t251 + t255 * t273 + t82 * t281 - t25 * t67 - t24 * t69 - t78 * t101 + Ifges(5,6) * t110 + Ifges(5,5) * t111 - t124 * t37 + t125 * t39 + t135 * t20 + t146 * t21 + t18 * (-mrSges(7,1) * t157 + mrSges(7,2) * t155); -t119 * t101 + t227 * t118 + (t262 * t324 - t264) * t157 + (-t261 * t324 + t263) * t155 + t196 + (t1 * t157 + t118 * t28 + t155 * t2 + t324 * (t12 * t157 - t155 * t9)) * m(7) + (t118 * t56 + t155 * t5 + t157 * t6 + t324 * (-t155 * t22 + t157 * t23)) * m(6) + (-t118 * t64 - t119 * t65 + t131) * m(5); t296 + (-(t11 - t9) * t12 + (-t28 * t96 + t1) * pkin(5)) * m(7) + (-t47 * t96 + t37) * pkin(5) - t5 * mrSges(6,2) + t6 * mrSges(6,1) - t11 * t66 - t22 * t67 + t12 * t68 + t23 * t69 + t1 * mrSges(7,1) - t2 * mrSges(7,2) - t28 * (mrSges(7,1) * t96 + mrSges(7,2) * t95) - t56 * (mrSges(6,1) * t96 + mrSges(6,2) * t95) + (t22 * t95 + t23 * t96) * mrSges(6,3) + (t12 * t96 + t9 * t95) * mrSges(7,3) + (t319 * t95 - t328) * t286 + t308 * t285 + (-t315 * t96 + t317 * t95) * t283 + (-t316 * t96 + t307 + t329) * t288; -t95 * t66 + t96 * t68 + 0.2e1 * (t18 / 0.2e1 + t12 * t288 + t9 * t285) * m(7) + t20;];
tauc  = t10(:);
