% Calculate vector of inverse dynamics joint torques for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:19
% EndTime: 2019-12-31 18:21:34
% DurationCPUTime: 8.85s
% Computational Cost: add. (3396->483), mult. (7286->672), div. (0->0), fcn. (4696->14), ass. (0->220)
t160 = sin(pkin(9));
t162 = cos(pkin(9));
t165 = sin(qJ(5));
t168 = cos(qJ(5));
t188 = t160 * t165 - t162 * t168;
t110 = t188 * qJD(5);
t169 = cos(qJ(3));
t183 = t188 * t169;
t92 = qJD(1) * t183;
t300 = t110 - t92;
t124 = t160 * t168 + t162 * t165;
t111 = t124 * qJD(5);
t184 = t124 * t169;
t91 = qJD(1) * t184;
t299 = t111 - t91;
t229 = qJD(1) * t169;
t141 = qJD(5) - t229;
t166 = sin(qJ(3));
t230 = qJD(1) * t166;
t120 = qJD(3) * t162 - t160 * t230;
t121 = qJD(3) * t160 + t162 * t230;
t203 = t120 * t168 - t121 * t165;
t71 = t120 * t165 + t121 * t168;
t262 = Ifges(6,4) * t71;
t26 = Ifges(6,2) * t203 + Ifges(6,6) * t141 + t262;
t277 = t26 / 0.2e1;
t65 = Ifges(6,4) * t203;
t27 = Ifges(6,1) * t71 + Ifges(6,5) * t141 + t65;
t276 = t27 / 0.2e1;
t268 = m(5) + m(6);
t290 = -m(4) - t268;
t161 = sin(pkin(8));
t145 = pkin(1) * t161 + pkin(6);
t308 = qJD(2) * qJD(3) + t145 * qJDD(1);
t159 = qJ(1) + pkin(8);
t152 = sin(t159);
t154 = cos(t159);
t307 = g(1) * t154 + g(2) * t152;
t224 = qJD(1) * qJD(3);
t131 = qJDD(1) * t166 + t169 * t224;
t93 = qJDD(3) * t162 - t131 * t160;
t94 = qJDD(3) * t160 + t131 * t162;
t23 = qJD(5) * t203 + t165 * t93 + t168 * t94;
t279 = t23 / 0.2e1;
t24 = -qJD(5) * t71 - t165 * t94 + t168 * t93;
t278 = t24 / 0.2e1;
t270 = t93 / 0.2e1;
t269 = t94 / 0.2e1;
t130 = -qJDD(1) * t169 + t166 * t224;
t122 = qJDD(5) + t130;
t267 = t122 / 0.2e1;
t306 = -t130 / 0.2e1;
t266 = t130 / 0.2e1;
t305 = t131 / 0.2e1;
t304 = qJD(3) / 0.2e1;
t231 = t162 * t169;
t187 = pkin(4) * t166 - pkin(7) * t231;
t136 = t145 * qJD(1);
t127 = t166 * t136;
t103 = qJD(2) * t169 - t127;
t192 = pkin(3) * t166 - qJ(4) * t169;
t129 = t192 * qJD(1);
t57 = -t103 * t160 + t129 * t162;
t34 = qJD(1) * t187 + t57;
t217 = t160 * t229;
t58 = t103 * t162 + t129 * t160;
t39 = -pkin(7) * t217 + t58;
t254 = pkin(7) + qJ(4);
t132 = t254 * t160;
t133 = t254 * t162;
t75 = -t132 * t168 - t133 * t165;
t303 = -qJD(4) * t188 + qJD(5) * t75 - t165 * t34 - t168 * t39;
t76 = -t132 * t165 + t133 * t168;
t302 = -qJD(4) * t124 - qJD(5) * t76 + t165 * t39 - t168 * t34;
t301 = Ifges(5,5) * t121 + Ifges(6,5) * t71 + Ifges(5,6) * t120 + Ifges(6,6) * t203 - Ifges(5,3) * t229 + Ifges(6,3) * t141;
t40 = -mrSges(5,1) * t93 + mrSges(5,2) * t94;
t298 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t131 + t40;
t221 = mrSges(4,3) * t230;
t297 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t120 + mrSges(5,2) * t121 + t221;
t296 = t166 * t307;
t199 = -mrSges(5,1) * t162 + mrSges(5,2) * t160;
t182 = m(5) * pkin(3) - t199;
t295 = t169 * t182;
t250 = Ifges(5,4) * t162;
t195 = -Ifges(5,2) * t160 + t250;
t251 = Ifges(5,4) * t160;
t197 = Ifges(5,1) * t162 - t251;
t294 = t120 * (Ifges(5,6) * t166 + t169 * t195) + t121 * (Ifges(5,5) * t166 + t169 * t197);
t219 = t166 * qJDD(2) + t169 * t308;
t227 = qJD(3) * t166;
t54 = -t136 * t227 + t219;
t226 = qJD(3) * t169;
t55 = qJDD(2) * t169 - t136 * t226 - t166 * t308;
t293 = -t166 * t55 + t169 * t54;
t43 = qJDD(3) * qJ(4) + (qJD(4) - t127) * qJD(3) + t219;
t163 = cos(pkin(8));
t261 = pkin(1) * t163;
t147 = -pkin(2) - t261;
t135 = t147 * qJDD(1);
t225 = qJD(4) * t166;
t50 = pkin(3) * t130 - qJ(4) * t131 - qJD(1) * t225 + t135;
t16 = -t160 * t43 + t162 * t50;
t17 = t160 * t50 + t162 * t43;
t191 = -t16 * t160 + t162 * t17;
t89 = mrSges(5,2) * t229 + mrSges(5,3) * t120;
t90 = -mrSges(5,1) * t229 - mrSges(5,3) * t121;
t292 = -t160 * t90 + t162 * t89;
t60 = -mrSges(5,2) * t130 + mrSges(5,3) * t93;
t61 = mrSges(5,1) * t130 - mrSges(5,3) * t94;
t291 = -t160 * t61 + t162 * t60;
t201 = mrSges(4,1) * t169 - mrSges(4,2) * t166;
t289 = t166 * mrSges(6,3) + mrSges(3,1) + t201;
t150 = Ifges(4,4) * t229;
t287 = t162 * (Ifges(5,1) * t121 + Ifges(5,4) * t120 - Ifges(5,5) * t229) + Ifges(4,1) * t230 + Ifges(4,5) * qJD(3) + t150;
t198 = mrSges(5,1) * t160 + mrSges(5,2) * t162;
t200 = mrSges(4,1) * t166 + mrSges(4,2) * t169;
t233 = t160 * t169;
t228 = qJD(2) * t166;
t104 = t136 * t169 + t228;
t88 = qJD(3) * qJ(4) + t104;
t207 = -qJ(4) * t166 - pkin(2);
t117 = -pkin(3) * t169 + t207 - t261;
t95 = t117 * qJD(1);
t35 = -t160 * t88 + t162 * t95;
t36 = t160 * t95 + t162 * t88;
t86 = -qJD(3) * pkin(3) + qJD(4) - t103;
t286 = -t169 * t86 * t198 - t36 * (-mrSges(5,2) * t166 - mrSges(5,3) * t233) - t35 * (mrSges(5,1) * t166 - mrSges(5,3) * t231) - t147 * qJD(1) * t200;
t10 = pkin(4) * t130 - pkin(7) * t94 + t16;
t11 = pkin(7) * t93 + t17;
t29 = -pkin(4) * t229 - pkin(7) * t121 + t35;
t31 = pkin(7) * t120 + t36;
t8 = -t165 * t31 + t168 * t29;
t1 = qJD(5) * t8 + t10 * t165 + t11 * t168;
t9 = t165 * t29 + t168 * t31;
t2 = -qJD(5) * t9 + t10 * t168 - t11 * t165;
t285 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t259 = pkin(4) * t160;
t284 = -m(6) * t259 + mrSges(3,2) - mrSges(4,3) - t198;
t253 = Ifges(4,4) * t166;
t196 = Ifges(4,2) * t169 + t253;
t283 = t8 * mrSges(6,1) - Ifges(4,6) * qJD(3) / 0.2e1 - qJD(1) * t196 / 0.2e1 - t9 * mrSges(6,2);
t281 = Ifges(6,4) * t279 + Ifges(6,2) * t278 + Ifges(6,6) * t267;
t280 = Ifges(6,1) * t279 + Ifges(6,4) * t278 + Ifges(6,5) * t267;
t275 = Ifges(5,1) * t269 + Ifges(5,4) * t270 + Ifges(5,5) * t266;
t274 = -t203 / 0.2e1;
t273 = t203 / 0.2e1;
t272 = -t71 / 0.2e1;
t271 = t71 / 0.2e1;
t265 = -t141 / 0.2e1;
t264 = t141 / 0.2e1;
t167 = sin(qJ(1));
t260 = pkin(1) * t167;
t170 = cos(qJ(1));
t157 = t170 * pkin(1);
t252 = Ifges(4,4) * t169;
t48 = -qJDD(3) * pkin(3) + qJDD(4) - t55;
t242 = t166 * t48;
t108 = qJD(3) * t192 - t225;
t216 = t145 * t227;
t66 = t108 * t162 + t160 * t216;
t236 = t152 * t169;
t235 = t154 * t169;
t234 = t160 * t166;
t232 = t162 * t166;
t73 = t117 * t160 + t145 * t231;
t222 = Ifges(6,5) * t23 + Ifges(6,6) * t24 + Ifges(6,3) * t122;
t220 = mrSges(4,3) * t229;
t214 = m(5) * qJ(4) + mrSges(5,3);
t213 = m(6) * t254 + mrSges(6,3);
t7 = -mrSges(6,1) * t24 + mrSges(6,2) * t23;
t206 = t145 + t259;
t205 = -t224 / 0.2e1;
t204 = t224 / 0.2e1;
t194 = Ifges(4,5) * t169 - Ifges(4,6) * t166;
t193 = Ifges(5,5) * t162 - Ifges(5,6) * t160;
t190 = -t160 * t35 + t162 * t36;
t102 = t162 * t117;
t51 = -pkin(7) * t232 + t102 + (-t145 * t160 - pkin(4)) * t169;
t56 = -pkin(7) * t234 + t73;
t18 = -t165 * t56 + t168 * t51;
t19 = t165 * t51 + t168 * t56;
t146 = pkin(4) * t162 + pkin(3);
t189 = t146 * t169 + t166 * t254;
t185 = t166 * (Ifges(4,1) * t169 - t253);
t158 = pkin(9) + qJ(5);
t151 = sin(t158);
t153 = cos(t158);
t181 = m(6) * t146 + mrSges(6,1) * t153 - mrSges(6,2) * t151;
t175 = t169 * (Ifges(5,3) * t166 + t169 * t193);
t174 = t166 * t214 + t295;
t139 = -qJD(3) * mrSges(4,2) + t220;
t107 = t206 * t166;
t105 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t130;
t100 = t188 * t166;
t99 = t124 * t166;
t98 = t206 * t226;
t96 = t160 * t108;
t84 = t151 * t152 + t153 * t235;
t83 = -t151 * t235 + t152 * t153;
t82 = t151 * t154 - t153 * t236;
t81 = t151 * t236 + t153 * t154;
t77 = t228 + (qJD(1) * t259 + t136) * t169;
t72 = -t145 * t233 + t102;
t67 = -t162 * t216 + t96;
t63 = Ifges(5,4) * t121 + Ifges(5,2) * t120 - Ifges(5,6) * t229;
t59 = -pkin(4) * t120 + t86;
t53 = -qJD(3) * t184 + t110 * t166;
t52 = -qJD(3) * t183 - t111 * t166;
t49 = t96 + (-pkin(7) * t233 - t145 * t232) * qJD(3);
t47 = mrSges(6,1) * t141 - mrSges(6,3) * t71;
t46 = -mrSges(6,2) * t141 + mrSges(6,3) * t203;
t38 = qJD(3) * t187 + t66;
t32 = Ifges(5,4) * t94 + Ifges(5,2) * t93 + Ifges(5,6) * t130;
t30 = -pkin(4) * t93 + t48;
t28 = -mrSges(6,1) * t203 + mrSges(6,2) * t71;
t15 = -mrSges(6,2) * t122 + mrSges(6,3) * t24;
t14 = mrSges(6,1) * t122 - mrSges(6,3) * t23;
t6 = -qJD(5) * t19 - t165 * t49 + t168 * t38;
t5 = qJD(5) * t18 + t165 * t38 + t168 * t49;
t3 = [(Ifges(2,3) + Ifges(3,3) + (0.2e1 * mrSges(3,1) * t163 - 0.2e1 * mrSges(3,2) * t161 + m(3) * (t161 ^ 2 + t163 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (m(4) * t147 - t201) * t135 + (t297 * t226 + t298 * t166 + t169 * t105 + m(5) * (t226 * t86 + t242) + (t293 + (-t103 * t169 - t104 * t166) * qJD(3)) * m(4)) * t145 + (-t104 * mrSges(4,3) + t301 / 0.2e1 + Ifges(6,5) * t271 + Ifges(6,6) * t273 + Ifges(6,3) * t264 + t283) * t227 + ((-Ifges(4,2) * t166 + t252) * t204 - Ifges(6,6) * t278 - Ifges(6,5) * t279 + Ifges(4,4) * t305 + Ifges(4,2) * t306 - Ifges(5,3) * t266 - Ifges(6,3) * t267 - Ifges(5,5) * t269 - Ifges(5,6) * t270 - t16 * mrSges(5,1) + t17 * mrSges(5,2) - t285) * t169 + t293 * mrSges(4,3) - (Ifges(5,5) * t94 + Ifges(5,6) * t93 + Ifges(5,3) * t130 + t222) * t169 / 0.2e1 + t185 * t204 + t175 * t205 + qJDD(3) * (Ifges(4,5) * t166 + Ifges(4,6) * t169) + t147 * (mrSges(4,1) * t130 + mrSges(4,2) * t131) + (Ifges(6,5) * t52 + Ifges(6,6) * t53) * t264 + t107 * t7 + t98 * t28 + t67 * t89 + t66 * t90 + t72 * t61 + t73 * t60 + t59 * (-mrSges(6,1) * t53 + mrSges(6,2) * t52) + t5 * t46 + t6 * t47 + t18 * t14 + t19 * t15 + (-Ifges(6,5) * t100 - Ifges(6,6) * t99) * t267 + t30 * (mrSges(6,1) * t99 - mrSges(6,2) * t100) + (-Ifges(6,1) * t100 - Ifges(6,4) * t99) * t279 + (-Ifges(6,4) * t100 - Ifges(6,2) * t99) * t278 + (-t1 * t99 + t100 * t2 - t52 * t8 + t53 * t9) * mrSges(6,3) + (Ifges(6,1) * t52 + Ifges(6,4) * t53) * t271 + (-t16 * t232 - t17 * t234) * mrSges(5,3) + (m(3) * t260 + mrSges(2,1) * t167 - t82 * mrSges(6,1) + mrSges(2,2) * t170 - t81 * mrSges(6,2) + t290 * (pkin(6) * t154 - t260) + t284 * t154 + (-m(6) * (-pkin(2) - t189) + m(4) * pkin(2) - m(5) * t207 + t166 * mrSges(5,3) + t295 + t289) * t152) * g(1) + (-t160 * t63 / 0.2e1 + t287 / 0.2e1 - t103 * mrSges(4,3)) * t226 + m(6) * (t1 * t19 + t107 * t30 + t18 * t2 + t5 * t9 + t59 * t98 + t6 * t8) + m(5) * (t16 * t72 + t17 * t73 + t35 * t66 + t36 * t67) + t53 * t277 - t100 * t280 - t99 * t281 + (Ifges(6,4) * t52 + Ifges(6,2) * t53) * t273 + (-m(3) * t157 - mrSges(2,1) * t170 - t84 * mrSges(6,1) + mrSges(2,2) * t167 - t83 * mrSges(6,2) + t290 * (pkin(2) * t154 + pkin(6) * t152 + t157) + t284 * t152 + (-m(6) * t189 - t174 - t289) * t154) * g(2) + (Ifges(4,1) * t131 + Ifges(4,4) * t306 + t193 * t266 + t195 * t270 + t197 * t269) * t166 + t232 * t275 + t52 * t276 + t198 * t242 + (t194 * t304 - t286) * qJD(3) + t294 * t304 + t252 * t305 + t196 * t306 - t139 * t216 - t32 * t234 / 0.2e1; m(3) * qJDD(2) - t100 * t15 - t99 * t14 + t52 * t46 + t53 * t47 + (-t7 - t298) * t169 + (t105 + t291) * t166 + ((t139 + t292) * t169 + (t28 + t297) * t166) * qJD(3) + m(6) * (-t1 * t100 - t169 * t30 - t2 * t99 + t227 * t59 + t52 * t9 + t53 * t8) + m(4) * (t166 * t54 + t169 * t55 + (-t103 * t166 + t104 * t169) * qJD(3)) + m(5) * (-t169 * t48 + t191 * t166 + (t166 * t86 + t169 * t190) * qJD(3)) + (-m(3) + t290) * g(3); (-Ifges(6,4) * t92 - Ifges(6,2) * t91) * t274 + t30 * (mrSges(6,1) * t188 + mrSges(6,2) * t124) + (Ifges(6,4) * t124 - Ifges(6,2) * t188) * t278 + (Ifges(6,1) * t124 - Ifges(6,4) * t188) * t279 + (Ifges(6,5) * t124 - Ifges(6,6) * t188) * t267 + (-t1 * t188 - t124 * t2 - t299 * t9 + t300 * t8) * mrSges(6,3) - t188 * t281 + (-t166 * t213 - t169 * t181 - t174 - t201) * g(3) + (-Ifges(6,5) * t92 - Ifges(6,6) * t91) * t265 + (-Ifges(6,1) * t92 - Ifges(6,4) * t91) * t272 - t299 * t277 - t300 * t276 + t194 * t205 + (-t139 + t220) * t103 + t162 * t32 / 0.2e1 - t146 * t7 - Ifges(4,6) * t130 + Ifges(4,5) * t131 - t58 * t89 - t57 * t90 + t75 * t14 + t76 * t15 - t77 * t28 - t54 * mrSges(4,2) + t55 * mrSges(4,1) - pkin(3) * t40 + (-Ifges(6,5) * t110 - Ifges(6,6) * t111) * t264 + (-Ifges(6,1) * t110 - Ifges(6,4) * t111) * t271 + (-Ifges(6,4) * t110 - Ifges(6,2) * t111) * t273 + (mrSges(6,1) * t299 - mrSges(6,2) * t300) * t59 + (-m(5) * t86 + t221 - t297) * t104 - (-Ifges(4,2) * t230 + t150 + t287) * t229 / 0.2e1 + t291 * qJ(4) + t292 * qJD(4) + t191 * mrSges(5,3) - t301 * t230 / 0.2e1 + t302 * t47 + (t1 * t76 - t146 * t30 + t2 * t75 + t302 * t8 + t303 * t9 - t59 * t77) * m(6) + t303 * t46 + t124 * t280 + t307 * (t200 + (-t213 - t214) * t169 + (t181 + t182) * t166) + (Ifges(6,5) * t272 + Ifges(6,6) * t274 + Ifges(6,3) * t265 - t283) * t230 + (Ifges(5,5) * t160 + Ifges(5,6) * t162) * t266 + (Ifges(5,1) * t160 + t250) * t269 + (Ifges(5,2) * t162 + t251) * t270 + t160 * t275 + (-pkin(3) * t48 + qJ(4) * t191 + qJD(4) * t190 - t35 * t57 - t36 * t58) * m(5) + (t286 - t294 / 0.2e1 + (t175 / 0.2e1 - t185 / 0.2e1) * qJD(1)) * qJD(1) + Ifges(4,3) * qJDD(3) + t48 * t199 + t63 * t217 / 0.2e1; t268 * t169 * g(3) - t120 * t89 + t121 * t90 - t203 * t46 + t71 * t47 + t40 + t7 + (-t203 * t9 + t71 * t8 - t296 + t30) * m(6) + (-t120 * t36 + t121 * t35 - t296 + t48) * m(5); -t59 * (mrSges(6,1) * t71 + mrSges(6,2) * t203) + (Ifges(6,1) * t203 - t262) * t272 + t26 * t271 + (Ifges(6,5) * t203 - Ifges(6,6) * t71) * t265 - t8 * t46 + t9 * t47 - g(1) * (mrSges(6,1) * t83 - mrSges(6,2) * t84) - g(2) * (-mrSges(6,1) * t81 + mrSges(6,2) * t82) - g(3) * (-mrSges(6,1) * t151 - mrSges(6,2) * t153) * t166 + (t203 * t8 + t71 * t9) * mrSges(6,3) + t222 + (-Ifges(6,2) * t71 + t27 + t65) * t274 + t285;];
tau = t3;
