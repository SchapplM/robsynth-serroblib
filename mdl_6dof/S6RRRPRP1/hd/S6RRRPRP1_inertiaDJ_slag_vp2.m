% Calculate time derivative of joint inertia matrix for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:30:46
% EndTime: 2019-03-09 16:30:57
% DurationCPUTime: 5.76s
% Computational Cost: add. (7243->402), mult. (16043->574), div. (0->0), fcn. (15392->8), ass. (0->181)
t304 = Ifges(6,4) + Ifges(7,4);
t303 = Ifges(6,1) + Ifges(7,1);
t302 = Ifges(6,2) + Ifges(7,2);
t257 = Ifges(6,5) + Ifges(7,5);
t169 = sin(qJ(3));
t170 = sin(qJ(2));
t172 = cos(qJ(3));
t173 = cos(qJ(2));
t227 = t172 * t173;
t131 = -t169 * t170 + t227;
t132 = t169 * t173 + t172 * t170;
t167 = sin(pkin(10));
t239 = cos(pkin(10));
t98 = -t131 * t239 + t132 * t167;
t301 = t257 * t98;
t290 = Ifges(6,6) + Ifges(7,6);
t300 = t290 * t98;
t171 = cos(qJ(5));
t299 = t304 * t171;
t168 = sin(qJ(5));
t298 = t304 * t168;
t296 = -t302 * t168 + t299;
t295 = t303 * t171 - t298;
t294 = t257 * t168 + t290 * t171;
t292 = t171 * t302 + t298;
t291 = t168 * t303 + t299;
t289 = Ifges(6,3) + Ifges(7,3);
t219 = qJD(5) * t168;
t279 = qJD(2) + qJD(3);
t102 = t279 * t131;
t103 = t279 * t132;
t77 = t102 * t239 - t167 * t103;
t243 = t171 * t77;
t99 = t167 * t131 + t132 * t239;
t178 = t99 * t219 - t243;
t218 = qJD(5) * t171;
t208 = t99 * t218;
t245 = t168 * t77;
t179 = t208 + t245;
t76 = t102 * t167 + t103 * t239;
t288 = -t178 * t304 - t179 * t302 + t290 * t76;
t287 = -t178 * t303 - t179 * t304 + t257 * t76;
t286 = t296 * t99 + t300;
t254 = t295 * t99 + t301;
t281 = -mrSges(6,1) * t171 + mrSges(6,2) * t168 - mrSges(5,1);
t285 = t296 * qJD(5);
t284 = t295 * qJD(5);
t224 = t257 * t218;
t230 = t167 * t169;
t247 = pkin(2) * qJD(3);
t116 = (t172 * t239 - t230) * t247;
t194 = (t168 ^ 2 + t171 ^ 2) * t116;
t267 = -pkin(8) - pkin(7);
t150 = t267 * t170;
t151 = t267 * t173;
t106 = t172 * t150 + t151 * t169;
t177 = -qJ(4) * t132 + t106;
t133 = t169 * t150;
t107 = -t172 * t151 + t133;
t94 = qJ(4) * t131 + t107;
t54 = t167 * t177 + t239 * t94;
t51 = t171 * t54;
t157 = -pkin(2) * t173 - pkin(1);
t112 = -t131 * pkin(3) + t157;
t52 = t98 * pkin(4) - t99 * pkin(9) + t112;
t27 = t168 * t52 + t51;
t280 = t27 * qJD(5);
t175 = t284 * t168 + t285 * t171 + t218 * t291 - t219 * t292;
t278 = 2 * m(5);
t277 = 2 * m(6);
t276 = 2 * m(7);
t181 = -t102 * qJ(4) - t132 * qJD(4);
t221 = qJD(3) * t172;
t222 = qJD(3) * t169;
t79 = qJD(2) * t132 * t267 + t150 * t221 + t151 * t222;
t46 = -qJ(4) * t103 + qJD(4) * t131 + t79;
t176 = (t227 * t267 - t133) * qJD(2);
t80 = -qJD(3) * t107 + t176;
t24 = t167 * t46 - t239 * (t181 + t80);
t275 = 0.2e1 * t24;
t53 = t167 * t94 - t239 * t177;
t274 = 0.2e1 * t53;
t96 = qJD(2) * t170 * pkin(2) + pkin(3) * t103;
t273 = 0.2e1 * t96;
t272 = 0.2e1 * t157;
t271 = m(5) * pkin(3);
t270 = m(7) * pkin(5);
t266 = mrSges(7,3) * pkin(5);
t263 = pkin(3) * t167;
t262 = t171 * pkin(5);
t25 = t239 * t46 + (-t150 * t222 + t151 * t221 + t176 + t181) * t167;
t34 = pkin(4) * t76 - pkin(9) * t77 + t96;
t216 = t168 * t34 + t171 * t25 + t52 * t218;
t5 = -t219 * t54 + t216;
t261 = t171 * t5;
t260 = t24 * t53;
t197 = -t168 * t25 + t171 * t34;
t6 = t197 - t280;
t259 = t6 * t168;
t258 = t76 * mrSges(5,3);
t253 = mrSges(6,2) * t171;
t196 = t239 * t169;
t115 = (t167 * t172 + t196) * t247;
t246 = t115 * t53;
t244 = t168 * t99;
t242 = t171 * t99;
t163 = t171 * qJD(6);
t156 = pkin(2) * t172 + pkin(3);
t118 = pkin(2) * t196 + t167 * t156;
t114 = pkin(9) + t118;
t226 = -qJ(6) - t114;
t193 = qJD(5) * t226;
t81 = t116 * t171 + t168 * t193 + t163;
t241 = t81 * t171;
t82 = (-qJD(6) - t116) * t168 + t171 * t193;
t240 = t82 * t168;
t214 = pkin(5) * t219;
t108 = t115 + t214;
t143 = -mrSges(7,1) * t171 + mrSges(7,2) * t168;
t238 = t108 * t143;
t117 = -pkin(2) * t230 + t156 * t239;
t113 = -pkin(4) - t117;
t109 = t113 - t262;
t136 = mrSges(7,1) * t219 + mrSges(7,2) * t218;
t237 = t109 * t136;
t154 = pkin(9) + t263;
t225 = -qJ(6) - t154;
t192 = qJD(5) * t225;
t110 = t168 * t192 + t163;
t236 = t110 * t171;
t111 = -qJD(6) * t168 + t171 * t192;
t235 = t111 * t168;
t187 = mrSges(6,1) * t168 + t253;
t137 = t187 * qJD(5);
t234 = t113 * t137;
t202 = t239 * pkin(3);
t155 = -t202 - pkin(4);
t142 = t155 - t262;
t233 = t142 * t136;
t232 = t154 * t171;
t231 = t155 * t137;
t164 = t171 * qJ(6);
t220 = qJD(5) * t114;
t217 = 0.2e1 * t173;
t215 = 0.2e1 * qJD(5);
t26 = -t168 * t54 + t171 * t52;
t213 = t26 * qJD(5) * mrSges(6,3);
t212 = mrSges(7,1) + t270;
t211 = mrSges(7,3) * t219;
t209 = mrSges(7,3) * t218;
t201 = t76 * mrSges(5,1) + t77 * mrSges(5,2);
t200 = t290 * t168;
t199 = -t219 / 0.2e1;
t195 = t243 * t257 + t289 * t76;
t191 = t143 * t214;
t188 = -(2 * Ifges(5,4)) - t200;
t182 = -qJ(6) * t77 - qJD(6) * t99;
t19 = mrSges(7,1) * t179 - mrSges(7,2) * t178;
t3 = -qJ(6) * t208 + (-qJD(5) * t54 + t182) * t168 + t216;
t36 = pkin(5) * t244 + t53;
t8 = pkin(5) * t179 + t24;
t174 = (-t219 * t290 + t224) * t98 / 0.2e1 + t254 * t218 / 0.2e1 + t286 * t199 + t287 * t168 / 0.2e1 + t281 * t24 + t284 * t242 / 0.2e1 - t285 * t244 / 0.2e1 + mrSges(6,3) * t261 - t25 * mrSges(5,2) + Ifges(5,5) * t77 - t79 * mrSges(4,2) + t80 * mrSges(4,1) + Ifges(4,5) * t102 - Ifges(4,6) * t103 + t36 * t136 + t53 * t137 + t8 * t143 + t292 * (-t208 / 0.2e1 - t245 / 0.2e1) + t291 * (t99 * t199 + t243 / 0.2e1) + (t294 / 0.2e1 - Ifges(5,6)) * t76 + (t288 / 0.2e1 + t3 * mrSges(7,3)) * t171;
t127 = t164 + t232;
t126 = t225 * t168;
t105 = t114 * t171 + t164;
t104 = t226 * t168;
t66 = mrSges(6,1) * t98 - mrSges(6,3) * t242;
t65 = mrSges(7,1) * t98 - mrSges(7,3) * t242;
t64 = -mrSges(6,2) * t98 - mrSges(6,3) * t244;
t63 = -mrSges(7,2) * t98 - mrSges(7,3) * t244;
t62 = t187 * t99;
t61 = (mrSges(7,1) * t168 + mrSges(7,2) * t171) * t99;
t31 = -mrSges(6,2) * t76 - mrSges(6,3) * t179;
t30 = -mrSges(7,2) * t76 - mrSges(7,3) * t179;
t29 = mrSges(6,1) * t76 + mrSges(6,3) * t178;
t28 = mrSges(7,1) * t76 + mrSges(7,3) * t178;
t21 = -qJ(6) * t244 + t27;
t20 = mrSges(6,1) * t179 - mrSges(6,2) * t178;
t17 = pkin(5) * t98 - t164 * t99 + t26;
t1 = pkin(5) * t76 + t182 * t171 + (-t51 + (qJ(6) * t99 - t52) * t168) * qJD(5) + t197;
t2 = [((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t173) * t217 + (m(4) * pkin(2) * t272 + 0.2e1 * pkin(2) * (-mrSges(4,1) * t131 + mrSges(4,2) * t132) - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t170 + (Ifges(3,1) - Ifges(3,2)) * t217) * t170) * qJD(2) - 0.2e1 * t54 * t258 - 0.2e1 * t131 * Ifges(4,2) * t103 + 0.2e1 * t132 * t102 * Ifges(4,1) + 0.2e1 * t112 * t201 + (mrSges(5,1) * t273 - 0.2e1 * mrSges(5,3) * t25 + t188 * t77 + ((2 * Ifges(5,2)) + t289) * t76 + t195) * t98 + (mrSges(5,3) * t274 - t168 * t286 + t171 * t254) * t77 + (t1 * t17 + t21 * t3 + t36 * t8) * t276 + (t26 * t6 + t27 * t5 + t260) * t277 + (t112 * t96 + t25 * t54 + t260) * t278 + t20 * t274 + t62 * t275 + (mrSges(5,2) * t273 + mrSges(5,3) * t275 + 0.2e1 * Ifges(5,1) * t77 + t287 * t171 - t288 * t168 + (t171 * t257 + t188) * t76 + ((-t286 - t300) * t171 + (-t254 - t301) * t168) * qJD(5)) * t99 + 0.2e1 * m(4) * (t106 * t80 + t107 * t79) + (mrSges(4,1) * t103 + mrSges(4,2) * t102) * t272 + 0.2e1 * (t102 * t131 - t103 * t132) * Ifges(4,4) + 0.2e1 * (-t102 * t106 - t103 * t107 + t131 * t79 - t132 * t80) * mrSges(4,3) + 0.2e1 * t17 * t28 + 0.2e1 * t26 * t29 + 0.2e1 * t21 * t30 + 0.2e1 * t27 * t31 + 0.2e1 * t36 * t19 + 0.2e1 * t8 * t61 + 0.2e1 * t3 * t63 + 0.2e1 * t5 * t64 + 0.2e1 * t1 * t65 + 0.2e1 * t6 * t66; t174 + (t115 * t99 - t116 * t98 - t117 * t77 - t118 * t76) * mrSges(5,3) + (m(4) * (t169 * t79 + t172 * t80 + (-t106 * t169 + t107 * t172) * qJD(3)) + (-t172 * t102 - t169 * t103 + (t131 * t172 + t132 * t169) * qJD(3)) * mrSges(4,3)) * pkin(2) + (m(6) * (-t114 * t6 - t116 * t26 - t220 * t27) - t114 * t29 - t116 * t66 - t64 * t220 + (-t21 * qJD(5) - t1) * mrSges(7,3) + (-t6 - t280) * mrSges(6,3)) * t168 + (m(6) * (t114 * t5 + t116 * t27 - t220 * t26) + t114 * t31 + t116 * t64 - t213 - t17 * qJD(5) * mrSges(7,3) - t66 * t220) * t171 + m(5) * (t116 * t54 - t117 * t24 + t118 * t25 + t246) + m(7) * (t1 * t104 + t105 * t3 + t108 * t36 + t109 * t8 + t17 * t82 + t21 * t81) + (Ifges(3,5) * t173 - Ifges(3,6) * t170 + (-mrSges(3,1) * t173 + mrSges(3,2) * t170) * pkin(7)) * qJD(2) + m(6) * (t113 * t24 + t246) + t81 * t63 + t82 * t65 + t104 * t28 + t105 * t30 + t108 * t61 + t109 * t19 + t113 * t20 + t115 * t62; 0.2e1 * t238 + 0.2e1 * t237 + 0.2e1 * t234 - 0.2e1 * t116 * mrSges(5,2) + (t113 * t115 + t114 * t194) * t277 + (-t115 * t117 + t116 * t118) * t278 + (t104 * t82 + t105 * t81 + t108 * t109) * t276 + (-0.2e1 * t240 + 0.2e1 * t241 + (-t104 * t171 - t105 * t168) * t215) * mrSges(7,3) + 0.2e1 * t281 * t115 + 0.2e1 * (-mrSges(4,1) * t169 - mrSges(4,2) * t172) * t247 + 0.2e1 * mrSges(6,3) * t194 + t175; t174 + m(7) * (t1 * t126 + t110 * t21 + t111 * t17 + t127 * t3 + t142 * t8 + t214 * t36) + (t167 * t25 - t239 * t24) * t271 + t31 * t232 - t1 * t168 * mrSges(7,3) - t77 * mrSges(5,3) * t202 - t258 * t263 - t21 * t211 + t61 * t214 - t171 * t213 - t17 * t209 + t110 * t63 + t111 * t65 + t126 * t28 + t127 * t30 + t142 * t19 + (m(6) * t24 + t20) * t155 + (-t219 * t27 - t259) * mrSges(6,3) + (m(6) * (-t259 + t261 + (-t168 * t27 - t171 * t26) * qJD(5)) - t168 * t29 - t64 * t219 - t66 * t218) * t154; t175 + m(7) * (t104 * t111 + t105 * t110 + t108 * t142 + t109 * t214 + t126 * t82 + t127 * t81) + t231 + t233 + t234 + t237 + t238 + t191 + (-t127 - t105) * t211 + (-t104 - t126) * t209 + (t167 * t271 - mrSges(5,2)) * t116 + (m(6) * t155 - t239 * t271 + t281) * t115 + (-mrSges(4,1) * t222 - mrSges(4,2) * t221) * pkin(2) + (-t240 - t235 + t236 + t241) * mrSges(7,3) + (m(6) * t154 + mrSges(6,3)) * t194; 0.2e1 * t231 + 0.2e1 * t191 + 0.2e1 * t233 + (t110 * t127 + t111 * t126 + t142 * t214) * t276 + (0.2e1 * t236 - 0.2e1 * t235 + (-t126 * t171 - t127 * t168) * t215) * mrSges(7,3) + t175; (t28 + t29) * t171 + (t30 + t31) * t168 + ((t63 + t64) * t171 + (-t65 - t66) * t168) * qJD(5) + m(6) * (t168 * t5 + t171 * t6 + (-t168 * t26 + t171 * t27) * qJD(5)) + m(7) * (t1 * t171 + t168 * t3 + (-t168 * t17 + t171 * t21) * qJD(5)) + m(5) * t96 + t201; m(7) * (t168 * t81 + t171 * t82 + (-t104 * t168 + t105 * t171) * qJD(5)); m(7) * (t110 * t168 + t111 * t171 + (-t126 * t168 + t127 * t171) * qJD(5)); 0; mrSges(6,1) * t6 + mrSges(7,1) * t1 - mrSges(6,2) * t5 - mrSges(7,2) * t3 - t77 * t200 + (m(7) * t1 + t28) * pkin(5) - t294 * t99 * qJD(5) + t195; -mrSges(7,2) * t81 + t212 * t82 - t187 * t116 + ((-mrSges(6,1) * t114 - t266) * t171 + (mrSges(6,2) * t114 - t290) * t168) * qJD(5) + t224; -mrSges(7,2) * t110 + t212 * t111 + ((-mrSges(6,1) * t154 - t266) * t171 + (mrSges(6,2) * t154 - t290) * t168) * qJD(5) + t224; (-t253 + (-mrSges(6,1) - t270) * t168) * qJD(5) - t136; 0; m(7) * t8 + t19; m(7) * t108 + t136; m(7) * t214 + t136; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
