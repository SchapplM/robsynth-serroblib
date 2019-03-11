% Calculate time derivative of joint inertia matrix for
% S6RRRPRP2
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
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:35:09
% EndTime: 2019-03-09 16:35:19
% DurationCPUTime: 4.58s
% Computational Cost: add. (7262->374), mult. (16083->533), div. (0->0), fcn. (15487->8), ass. (0->181)
t299 = Ifges(6,1) + Ifges(7,1);
t257 = Ifges(7,4) + Ifges(6,5);
t166 = sin(qJ(3));
t167 = sin(qJ(2));
t169 = cos(qJ(3));
t170 = cos(qJ(2));
t221 = t169 * t170;
t129 = -t166 * t167 + t221;
t130 = t166 * t170 + t169 * t167;
t164 = sin(pkin(10));
t237 = cos(pkin(10));
t99 = -t129 * t237 + t130 * t164;
t298 = t257 * t99;
t165 = sin(qJ(5));
t168 = cos(qJ(5));
t297 = t165 ^ 2 + t168 ^ 2;
t246 = Ifges(7,5) * t165;
t248 = Ifges(6,4) * t165;
t296 = t299 * t168 + t246 - t248;
t245 = Ifges(7,5) * t168;
t247 = Ifges(6,4) * t168;
t295 = t165 * t299 - t245 + t247;
t215 = qJD(5) * t168;
t216 = qJD(5) * t165;
t268 = -pkin(8) - pkin(7);
t147 = t268 * t167;
t148 = t268 * t170;
t131 = t166 * t147;
t174 = (t221 * t268 - t131) * qJD(2);
t280 = qJD(2) + qJD(3);
t103 = t280 * t129;
t179 = -t103 * qJ(4) - t130 * qJD(4);
t217 = qJD(3) * t169;
t218 = qJD(3) * t166;
t104 = t280 * t130;
t79 = qJD(2) * t130 * t268 + t147 * t217 + t148 * t218;
t47 = -qJ(4) * t104 + qJD(4) * t129 + t79;
t27 = t237 * t47 + (-t147 * t218 + t148 * t217 + t174 + t179) * t164;
t76 = t103 * t164 + t104 * t237;
t77 = t103 * t237 - t164 * t104;
t96 = qJD(2) * t167 * pkin(2) + pkin(3) * t104;
t35 = pkin(4) * t76 - pkin(9) * t77 + t96;
t100 = t164 * t129 + t130 * t237;
t154 = -pkin(2) * t170 - pkin(1);
t110 = -t129 * pkin(3) + t154;
t52 = t99 * pkin(4) - t100 * pkin(9) + t110;
t108 = t169 * t147 + t148 * t166;
t176 = -qJ(4) * t130 + t108;
t109 = -t169 * t148 + t131;
t92 = qJ(4) * t129 + t109;
t54 = t164 * t176 + t237 * t92;
t6 = t165 * t35 + t168 * t27 + t52 * t215 - t216 * t54;
t2 = qJ(6) * t76 + qJD(6) * t99 + t6;
t252 = t165 * t52 + t168 * t54;
t281 = qJD(5) * t252;
t7 = -t165 * t27 + t168 * t35 - t281;
t4 = -pkin(5) * t76 - t7;
t192 = t165 * t4 + t168 * t2;
t28 = -t165 * t54 + t168 * t52;
t25 = -pkin(5) * t99 - t28;
t294 = t25 * t215 + t192;
t293 = Ifges(6,6) * t168 + t257 * t165;
t292 = mrSges(7,2) + mrSges(6,3);
t291 = Ifges(7,2) + Ifges(6,3);
t239 = t168 * t77;
t177 = t100 * t216 - t239;
t203 = t100 * t215;
t240 = t165 * t77;
t178 = t203 + t240;
t290 = t257 * t76 + (-Ifges(6,4) + Ifges(7,5)) * t178 - t299 * t177;
t253 = t296 * t100 + t298;
t141 = -t168 * mrSges(6,1) + t165 * mrSges(6,2);
t284 = -mrSges(5,1) + t141;
t289 = t296 * qJD(5);
t288 = Ifges(7,6) * t216 + t257 * t215;
t224 = t164 * t166;
t242 = pkin(2) * qJD(3);
t114 = (t169 * t237 - t224) * t242;
t286 = t297 * t114;
t32 = -mrSges(6,2) * t76 - mrSges(6,3) * t178;
t33 = -mrSges(7,2) * t178 + mrSges(7,3) * t76;
t255 = t32 + t33;
t30 = mrSges(6,1) * t76 + mrSges(6,3) * t177;
t31 = -t76 * mrSges(7,1) - mrSges(7,2) * t177;
t256 = t30 - t31;
t285 = -t256 * t165 + t255 * t168;
t236 = t100 * t165;
t63 = -mrSges(6,2) * t99 - mrSges(6,3) * t236;
t66 = -mrSges(7,2) * t236 + mrSges(7,3) * t99;
t251 = t63 + t66;
t235 = t100 * t168;
t64 = mrSges(6,1) * t99 - mrSges(6,3) * t235;
t65 = -mrSges(7,1) * t99 + mrSges(7,2) * t235;
t250 = -t64 + t65;
t283 = -t28 * t215 - t252 * t216;
t279 = 2 * m(5);
t278 = 2 * m(6);
t277 = 2 * m(7);
t276 = -2 * Ifges(5,4);
t80 = -qJD(3) * t109 + t174;
t26 = t164 * t47 - t237 * (t179 + t80);
t275 = 0.2e1 * t26;
t53 = t164 * t92 - t237 * t176;
t274 = 0.2e1 * t53;
t273 = 0.2e1 * t96;
t272 = 0.2e1 * t154;
t271 = m(5) * pkin(3);
t270 = t76 / 0.2e1;
t143 = Ifges(6,2) * t168 + t248;
t267 = -t143 / 0.2e1;
t265 = Ifges(6,6) * t99;
t264 = pkin(3) * t164;
t261 = t168 * t6;
t260 = t26 * t53;
t259 = t7 * t165;
t258 = t76 * mrSges(5,3);
t185 = Ifges(7,3) * t165 + t245;
t42 = Ifges(7,6) * t99 + t100 * t185;
t186 = -Ifges(6,2) * t165 + t247;
t43 = t100 * t186 + t265;
t254 = t42 - t43;
t153 = pkin(2) * t169 + pkin(3);
t196 = t237 * t166;
t116 = pkin(2) * t196 + t164 * t153;
t112 = pkin(9) + t116;
t249 = t286 * t112;
t244 = Ifges(6,6) * t165;
t113 = (t164 * t169 + t196) * t242;
t241 = t113 * t53;
t140 = -t168 * mrSges(7,1) - t165 * mrSges(7,3);
t117 = -pkin(5) * t216 + qJ(6) * t215 + t165 * qJD(6);
t97 = t113 - t117;
t238 = t97 * t140;
t115 = -pkin(2) * t224 + t153 * t237;
t111 = -pkin(4) - t115;
t184 = t168 * pkin(5) + t165 * qJ(6);
t107 = t111 - t184;
t189 = t165 * mrSges(7,1) - t168 * mrSges(7,3);
t134 = t189 * qJD(5);
t234 = t107 * t134;
t190 = t165 * mrSges(6,1) + t168 * mrSges(6,2);
t135 = t190 * qJD(5);
t233 = t111 * t135;
t230 = t114 * t165;
t229 = t114 * t168;
t205 = t237 * pkin(3);
t152 = -t205 - pkin(4);
t124 = -t184 + t152;
t228 = t124 * t134;
t225 = t152 * t135;
t151 = pkin(9) + t264;
t220 = t286 * t151;
t214 = qJD(6) * t168;
t213 = 0.2e1 * t170;
t212 = m(7) * t214;
t24 = qJ(6) * t99 + t252;
t211 = t24 * t216;
t201 = t151 * t215;
t199 = t76 * mrSges(5,1) + t77 * mrSges(5,2);
t198 = -t216 / 0.2e1;
t195 = mrSges(7,2) * t214 + t288;
t191 = -t259 + t261;
t183 = pkin(5) * t165 - qJ(6) * t168;
t180 = Ifges(7,6) * t178 + t239 * t257 + t291 * t76;
t175 = -mrSges(7,2) * t184 - t244;
t136 = t185 * qJD(5);
t137 = t186 * qJD(5);
t142 = -Ifges(7,3) * t168 + t246;
t173 = (t142 - t143) * t216 + t295 * t215 + (-t136 + t137) * t168 + t289 * t165;
t172 = -m(7) * t184 + t140 + t141;
t14 = -Ifges(7,5) * t177 + Ifges(7,6) * t76 + Ifges(7,3) * t178;
t15 = -Ifges(6,4) * t177 - Ifges(6,2) * t178 + Ifges(6,6) * t76;
t37 = t100 * t183 + t53;
t9 = t183 * t77 + (qJD(5) * t184 - t214) * t100 + t26;
t171 = t43 * t198 + (t267 + t142 / 0.2e1) * t240 + (-t137 / 0.2e1 + t136 / 0.2e1) * t236 + t295 * (t100 * t198 + t239 / 0.2e1) + t294 * mrSges(7,2) + t293 * t270 + (-Ifges(7,6) * t270 + t15 / 0.2e1 - t14 / 0.2e1) * t168 + mrSges(6,3) * t261 + t203 * t267 + t42 * t216 / 0.2e1 + (t100 * t142 + t253) * t215 / 0.2e1 + t290 * t165 / 0.2e1 - t27 * mrSges(5,2) - Ifges(5,6) * t76 + Ifges(5,5) * t77 - t79 * mrSges(4,2) + t80 * mrSges(4,1) + Ifges(4,5) * t103 - Ifges(4,6) * t104 + t37 * t134 + t53 * t135 + t9 * t140 + (-Ifges(6,6) * t216 + t288) * t99 / 0.2e1 + t289 * t235 / 0.2e1 + t284 * t26;
t156 = mrSges(7,2) * t215;
t62 = t190 * t100;
t61 = t189 * t100;
t21 = mrSges(6,1) * t178 - mrSges(6,2) * t177;
t20 = mrSges(7,1) * t178 + mrSges(7,3) * t177;
t1 = [-0.2e1 * t129 * Ifges(4,2) * t104 + 0.2e1 * t130 * t103 * Ifges(4,1) + t21 * t274 + t62 * t275 + (t2 * t24 + t25 * t4 + t37 * t9) * t277 + (t110 * t96 + t27 * t54 + t260) * t279 + (mrSges(5,2) * t273 + mrSges(5,3) * t275 + 0.2e1 * Ifges(5,1) * t77 + t290 * t168 + (t14 - t15) * t165 + (t276 + t257 * t168 + (-Ifges(6,6) + Ifges(7,6)) * t165) * t76 + ((t254 - t265) * t168 + (-t253 - t298) * t165) * qJD(5)) * t100 + (mrSges(4,1) * t104 + mrSges(4,2) * t103) * t272 + 0.2e1 * (t103 * t129 - t104 * t130) * Ifges(4,4) + 0.2e1 * (-t103 * t108 - t104 * t109 + t129 * t79 - t130 * t80) * mrSges(4,3) + (t252 * t6 + t28 * t7 + t260) * t278 + 0.2e1 * t252 * t32 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t170) * t213 + (m(4) * pkin(2) * t272 + 0.2e1 * pkin(2) * (-mrSges(4,1) * t129 + mrSges(4,2) * t130) - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t167 + (Ifges(3,1) - Ifges(3,2)) * t213) * t167) * qJD(2) + (mrSges(5,3) * t274 + t254 * t165 + t253 * t168) * t77 + 0.2e1 * m(4) * (t108 * t80 + t109 * t79) - 0.2e1 * t54 * t258 + (mrSges(5,1) * t273 - 0.2e1 * mrSges(5,3) * t27 + (t276 - t244) * t77 + ((2 * Ifges(5,2)) + t291) * t76 + t180) * t99 + 0.2e1 * t28 * t30 + 0.2e1 * t25 * t31 + 0.2e1 * t24 * t33 + 0.2e1 * t37 * t20 + 0.2e1 * t9 * t61 + 0.2e1 * t6 * t63 + 0.2e1 * t7 * t64 + 0.2e1 * t4 * t65 + 0.2e1 * t2 * t66 + 0.2e1 * t110 * t199; t171 + (t100 * t113 - t114 * t99 - t115 * t77 - t116 * t76) * mrSges(5,3) + (m(4) * (t166 * t79 + t169 * t80 + (-t108 * t166 + t109 * t169) * qJD(3)) + (-t169 * t103 - t166 * t104 + (t129 * t169 + t130 * t166) * qJD(3)) * mrSges(4,3)) * pkin(2) + ((-t165 * t251 + t168 * t250) * qJD(5) + m(6) * (t191 + t283) + m(7) * (-t211 + t294) + t285) * t112 + m(7) * (t107 * t9 + t24 * t229 + t25 * t230 + t37 * t97) + m(6) * (t111 * t26 + t229 * t252 - t28 * t230 + t241) + m(5) * (t114 * t54 - t115 * t26 + t116 * t27 + t241) + t97 * t61 + t107 * t20 + t111 * t21 + t113 * t62 + (-t28 * qJD(5) * mrSges(6,3) + t251 * t114) * t168 + (-t24 * qJD(5) * mrSges(7,2) + t250 * t114 + (-t7 - t281) * mrSges(6,3)) * t165 + (Ifges(3,5) * t170 - Ifges(3,6) * t167 + (-mrSges(3,1) * t170 + mrSges(3,2) * t167) * pkin(7)) * qJD(2); 0.2e1 * (t292 * t297 - mrSges(5,2)) * t114 + (t107 * t97 + t249) * t277 + (t111 * t113 + t249) * t278 + (-t113 * t115 + t114 * t116) * t279 + 0.2e1 * t234 + 0.2e1 * t233 + 0.2e1 * t238 + t173 + 0.2e1 * t284 * t113 + 0.2e1 * (-mrSges(4,1) * t166 - mrSges(4,2) * t169) * t242; -t77 * mrSges(5,3) * t205 + (t164 * t27 - t237 * t26) * t271 + t171 + m(7) * (-t117 * t37 + t124 * t9) - t258 * t264 - mrSges(7,2) * t211 - t117 * t61 + t124 * t20 + t250 * t201 + (-t259 + t283) * mrSges(6,3) + (m(6) * t26 + t21) * t152 + (m(7) * ((-t165 * t24 + t168 * t25) * qJD(5) + t192) + m(6) * ((-t165 * t252 - t168 * t28) * qJD(5) + t191) - t251 * t216 + t285) * t151; m(6) * t220 + m(7) * (-t107 * t117 + t124 * t97 + t220) + t228 + t234 + t233 - t117 * t140 + t238 + t225 + t173 + (t164 * t271 - mrSges(5,2)) * t114 + (m(6) * t152 - t237 * t271 + t284) * t113 + (-mrSges(4,1) * t218 - mrSges(4,2) * t217) * pkin(2) + t292 * t286; 0.2e1 * t228 + 0.2e1 * t225 + 0.2e1 * (-m(7) * t124 - t140) * t117 + t173; t256 * t168 + t255 * t165 + (t250 * t165 + t251 * t168) * qJD(5) + m(7) * (t165 * t2 - t168 * t4 + (t165 * t25 + t168 * t24) * qJD(5)) + m(6) * (t165 * t6 + t168 * t7 + (-t165 * t28 + t168 * t252) * qJD(5)) + m(5) * t96 + t199; 0; 0; 0; -Ifges(6,6) * t240 + m(7) * (-pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t24) + t2 * mrSges(7,3) + qJD(6) * t66 + qJ(6) * t33 - t6 * mrSges(6,2) + t7 * mrSges(6,1) - t4 * mrSges(7,1) - pkin(5) * t31 - t293 * t100 * qJD(5) + t180; t112 * t212 + (-m(7) * t183 - t189 - t190) * t114 + (t112 * t172 + t175) * qJD(5) + t195; t175 * qJD(5) + (qJD(5) * t172 + t212) * t151 + t195; m(7) * t117 + ((-mrSges(6,2) + mrSges(7,3)) * t168 + (-mrSges(6,1) - mrSges(7,1)) * t165) * qJD(5); 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t4 + t31; t156 + (t112 * t215 + t230) * m(7); m(7) * t201 + t156; m(7) * t216; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
