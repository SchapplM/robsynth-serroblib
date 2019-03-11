% Calculate time derivative of joint inertia matrix for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:39:02
% EndTime: 2019-03-09 11:39:14
% DurationCPUTime: 5.24s
% Computational Cost: add. (7016->385), mult. (15236->546), div. (0->0), fcn. (15247->8), ass. (0->171)
t291 = Ifges(6,4) + Ifges(7,4);
t290 = Ifges(6,1) + Ifges(7,1);
t289 = Ifges(6,2) + Ifges(7,2);
t165 = sin(pkin(10));
t166 = cos(pkin(10));
t169 = sin(qJ(2));
t172 = cos(qJ(2));
t129 = -t165 * t169 + t166 * t172;
t130 = t165 * t172 + t166 * t169;
t168 = sin(qJ(4));
t171 = cos(qJ(4));
t178 = t171 * t129 - t130 * t168;
t240 = Ifges(6,5) + Ifges(7,5);
t288 = t178 * t240;
t277 = Ifges(6,6) + Ifges(7,6);
t287 = t178 * t277;
t170 = cos(qJ(5));
t286 = t291 * t170;
t167 = sin(qJ(5));
t285 = t291 * t167;
t283 = -t289 * t167 + t286;
t282 = t290 * t170 - t285;
t281 = t240 * t167 + t277 * t170;
t279 = t170 * t289 + t285;
t278 = t167 * t290 + t286;
t276 = Ifges(6,3) + Ifges(7,3);
t214 = qJD(5) * t167;
t121 = t130 * qJD(2);
t122 = t129 * qJD(2);
t70 = qJD(4) * t178 - t121 * t168 + t122 * t171;
t225 = t170 * t70;
t97 = t129 * t168 + t130 * t171;
t175 = t97 * t214 - t225;
t213 = qJD(5) * t170;
t205 = t97 * t213;
t228 = t167 * t70;
t176 = t205 + t228;
t71 = qJD(4) * t97 + t171 * t121 + t122 * t168;
t275 = -t175 * t291 - t176 * t289 + t277 * t71;
t274 = -t175 * t290 - t176 * t291 + t240 * t71;
t273 = t283 * t97 - t287;
t235 = t282 * t97 - t288;
t272 = -mrSges(6,1) * t170 + mrSges(6,2) * t167 - mrSges(5,1);
t271 = t283 * qJD(5);
t270 = t282 * qJD(5);
t217 = t240 * t213;
t177 = t270 * t167 + t271 * t170 + t213 * t278;
t269 = -t214 * t279 + t177;
t238 = -qJ(3) - pkin(7);
t143 = t238 * t169;
t145 = t238 * t172;
t98 = t166 * t143 + t145 * t165;
t93 = -pkin(8) * t130 + t98;
t99 = t165 * t143 - t166 * t145;
t94 = pkin(8) * t129 + t99;
t266 = -t168 * t94 + t171 * t93;
t52 = t168 * t93 + t171 * t94;
t49 = t170 * t52;
t155 = -pkin(2) * t172 - pkin(1);
t105 = -t129 * pkin(3) + t155;
t50 = -pkin(4) * t178 - t97 * pkin(9) + t105;
t31 = t167 * t50 + t49;
t265 = qJD(5) * t31;
t153 = pkin(2) * t166 + pkin(3);
t245 = pkin(2) * t165;
t114 = t153 * t171 - t168 * t245;
t115 = t168 * t153 + t171 * t245;
t263 = 2 * m(5);
t262 = 2 * m(6);
t261 = 2 * m(7);
t260 = -2 * mrSges(5,3);
t195 = qJD(2) * t238;
t118 = qJD(3) * t172 + t169 * t195;
t119 = -t169 * qJD(3) + t172 * t195;
t91 = -t118 * t165 + t119 * t166;
t174 = -pkin(8) * t122 + t91;
t92 = t166 * t118 + t165 * t119;
t77 = -pkin(8) * t121 + t92;
t25 = qJD(4) * t52 + t168 * t77 - t171 * t174;
t259 = 0.2e1 * t25;
t258 = -0.2e1 * t266;
t160 = qJD(2) * t169 * pkin(2);
t103 = pkin(3) * t121 + t160;
t257 = 0.2e1 * t103;
t134 = mrSges(7,1) * t214 + mrSges(7,2) * t213;
t256 = 0.2e1 * t134;
t255 = 0.2e1 * t155;
t254 = -0.2e1 * t167;
t253 = 0.2e1 * t170;
t252 = m(4) * pkin(2);
t251 = m(7) * pkin(5);
t248 = mrSges(7,3) * pkin(5);
t244 = pkin(5) * t170;
t24 = qJD(4) * t266 + t168 * t174 + t171 * t77;
t34 = pkin(4) * t71 - pkin(9) * t70 + t103;
t196 = -t167 * t24 + t170 * t34;
t6 = t196 - t265;
t243 = t167 * t6;
t211 = t167 * t34 + t170 * t24 + t50 * t213;
t5 = -t214 * t52 + t211;
t242 = t170 * t5;
t241 = t25 * t266;
t237 = -qJ(6) - pkin(9);
t234 = mrSges(6,2) * t170;
t109 = t115 * qJD(4);
t229 = t109 * t266;
t227 = t167 * t97;
t224 = t170 * t97;
t209 = pkin(5) * t214;
t102 = t109 + t209;
t141 = -mrSges(7,1) * t170 + mrSges(7,2) * t167;
t223 = t102 * t141;
t108 = t114 * qJD(4);
t222 = t108 * t170;
t112 = pkin(9) + t115;
t221 = t112 * t170;
t162 = t170 * qJ(6);
t219 = -qJ(6) - t112;
t216 = t167 ^ 2 + t170 ^ 2;
t215 = qJD(5) * t112;
t212 = 0.2e1 * t172;
t210 = 0.2e1 * qJD(5);
t208 = mrSges(7,1) + t251;
t30 = -t167 * t52 + t170 * t50;
t206 = t30 * t213;
t200 = t71 * mrSges(5,1) + t70 * mrSges(5,2);
t199 = t277 * t167;
t198 = -t214 / 0.2e1;
t194 = qJD(5) * t237;
t193 = t272 * t109;
t192 = t216 * mrSges(6,3);
t191 = t225 * t240 + t276 * t71;
t190 = t121 * mrSges(4,1) + t122 * mrSges(4,2);
t189 = t216 * t108;
t188 = qJD(5) * t219;
t111 = -pkin(4) - t114;
t185 = -(2 * Ifges(5,4)) - t199;
t184 = mrSges(6,1) * t167 + t234;
t179 = -qJ(6) * t70 - qJD(6) * t97;
t18 = mrSges(7,1) * t176 - mrSges(7,2) * t175;
t135 = t184 * qJD(5);
t3 = -qJ(6) * t205 + (-qJD(5) * t52 + t179) * t167 + t211;
t36 = pkin(5) * t227 - t266;
t8 = pkin(5) * t176 + t25;
t173 = -t24 * mrSges(5,2) + mrSges(6,3) * t242 + Ifges(5,5) * t70 + t36 * t134 - t266 * t135 + t8 * t141 + t272 * t25 - (-t214 * t277 + t217) * t178 / 0.2e1 + t274 * t167 / 0.2e1 - t271 * t227 / 0.2e1 + t270 * t224 / 0.2e1 + t273 * t198 + t235 * t213 / 0.2e1 + t279 * (-t205 / 0.2e1 - t228 / 0.2e1) + t278 * (t97 * t198 + t225 / 0.2e1) + (-Ifges(5,6) + t281 / 0.2e1) * t71 + (t3 * mrSges(7,3) + t275 / 0.2e1) * t170;
t161 = t170 * qJD(6);
t154 = -pkin(4) - t244;
t144 = pkin(9) * t170 + t162;
t140 = t237 * t167;
t117 = -qJD(6) * t167 + t170 * t194;
t116 = t167 * t194 + t161;
t104 = t111 - t244;
t101 = t162 + t221;
t100 = t219 * t167;
t79 = (-qJD(6) - t108) * t167 + t170 * t188;
t78 = t167 * t188 + t161 + t222;
t75 = -mrSges(6,1) * t178 - mrSges(6,3) * t224;
t74 = -mrSges(7,1) * t178 - mrSges(7,3) * t224;
t73 = mrSges(6,2) * t178 - mrSges(6,3) * t227;
t72 = mrSges(7,2) * t178 - mrSges(7,3) * t227;
t64 = t184 * t97;
t63 = (mrSges(7,1) * t167 + mrSges(7,2) * t170) * t97;
t29 = -mrSges(6,2) * t71 - mrSges(6,3) * t176;
t28 = -mrSges(7,2) * t71 - mrSges(7,3) * t176;
t27 = mrSges(6,1) * t71 + mrSges(6,3) * t175;
t26 = mrSges(7,1) * t71 + mrSges(7,3) * t175;
t21 = -qJ(6) * t227 + t31;
t19 = mrSges(6,1) * t176 - mrSges(6,2) * t175;
t17 = -pkin(5) * t178 - t162 * t97 + t30;
t1 = pkin(5) * t71 + t179 * t170 + (-t49 + (qJ(6) * t97 - t50) * t167) * qJD(5) + t196;
t2 = [-(mrSges(5,1) * t257 + t24 * t260 + ((2 * Ifges(5,2)) + t276) * t71 + t185 * t70 + t191) * t178 - 0.2e1 * t129 * Ifges(4,2) * t121 + t52 * t71 * t260 + 0.2e1 * t130 * t122 * Ifges(4,1) + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t172) * t212 + (0.2e1 * pkin(2) * (-mrSges(4,1) * t129 + mrSges(4,2) * t130) + t252 * t255 - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t169 + (Ifges(3,1) - Ifges(3,2)) * t212) * t169) * qJD(2) + 0.2e1 * t105 * t200 + t190 * t255 + t19 * t258 + t64 * t259 + (t1 * t17 + t21 * t3 + t36 * t8) * t261 + (t30 * t6 + t31 * t5 - t241) * t262 + (t103 * t105 + t24 * t52 - t241) * t263 + (mrSges(5,3) * t258 - t167 * t273 + t235 * t170) * t70 + 0.2e1 * t17 * t26 + 0.2e1 * t21 * t28 + 0.2e1 * t30 * t27 + 0.2e1 * t31 * t29 + 0.2e1 * t36 * t18 + 0.2e1 * t8 * t63 + 0.2e1 * t3 * t72 + 0.2e1 * t5 * t73 + 0.2e1 * t1 * t74 + 0.2e1 * t6 * t75 + 0.2e1 * (-t121 * t130 + t122 * t129) * Ifges(4,4) + 0.2e1 * (-t121 * t99 - t122 * t98 + t129 * t92 - t130 * t91) * mrSges(4,3) + 0.2e1 * m(4) * (t91 * t98 + t92 * t99) + (mrSges(5,2) * t257 + mrSges(5,3) * t259 + 0.2e1 * Ifges(5,1) * t70 + t274 * t170 - t275 * t167 + (t240 * t170 + t185) * t71 + ((-t273 + t287) * t170 + (-t235 + t288) * t167) * qJD(5)) * t97; m(6) * (t111 * t25 - t112 * t206 + t5 * t221 + t31 * t222 - t229) + (t108 * t73 + t112 * t29 + (-mrSges(6,3) * t30 - mrSges(7,3) * t17 - t112 * t75) * qJD(5)) * t170 + m(5) * (t108 * t52 - t114 * t25 + t115 * t24 - t229) + m(7) * (t1 * t100 + t101 * t3 + t102 * t36 + t104 * t8 + t17 * t79 + t21 * t78) + t173 + (t108 * t178 + t109 * t97 - t114 * t70 - t115 * t71) * mrSges(5,3) + (t165 * t92 + t166 * t91) * t252 + (-t121 * t165 - t122 * t166) * pkin(2) * mrSges(4,3) + (-t73 * t215 + m(6) * (-t108 * t30 - t112 * t6 - t215 * t31) - t108 * t75 - t112 * t27 + (-qJD(5) * t21 - t1) * mrSges(7,3) + (-t6 - t265) * mrSges(6,3)) * t167 + t78 * t72 + t79 * t74 + t91 * mrSges(4,1) - t92 * mrSges(4,2) + t100 * t26 + t101 * t28 + t102 * t63 + t104 * t18 + t109 * t64 + t111 * t19 - Ifges(4,6) * t121 + Ifges(4,5) * t122 + (Ifges(3,5) * t172 - Ifges(3,6) * t169 + (-mrSges(3,1) * t172 + mrSges(3,2) * t169) * pkin(7)) * qJD(2); 0.2e1 * t223 + t104 * t256 + 0.2e1 * t111 * t135 + (t109 * t111 + t112 * t189) * t262 - t109 * t114 * t263 + (t100 * t79 + t101 * t78 + t102 * t104) * t261 + (t79 * t254 + t78 * t253 + (-t100 * t170 - t101 * t167) * t210) * mrSges(7,3) + 0.2e1 * t193 + (t115 * t263 - 0.2e1 * mrSges(5,2) + 0.2e1 * t192) * t108 + t269; m(4) * t160 + (t26 + t27) * t170 + (t28 + t29) * t167 + ((t72 + t73) * t170 + (-t74 - t75) * t167) * qJD(5) + m(6) * (t167 * t5 + t170 * t6 + (-t167 * t30 + t170 * t31) * qJD(5)) + m(7) * (t1 * t170 + t167 * t3 + (-t167 * t17 + t170 * t21) * qJD(5)) + m(5) * t103 + t190 + t200; m(7) * (t167 * t78 + t170 * t79 + (-t100 * t167 + t101 * t170) * qJD(5)); 0; t63 * t209 + (-t243 + (-t167 * t31 - t170 * t30) * qJD(5)) * mrSges(6,3) + (-t1 * t167 + (-t167 * t21 - t17 * t170) * qJD(5)) * mrSges(7,3) + m(7) * (t1 * t140 + t116 * t21 + t117 * t17 + t144 * t3 + t154 * t8 + t209 * t36) + t173 + (-t75 * t213 - t73 * t214 + m(6) * (-t214 * t31 - t206 + t242 - t243) + t170 * t29 - t167 * t27) * pkin(9) + t116 * t72 + t117 * t74 + t140 * t26 + t144 * t28 + t154 * t18 + (-m(6) * t25 - t19) * pkin(4); t223 + (-pkin(4) + t111) * t135 + (t104 + t154) * t134 + t193 + (-mrSges(5,2) + t192) * t108 + (pkin(5) * t141 - t279) * t214 + m(7) * (t100 * t117 + t101 * t116 + t102 * t154 + t104 * t209 + t140 * t79 + t144 * t78) + m(6) * (-pkin(4) * t109 + pkin(9) * t189) + ((t116 + t78) * t170 + (-t117 - t79) * t167 + ((-t100 - t140) * t170 + (-t101 - t144) * t167) * qJD(5)) * mrSges(7,3) + t177; m(7) * (t116 * t167 + t117 * t170 + (-t140 * t167 + t144 * t170) * qJD(5)); -0.2e1 * pkin(4) * t135 + 0.2e1 * t141 * t209 + t154 * t256 + (t116 * t144 + t117 * t140 + t154 * t209) * t261 + (t116 * t253 + t117 * t254 + (-t140 * t170 - t144 * t167) * t210) * mrSges(7,3) + t269; mrSges(6,1) * t6 + mrSges(7,1) * t1 - mrSges(6,2) * t5 - mrSges(7,2) * t3 - t70 * t199 + (m(7) * t1 + t26) * pkin(5) - t281 * t97 * qJD(5) + t191; -mrSges(7,2) * t78 + t208 * t79 - t184 * t108 + ((-mrSges(6,1) * t112 - t248) * t170 + (mrSges(6,2) * t112 - t277) * t167) * qJD(5) + t217; (-t234 + (-mrSges(6,1) - t251) * t167) * qJD(5) - t134; -mrSges(7,2) * t116 + t208 * t117 + ((-mrSges(6,1) * pkin(9) - t248) * t170 + (mrSges(6,2) * pkin(9) - t277) * t167) * qJD(5) + t217; 0; m(7) * t8 + t18; m(7) * t102 + t134; 0; m(7) * t209 + t134; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
