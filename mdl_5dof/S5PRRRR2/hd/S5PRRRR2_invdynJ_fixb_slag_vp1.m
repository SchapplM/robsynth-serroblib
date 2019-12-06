% Calculate vector of inverse dynamics joint torques for
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:43
% EndTime: 2019-12-05 17:04:49
% DurationCPUTime: 3.14s
% Computational Cost: add. (7881->323), mult. (5728->422), div. (0->0), fcn. (4324->8), ass. (0->197)
t160 = qJ(2) + qJ(3);
t155 = qJ(4) + t160;
t146 = cos(t155);
t159 = qJD(2) + qJD(3);
t151 = qJD(4) + t159;
t231 = t146 * t151;
t114 = pkin(6) * t231;
t161 = sin(qJ(5));
t163 = cos(qJ(5));
t248 = rSges(6,2) * t163;
t130 = rSges(6,1) * t161 + t248;
t223 = qJD(5) * t146;
t145 = sin(t155);
t233 = t145 * t161;
t120 = rSges(6,2) * t233;
t232 = t145 * t163;
t74 = rSges(6,1) * t232 - t146 * rSges(6,3) - t120;
t66 = t151 * t74;
t274 = t130 * t223 - t114 + t66;
t162 = sin(qJ(2));
t247 = pkin(2) * qJD(2);
t217 = t162 * t247;
t152 = sin(t160);
t228 = t152 * t159;
t220 = pkin(3) * t228;
t100 = rSges(5,1) * t145 + rSges(5,2) * t146;
t92 = t151 * t100;
t60 = -t217 - t220 - t92;
t273 = t220 + t274;
t153 = cos(t160);
t104 = rSges(4,1) * t152 + rSges(4,2) * t153;
t238 = t104 * t159;
t80 = -t217 - t238;
t249 = rSges(6,1) * t163;
t132 = -rSges(6,2) * t161 + t249;
t111 = t132 * qJD(5);
t149 = t151 ^ 2;
t158 = qJDD(2) + qJDD(3);
t150 = qJDD(4) + t158;
t157 = t159 ^ 2;
t164 = cos(qJ(2));
t165 = qJD(2) ^ 2;
t187 = (-qJDD(2) * t162 - t164 * t165) * pkin(2);
t172 = t187 + (-t152 * t158 - t153 * t157) * pkin(3);
t214 = qJD(5) * t248;
t230 = t146 * t161;
t218 = rSges(6,2) * t230;
t221 = qJD(5) * t161;
t213 = t151 * t218 + (rSges(6,1) * t221 + t214) * t145;
t229 = t146 * t163;
t121 = rSges(6,1) * t229;
t265 = t145 * rSges(6,3) + t121;
t45 = t265 * t151 - t213;
t222 = qJD(5) * t151;
t85 = -qJDD(5) * t146 + t145 * t222;
t14 = -t111 * t223 + t130 * t85 - t150 * t74 - t151 * t45 + (-t145 * t149 + t146 * t150) * pkin(6) + t172;
t272 = -g(1) + t14;
t234 = t145 * t151;
t79 = rSges(5,1) * t231 - rSges(5,2) * t234;
t271 = -t100 * t150 - t151 * t79 - g(1) + t172;
t227 = t153 * t159;
t96 = rSges(4,1) * t227 - rSges(4,2) * t228;
t270 = -t104 * t158 - t159 * t96 - g(1) + t187;
t144 = pkin(3) * t153;
t156 = t164 * pkin(2);
t257 = pkin(2) * t162;
t204 = qJDD(2) * t156 - t165 * t257;
t256 = pkin(3) * t152;
t181 = t158 * t144 - t157 * t256 + t204;
t224 = qJD(5) * t145;
t191 = rSges(6,3) * t231 + t151 * t120 - t146 * t214;
t212 = t146 * t221;
t44 = (-t151 * t232 - t212) * rSges(6,1) + t191;
t75 = -t218 + t265;
t84 = qJDD(5) * t145 + t146 * t222;
t15 = -t111 * t224 - t130 * t84 + t150 * t75 + t151 * t44 + (t145 * t150 + t146 * t149) * pkin(6) + t181;
t269 = -g(2) + t15;
t101 = t146 * rSges(5,1) - rSges(5,2) * t145;
t268 = t101 * t150 - t151 * t92 - g(2) + t181;
t105 = t153 * rSges(4,1) - rSges(4,2) * t152;
t267 = t105 * t158 - t159 * t238 - g(2) + t204;
t266 = t145 * pkin(6) + t75;
t32 = -t217 - t273;
t154 = Icges(6,4) * t163;
t195 = -Icges(6,2) * t161 + t154;
t127 = Icges(6,1) * t161 + t154;
t124 = Icges(6,5) * t163 - Icges(6,6) * t161;
t123 = Icges(6,5) * t161 + Icges(6,6) * t163;
t177 = Icges(6,3) * t151 - qJD(5) * t123;
t189 = t195 * t146;
t71 = Icges(6,6) * t145 + t189;
t244 = t161 * t71;
t241 = Icges(6,4) * t161;
t128 = Icges(6,1) * t163 - t241;
t190 = t128 * t146;
t73 = Icges(6,5) * t145 + t190;
t197 = -t163 * t73 + t244;
t264 = -t124 * t234 + t146 * t177 + t151 * t197;
t188 = t124 * t146;
t70 = Icges(6,4) * t232 - Icges(6,2) * t233 - Icges(6,6) * t146;
t245 = t161 * t70;
t119 = Icges(6,4) * t233;
t72 = Icges(6,1) * t232 - Icges(6,5) * t146 - t119;
t198 = -t163 * t72 + t245;
t263 = t145 * t177 + (t188 + t198) * t151;
t192 = t130 * t224 - t151 * t266;
t125 = Icges(6,2) * t163 + t241;
t193 = t125 * t161 - t127 * t163;
t262 = t124 * qJD(5) + t151 * t193;
t68 = Icges(6,5) * t232 - Icges(6,6) * t233 - Icges(6,3) * t146;
t24 = -t145 * t198 - t146 * t68;
t216 = t164 * t247;
t219 = pkin(3) * t227;
t185 = t216 + t219;
t33 = t185 - t192;
t168 = (-t32 * t121 + (t32 * (-rSges(6,3) - pkin(6)) - t33 * t249) * t145) * t151;
t261 = t168 + (t213 - t192) * t32;
t251 = -Icges(6,2) * t232 - t119 + t72;
t253 = t127 * t145 + t70;
t260 = -t161 * t251 - t163 * t253;
t259 = t84 / 0.2e1;
t258 = t85 / 0.2e1;
t255 = -t145 * t68 - t72 * t229;
t69 = Icges(6,3) * t145 + t188;
t254 = t145 * t69 + t73 * t229;
t252 = -t127 * t146 - t71;
t250 = -t125 * t146 + t73;
t236 = t123 * t146;
t49 = -t145 * t193 - t236;
t243 = t49 * t151;
t239 = t101 * t151;
t237 = t123 * t145;
t235 = t124 * t151;
t226 = -t125 + t128;
t225 = t127 + t195;
t215 = m(2) + m(3) + m(4) + m(5);
t210 = -t224 / 0.2e1;
t209 = t224 / 0.2e1;
t208 = -t223 / 0.2e1;
t207 = t223 / 0.2e1;
t57 = t73 * t232;
t206 = t146 * t69 - t57;
t205 = -t68 + t244;
t83 = t101 + t144;
t133 = rSges(3,1) * t164 - rSges(3,2) * t162;
t131 = rSges(3,1) * t162 + rSges(3,2) * t164;
t25 = -t233 * t71 - t206;
t202 = t145 * t25 - t146 * t24;
t26 = -t230 * t70 - t255;
t27 = -t230 * t71 + t254;
t201 = t145 * t27 - t146 * t26;
t200 = -t145 * t33 - t146 * t32;
t199 = t145 * t74 + t146 * t75;
t46 = t161 * t72 + t163 * t70;
t47 = t161 * t73 + t163 * t71;
t194 = t125 * t163 + t127 * t161;
t63 = t146 * pkin(6) - t74;
t82 = -t100 - t256;
t56 = t144 + t266;
t184 = -t161 * t250 + t163 * t252;
t183 = -t79 - t219;
t55 = t63 - t256;
t180 = (-t161 * t225 + t163 * t226) * t151;
t179 = Icges(6,5) * t151 - qJD(5) * t127;
t178 = Icges(6,6) * t151 - qJD(5) * t125;
t174 = -rSges(6,1) * t212 + t114 + t191;
t50 = -t146 * t193 + t237;
t48 = t50 * t151;
t10 = qJD(5) * t201 + t48;
t107 = t195 * qJD(5);
t108 = t128 * qJD(5);
t40 = t145 * t178 + t151 * t189;
t42 = t145 * t179 + t151 * t190;
t16 = -qJD(5) * t198 + t161 * t42 + t163 * t40;
t39 = t146 * t178 - t195 * t234;
t41 = -t128 * t234 + t146 * t179;
t17 = -qJD(5) * t197 + t161 * t41 + t163 * t39;
t167 = -qJD(5) * t194 - t107 * t161 + t108 * t163 + t123 * t151;
t20 = t262 * t145 + t167 * t146;
t21 = t167 * t145 - t262 * t146;
t9 = qJD(5) * t202 + t243;
t173 = (t48 + ((t25 - t57 + (t69 + t245) * t146 + t255) * t146 + t254 * t145) * qJD(5)) * t207 + (-qJD(5) * t193 + t107 * t163 + t108 * t161) * t151 + (t47 + t50) * t259 + (t46 + t49) * t258 + (-t243 + ((t146 * t205 - t254 + t27) * t146 + (t145 * t205 + t206 + t26) * t145) * qJD(5) + t9) * t210 + (t17 + t20) * t209 + (Icges(5,3) + t194) * t150 + (t16 + t21 + t10) * t208;
t171 = Icges(4,3) * t158 + t173;
t170 = -qJD(5) * t46 + t151 * t68 - t161 * t40 + t163 * t42;
t169 = -qJD(5) * t47 + t151 * t69 - t161 * t39 + t163 * t41;
t166 = t174 - t220;
t94 = t130 * t146;
t93 = t130 * t145;
t81 = t105 * t159 + t216;
t61 = t185 + t239;
t36 = qJD(5) * t199 + qJD(1);
t11 = t74 * t84 - t75 * t85 + qJDD(1) + (t145 * t45 + t146 * t44) * qJD(5);
t6 = t169 * t145 - t264 * t146;
t5 = t170 * t145 - t263 * t146;
t4 = t264 * t145 + t169 * t146;
t3 = t263 * t145 + t170 * t146;
t1 = [m(6) * t11 + t215 * qJDD(1) + (-m(6) - t215) * g(3); Icges(3,3) * qJDD(2) + t171 + (t267 * (t105 + t156) + t270 * (-t104 - t257) + (-t96 - t216 + t81) * t80) * m(4) + (g(1) * t131 - g(2) * t133 + (t131 ^ 2 + t133 ^ 2) * qJDD(2)) * m(3) + (t32 * (-t185 + t213) + t168 + t269 * (t156 + t56) + t272 * (t55 - t257) + t33 * (t166 - t217)) * m(6) + (t268 * (t156 + t83) + t271 * (t82 - t257) + (t61 + t183 - t216) * t60) * m(5); t171 + (t269 * t56 + t272 * t55 + (t166 + t273) * t33 + t261) * m(6) + (t268 * t83 + t271 * t82 + (t183 + t219 + t239) * t60) * m(5) + (-t80 * t96 - t81 * t238 + (t80 * t159 + t267) * t105 + (t81 * t159 - t270) * t104) * m(4); t173 + (t272 * t63 + t269 * t266 + (t174 + t274) * t33 + t261) * m(6) + (-t60 * t79 - t61 * t92 + (t151 * t60 + t268) * t101 + (t151 * t61 - t271) * t100) * m(5); t10 * t231 / 0.2e1 + t145 * (t150 * t50 + t151 * t20 + t26 * t85 + t27 * t84 + (t145 * t4 - t146 * t3) * qJD(5)) / 0.2e1 + t201 * t259 + ((t151 * t27 - t3) * t146 + (t151 * t26 + t4) * t145) * t209 + t9 * t234 / 0.2e1 - t146 * (t150 * t49 + t151 * t21 + t24 * t85 + t25 * t84 + (t145 * t6 - t146 * t5) * qJD(5)) / 0.2e1 + t202 * t258 + ((t151 * t25 - t5) * t146 + (t151 * t24 + t6) * t145) * t208 + t150 * (t145 * t47 - t146 * t46) / 0.2e1 + t151 * ((t151 * t47 - t16) * t146 + (t151 * t46 + t17) * t145) / 0.2e1 + ((-t224 * t236 + t235) * t145 + (t180 + (-t260 * t146 + (t237 + t184) * t145) * qJD(5)) * t146) * t210 + ((-t223 * t237 - t235) * t146 + (t180 + (t184 * t145 + (-t260 + t236) * t146) * qJD(5)) * t145) * t207 - t151 * ((t161 * t226 + t163 * t225) * t151 + ((t145 * t250 - t146 * t251) * t163 + (t145 * t252 + t146 * t253) * t161) * qJD(5)) / 0.2e1 + (t11 * t199 + t36 * ((t44 + t66) * t146 + (-t151 * t75 + t45) * t145) + t200 * t111 + ((-t151 * t33 - t14) * t146 + (t151 * t32 - t15) * t145) * t130 - (t32 * t93 - t33 * t94) * t151 - (t36 * (-t145 * t93 - t146 * t94) + t200 * t132) * qJD(5) + g(1) * t94 + g(2) * t93 - g(3) * t132) * m(6);];
tau = t1;
