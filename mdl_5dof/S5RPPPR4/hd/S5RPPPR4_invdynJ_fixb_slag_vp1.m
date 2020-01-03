% Calculate vector of inverse dynamics joint torques for
% S5RPPPR4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:06
% EndTime: 2019-12-31 17:45:13
% DurationCPUTime: 5.04s
% Computational Cost: add. (5408->379), mult. (4897->469), div. (0->0), fcn. (3546->8), ass. (0->196)
t170 = cos(qJ(1));
t163 = t170 * pkin(1);
t165 = qJ(1) + pkin(7);
t160 = sin(t165);
t162 = cos(t165);
t108 = rSges(3,1) * t160 + rSges(3,2) * t162;
t169 = sin(qJ(1));
t276 = pkin(1) * t169;
t89 = -t108 - t276;
t171 = qJD(1) ^ 2;
t295 = t171 * t163;
t111 = -rSges(4,2) * t162 + t160 * rSges(4,3);
t288 = t162 * pkin(2) + t160 * qJ(3);
t217 = t163 + t288;
t62 = t111 + t217;
t164 = pkin(8) + qJ(5);
t159 = sin(t164);
t161 = cos(t164);
t202 = rSges(6,1) * t159 + rSges(6,2) * t161;
t247 = t160 * t161;
t120 = Icges(6,4) * t247;
t248 = t159 * t160;
t253 = Icges(6,5) * t162;
t67 = Icges(6,1) * t248 + t120 + t253;
t254 = Icges(6,4) * t161;
t192 = Icges(6,1) * t159 + t254;
t68 = -Icges(6,5) * t160 + t162 * t192;
t102 = -Icges(6,2) * t159 + t254;
t79 = t102 * t162;
t177 = t160 * (t68 + t79) - t162 * (-Icges(6,2) * t248 + t120 + t67);
t255 = Icges(6,4) * t159;
t191 = Icges(6,2) * t161 + t255;
t65 = Icges(6,6) * t162 + t160 * t191;
t66 = -Icges(6,6) * t160 + t162 * t191;
t104 = Icges(6,1) * t161 - t255;
t80 = t104 * t160;
t81 = t104 * t162;
t178 = t160 * (t66 - t81) - t162 * (t65 - t80);
t294 = -t178 * t159 + t177 * t161;
t239 = t102 + t192;
t240 = -t191 + t104;
t293 = (t159 * t239 - t161 * t240) * qJD(1);
t195 = t159 * t68 + t161 * t66;
t290 = t195 * t162;
t112 = t162 * rSges(3,1) - rSges(3,2) * t160;
t90 = t112 + t163;
t166 = sin(pkin(8));
t246 = t160 * t166;
t131 = pkin(4) * t246;
t168 = -pkin(6) - qJ(4);
t287 = t162 * t168 - t131;
t227 = qJD(1) * qJD(3);
t232 = qJD(1) * t160;
t231 = qJD(1) * t162;
t137 = qJ(3) * t231;
t144 = qJD(3) * t160;
t235 = t137 + t144;
t286 = t160 * t227 + qJD(1) * (-pkin(2) * t232 + t235) + qJDD(1) * t288;
t167 = cos(pkin(8));
t260 = rSges(5,2) * t167;
t71 = rSges(5,1) * t246 + t162 * rSges(5,3) + t160 * t260;
t190 = Icges(6,5) * t159 + Icges(6,6) * t161;
t64 = -Icges(6,3) * t160 + t162 * t190;
t250 = qJD(1) * t64;
t32 = t159 * t66 - t161 * t68;
t42 = qJD(1) * t65 - qJD(5) * t79;
t44 = -qJD(5) * t81 + (t160 * t192 + t253) * qJD(1);
t285 = qJD(5) * t32 + t159 * t44 + t161 * t42 + t250;
t100 = Icges(6,5) * t161 - Icges(6,6) * t159;
t188 = t102 * t159 - t161 * t104;
t92 = t191 * qJD(5);
t93 = t192 * qJD(5);
t284 = qJD(1) * t100 + qJD(5) * t188 + t159 * t93 + t161 * t92;
t196 = t159 * t65 - t161 * t67;
t63 = Icges(6,3) * t162 + t160 * t190;
t251 = qJD(1) * t63;
t230 = qJD(5) * t160;
t43 = qJD(1) * t66 + t102 * t230;
t45 = qJD(1) * t68 + qJD(5) * t80;
t283 = qJD(5) * t196 - t159 * t45 - t161 * t43 + t251;
t225 = qJD(1) * qJD(5);
t95 = qJDD(5) * t160 + t162 * t225;
t281 = t95 / 0.2e1;
t96 = qJDD(5) * t162 - t160 * t225;
t280 = t96 / 0.2e1;
t279 = -m(5) - m(6);
t278 = t160 / 0.2e1;
t277 = -t162 / 0.2e1;
t275 = pkin(2) * t160;
t274 = pkin(4) * t166;
t259 = rSges(6,2) * t159;
t263 = rSges(6,1) * t161;
t109 = -t259 + t263;
t252 = qJ(4) * t162;
t199 = -t252 - t163;
t226 = qJD(1) * qJD(4);
t236 = qJDD(3) * t160 + t162 * t227;
t173 = qJDD(4) * t162 - 0.2e1 * t160 * t226 + t171 * t199 + t236;
t147 = t162 * qJ(3);
t106 = -t147 + t275;
t200 = -qJ(4) * t160 - t276;
t241 = qJ(4) + t168;
t151 = t160 * rSges(6,3);
t70 = t162 * t202 - t151;
t185 = t160 * t241 + t162 * t274 + t200 + t70;
t179 = -t106 + t185;
t154 = t162 * rSges(6,3);
t84 = t109 * t162;
t46 = -qJD(5) * t84 + (t160 * t202 + t154) * qJD(1);
t145 = qJD(3) * t162;
t82 = qJD(1) * t288 - t145;
t94 = t202 * qJD(5);
t7 = -t94 * t230 + t95 * t109 + t179 * qJDD(1) + (-t82 - t46 + (t252 + t287) * qJD(1)) * qJD(1) + t173;
t273 = t7 * t160;
t158 = qJDD(1) * t163;
t172 = qJDD(1) * t252 + qJDD(4) * t160 + 0.2e1 * t162 * t226 + t171 * t200 + t158 + t286;
t216 = t166 * t231;
t237 = pkin(4) * t216 + t168 * t232;
t69 = rSges(6,1) * t248 + rSges(6,2) * t247 + t154;
t265 = -t162 * t241 + t131 + t69;
t221 = qJD(5) * t263;
t219 = t160 * t221 + t202 * t231;
t220 = qJD(5) * t259;
t47 = (-rSges(6,3) * qJD(1) - t220) * t160 + t219;
t8 = -t96 * t109 + (qJD(5) * t94 - qJDD(3)) * t162 + t265 * qJDD(1) + (qJ(4) * t232 + t237 + t47) * qJD(1) + t172;
t272 = t8 * t162;
t271 = -pkin(2) - qJ(4);
t270 = -pkin(2) + t168;
t257 = rSges(4,3) * t162;
t233 = qJD(4) * t162 + t144;
t87 = t109 * t230;
t19 = qJD(1) * t179 + t233 + t87;
t256 = t19 * t162;
t249 = t100 * t162;
t76 = t160 * t100;
t189 = t161 * t102 + t159 * t104;
t37 = t162 * t189 - t76;
t243 = t37 * qJD(1);
t242 = t190 * qJD(1);
t238 = rSges(5,1) * t216 + t231 * t260;
t234 = rSges(4,2) * t232 + rSges(4,3) * t231;
t229 = qJD(5) * t162;
t228 = -m(4) + t279;
t224 = qJDD(3) * t162;
t223 = -rSges(5,3) + t271;
t21 = t162 * t63 + t65 * t247 + t67 * t248;
t22 = -t162 * t64 - t66 * t247 - t68 * t248;
t98 = qJD(1) * t106;
t222 = -t98 + t233;
t218 = t137 + t233;
t215 = -t230 / 0.2e1;
t214 = t230 / 0.2e1;
t213 = -t229 / 0.2e1;
t212 = t229 / 0.2e1;
t211 = rSges(4,2) * t160 + t257 - t276;
t210 = t147 - t276;
t209 = qJD(4) * t160 - t145;
t208 = -t171 * t276 + t158;
t197 = t159 * t67 + t161 * t65;
t175 = qJD(1) * t197 + qJD(5) * t76 + t250;
t176 = -qJD(1) * t195 - qJD(5) * t249 + t251;
t206 = (t175 * t160 + t162 * t283) * t162 + t160 * (t176 * t160 - t162 * t285);
t205 = t160 * (t160 * t285 + t176 * t162) + t162 * (-t160 * t283 + t175 * t162);
t130 = rSges(2,1) * t170 - rSges(2,2) * t169;
t129 = rSges(2,1) * t169 + rSges(2,2) * t170;
t203 = rSges(5,1) * t166 + t260;
t194 = t160 * t22 + t162 * t21;
t58 = t160 * t63;
t23 = -t197 * t162 + t58;
t24 = -t160 * t64 + t290;
t193 = t160 * t24 + t162 * t23;
t72 = -t160 * rSges(5,3) + t162 * t203;
t187 = t200 + t72;
t186 = t288 - t199;
t184 = -t106 + t187;
t183 = t202 + t274;
t174 = t189 * qJD(1) - t190 * qJD(5);
t83 = t109 * t160;
t36 = t160 * t189 + t249;
t35 = t36 * qJD(1);
t34 = (t186 + t71) * qJD(1) + t209;
t33 = qJD(1) * t184 + t233;
t30 = qJD(2) + (-t160 * t69 - t162 * t70) * qJD(5);
t26 = qJD(1) * t234 + qJDD(1) * t111 + t208 - t224 + t286;
t25 = -t295 + (-t106 + t211) * qJDD(1) + (-qJD(1) * t111 - t82) * qJD(1) + t236;
t20 = -t109 * t229 + (t186 + t265) * qJD(1) + t209;
t15 = -t224 + qJDD(1) * t71 + qJD(1) * (-rSges(5,3) * t232 + t238) + t172;
t14 = t184 * qJDD(1) + (-qJD(1) * t71 - t82) * qJD(1) + t173;
t13 = -t160 * t284 + t174 * t162;
t12 = t174 * t160 + t162 * t284;
t11 = t195 * qJD(5) - t159 * t42 + t161 * t44;
t10 = -qJD(5) * t197 - t159 * t43 + t161 * t45;
t9 = -t69 * t95 - t70 * t96 + qJDD(2) + (-t160 * t47 + t162 * t46) * qJD(5);
t6 = qJD(5) * t193 - t243;
t5 = qJD(5) * t194 + t35;
t1 = [-m(2) * (-g(1) * t129 + g(2) * t130) - t95 * t37 / 0.2e1 + t32 * t281 + (-qJD(5) * t189 + t159 * t92 - t161 * t93) * qJD(1) + (t35 + ((-t23 + t58 + t22) * t160 + (t24 - t290 + (-t197 + t64) * t160 + t21) * t162) * qJD(5)) * t215 + (t36 - t196) * t280 + (t6 + t243 + (t160 ^ 2 * t64 + (-t58 + t22 + (t197 + t64) * t162) * t162) * qJD(5)) * t213 + (t10 + t13) * t212 + ((-t108 * t171 - g(2) + t208) * t90 + (-g(1) - t295 + (-0.2e1 * t112 - t163 + t90) * t171) * t89) * m(3) + (t11 + t12 + t5) * t214 + (t19 * (-t209 + (-t220 + t221) * t162) + t20 * (-t160 * t220 + t218 + t219 + t237) + ((-t20 * t169 - t19 * t170) * pkin(1) + (-rSges(6,3) + t270) * t256 + (t19 * (-qJ(3) - t183) + t20 * (-rSges(6,3) - pkin(2))) * t160) * qJD(1) - (t185 * qJD(1) - t19 + t222 + t87) * t20 + (-g(2) + t8) * (t217 + t69 - t287) + (-g(1) + t7) * (t160 * t270 + t162 * t183 - t151 + t210)) * m(6) + (-(qJD(1) * t187 + t222 - t33) * t34 - t33 * t209 + t34 * (t218 + t238) + ((-t34 * t169 - t33 * t170) * pkin(1) + t33 * t223 * t162 + (t33 * (-qJ(3) - t203) + t34 * t223) * t160) * qJD(1) + (-g(2) + t15) * (t217 + t71 + t252) + (-g(1) + t14) * (t160 * t271 + t210 + t72)) * m(5) + ((-g(1) + t25) * (t257 + (rSges(4,2) - pkin(2)) * t160 + t210) + (-g(2) + t26) * t62 + (-t144 + t98 + t234 + t235 + (-t211 - t275 - t276) * qJD(1)) * (qJD(1) * t62 - t145)) * m(4) + (Icges(5,1) * t167 ^ 2 + (-0.2e1 * Icges(5,4) * t167 + Icges(5,2) * t166) * t166 - t188 + m(2) * (t129 ^ 2 + t130 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,1) + m(3) * (t112 * t90 + t89 ^ 2)) * qJDD(1); m(6) * t9 + (m(3) + m(4) + m(5)) * qJDD(2) + (-m(3) + t228) * g(3); t228 * (g(1) * t160 - g(2) * t162) + 0.2e1 * (t273 / 0.2e1 - t272 / 0.2e1) * m(6) + 0.2e1 * (t14 * t278 + t15 * t277) * m(5) + 0.2e1 * (t25 * t278 + t26 * t277) * m(4); t279 * (g(1) * t162 + g(2) * t160) + m(5) * (t14 * t162 + t15 * t160) + m(6) * (t160 * t8 + t162 * t7); -t5 * t232 / 0.2e1 + t162 * (qJD(1) * t13 + t205 * qJD(5) + qJDD(1) * t36 + t21 * t96 + t22 * t95) / 0.2e1 + t194 * t280 + ((-t21 * t160 + t22 * t162) * qJD(1) + t205) * t212 + t6 * t231 / 0.2e1 + (qJD(1) * t12 + qJD(5) * t206 - qJDD(1) * t37 + t23 * t96 + t24 * t95) * t278 + t193 * t281 + ((-t23 * t160 + t24 * t162) * qJD(1) + t206) * t214 + qJDD(1) * (t160 * t32 - t162 * t196) / 0.2e1 + qJD(1) * (t10 * t162 + t11 * t160 + (t160 * t196 + t32 * t162) * qJD(1)) / 0.2e1 + ((t76 * t229 - t242) * t162 + (-t293 + (-t162 * t249 - t294) * qJD(5)) * t160) * t213 + ((-t230 * t249 - t242) * t160 + (t293 + (t160 * t76 + t294) * qJD(5)) * t162) * t215 - qJD(1) * ((-t240 * t159 - t239 * t161) * qJD(1) + (t159 * t177 + t161 * t178) * qJD(5)) / 0.2e1 + ((t20 * t94 - t9 * t70 + t30 * (-qJD(1) * t69 + t46)) * t162 + (-t19 * t94 - t9 * t69 + t30 * (qJD(1) * t70 - t47)) * t160 + (t273 - t272 + (t20 * t160 + t256) * qJD(1)) * t109 - (t19 * t84 + t20 * t83) * qJD(1) - (t30 * (-t160 * t83 - t162 * t84) - (t19 * t160 - t162 * t20) * t202) * qJD(5) - g(1) * t83 + g(2) * t84 + g(3) * t202) * m(6);];
tau = t1;
