% Calculate vector of inverse dynamics joint torques for
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:52
% EndTime: 2019-12-31 17:43:59
% DurationCPUTime: 6.68s
% Computational Cost: add. (6297->432), mult. (8519->546), div. (0->0), fcn. (8084->8), ass. (0->196)
t189 = cos(qJ(1));
t183 = t189 * pkin(1);
t184 = qJ(1) + pkin(7);
t181 = sin(t184);
t182 = cos(t184);
t139 = rSges(3,1) * t181 + rSges(3,2) * t182;
t188 = sin(qJ(1));
t277 = pkin(1) * t188;
t127 = -t139 - t277;
t190 = qJD(1) ^ 2;
t299 = t190 * t183;
t186 = cos(pkin(8));
t258 = t181 * t186;
t275 = pkin(6) * t182;
t129 = pkin(4) * t258 + t275;
t185 = sin(pkin(8));
t187 = sin(qJ(5));
t278 = cos(qJ(5));
t144 = t185 * t278 - t186 * t187;
t117 = t144 * t181;
t200 = t185 * t187 + t186 * t278;
t118 = t200 * t181;
t70 = -t118 * rSges(6,1) - t117 * rSges(6,2) - rSges(6,3) * t182;
t298 = -t129 + t70;
t119 = t144 * t182;
t120 = t200 * t182;
t109 = Icges(6,4) * t118;
t63 = Icges(6,2) * t117 + Icges(6,6) * t182 + t109;
t108 = Icges(6,4) * t117;
t67 = -Icges(6,1) * t118 - Icges(6,5) * t182 - t108;
t273 = t119 * t63 - t120 * t67;
t60 = Icges(6,5) * t118 + Icges(6,6) * t117 + Icges(6,3) * t182;
t16 = -t181 * t60 + t273;
t256 = t182 * t186;
t257 = t182 * t185;
t105 = rSges(5,1) * t256 + t181 * rSges(5,2) + rSges(5,3) * t257;
t126 = pkin(3) * t256 + qJ(4) * t257;
t140 = t182 * pkin(2) + t181 * qJ(3);
t290 = t183 + t140;
t294 = t126 + t290;
t56 = t294 + t105;
t297 = t117 * t63 - t118 * t67;
t28 = -t144 * t67 - t200 * t63;
t14 = t182 * t60 + t297;
t263 = Icges(6,4) * t120;
t65 = Icges(6,2) * t119 - Icges(6,6) * t181 + t263;
t110 = Icges(6,4) * t119;
t68 = Icges(6,1) * t120 - Icges(6,5) * t181 + t110;
t272 = t119 * t65 + t120 * t68;
t62 = Icges(6,5) * t120 + Icges(6,6) * t119 - Icges(6,3) * t181;
t223 = t181 * t62 - t272;
t295 = t223 - t14;
t106 = rSges(4,1) * t256 - rSges(4,2) * t257 + t181 * rSges(4,3);
t87 = t106 + t290;
t293 = t117 * t65 + t118 * t68;
t141 = t182 * rSges(3,1) - rSges(3,2) * t181;
t128 = t141 + t183;
t160 = pkin(4) * t256;
t291 = pkin(6) * t181 - t160;
t289 = qJD(5) * t144;
t286 = t181 * (-Icges(6,2) * t120 + t110 + t68) - t182 * (-Icges(6,2) * t118 + t108 - t67);
t284 = -m(5) - m(6);
t240 = qJD(1) * qJD(5);
t131 = -qJDD(5) * t181 - t182 * t240;
t283 = t131 / 0.2e1;
t132 = qJDD(5) * t182 - t181 * t240;
t282 = t132 / 0.2e1;
t281 = t181 / 0.2e1;
t280 = -t182 / 0.2e1;
t279 = -rSges(6,3) - pkin(6);
t274 = g(1) * t181;
t134 = t200 * qJD(5);
t247 = qJD(1) * t181;
t74 = -t134 * t182 - t144 * t247;
t75 = t182 * t289 - t200 * t247;
t269 = t75 * rSges(6,1) + t74 * rSges(6,2);
t262 = Icges(6,4) * t144;
t98 = -Icges(6,2) * t200 + t262;
t268 = -Icges(6,1) * t200 - t262 - t98;
t137 = Icges(6,4) * t200;
t100 = Icges(6,1) * t144 - t137;
t265 = -Icges(6,2) * t144 + t100 - t137;
t253 = t120 * rSges(6,1) + t119 * rSges(6,2);
t71 = -rSges(6,3) * t181 + t253;
t264 = t71 - t291;
t261 = qJ(4) * t185;
t96 = Icges(6,5) * t144 - Icges(6,6) * t200;
t26 = t100 * t118 + t117 * t98 + t182 * t96;
t260 = qJD(1) * t26;
t259 = t181 * t185;
t154 = rSges(4,2) * t259;
t246 = qJD(1) * t182;
t252 = rSges(4,3) * t246 + qJD(1) * t154;
t245 = qJD(4) * t185;
t152 = t182 * t245;
t169 = qJD(3) * t181;
t251 = t152 + t169;
t250 = t182 * rSges(4,3) + t154;
t249 = qJ(3) * t246 + t169;
t172 = t182 * qJ(3);
t138 = pkin(2) * t181 - t172;
t248 = -qJD(1) * t138 + t169;
t244 = qJD(5) * t181;
t243 = qJD(5) * t182;
t242 = -m(4) + t284;
t241 = qJD(1) * qJD(3);
t239 = qJDD(4) * t185;
t15 = t182 * t62 + t293;
t238 = rSges(4,1) * t258;
t236 = t152 + t249;
t233 = t181 * t245;
t232 = -rSges(4,1) * t186 - pkin(2);
t231 = -t244 / 0.2e1;
t230 = t244 / 0.2e1;
t229 = -t243 / 0.2e1;
t228 = t243 / 0.2e1;
t227 = -t138 - t277;
t225 = t172 - t277;
t212 = pkin(3) * t186 + t261;
t125 = t212 * t181;
t224 = -qJD(1) * t125 + t152 + t248;
t167 = -qJDD(4) * t186 + qJDD(2);
t222 = qJDD(1) * t183 - t190 * t277;
t104 = t238 - t250;
t221 = -t104 + t227;
t220 = -t125 + t227;
t76 = qJD(1) * t119 - t134 * t181;
t77 = qJD(1) * t120 + t181 * t289;
t218 = -rSges(6,1) * t77 - rSges(6,2) * t76;
t33 = Icges(6,5) * t75 + Icges(6,6) * t74 - Icges(6,3) * t246;
t34 = Icges(6,5) * t77 + Icges(6,6) * t76 - Icges(6,3) * t247;
t35 = Icges(6,4) * t75 + Icges(6,2) * t74 - Icges(6,6) * t246;
t36 = Icges(6,4) * t77 + Icges(6,2) * t76 - Icges(6,6) * t247;
t37 = Icges(6,1) * t75 + Icges(6,4) * t74 - Icges(6,5) * t246;
t38 = Icges(6,1) * t77 + Icges(6,4) * t76 - Icges(6,5) * t247;
t217 = (t119 * t36 + t120 * t38 - t181 * t34 - t246 * t60 + t63 * t74 - t67 * t75) * t182 - t181 * (t119 * t35 + t120 * t37 - t181 * t33 - t246 * t62 + t65 * t74 + t68 * t75);
t216 = -t181 * (t117 * t35 + t118 * t37 + t182 * t33 - t247 * t62 + t65 * t76 + t68 * t77) + t182 * (t117 * t36 + t118 * t38 + t182 * t34 - t247 * t60 + t63 * t76 - t67 * t77);
t170 = qJD(3) * t182;
t215 = t170 - t233;
t159 = rSges(2,1) * t189 - rSges(2,2) * t188;
t158 = rSges(2,1) * t188 + rSges(2,2) * t189;
t213 = rSges(5,1) * t186 + rSges(5,3) * t185;
t211 = t14 * t182 - t15 * t181;
t210 = t16 * t182 + t181 * t223;
t209 = t181 * (Icges(6,5) * t119 - Icges(6,6) * t120) - t182 * (Icges(6,5) * t117 - Icges(6,6) * t118);
t205 = qJDD(3) * t181 + t182 * t241 - t299;
t176 = t182 * rSges(5,2);
t103 = t181 * t213 - t176;
t204 = -t103 + t220;
t203 = -pkin(2) - t212;
t116 = qJD(1) * t140 - t170;
t202 = -t212 * t246 - t116 - 0.2e1 * t233;
t201 = t220 + t298;
t199 = t182 * t239 + t205;
t198 = -(Icges(6,1) * t119 - t263 - t65) * t181 + (Icges(6,1) * t117 - t109 - t63) * t182;
t196 = -pkin(4) * t186 + t203;
t193 = -pkin(2) + (-rSges(5,1) - pkin(3)) * t186 + (-rSges(5,3) - qJ(4)) * t185;
t192 = -qJDD(3) * t182 + qJD(1) * (-pkin(2) * t247 + t249) + qJDD(1) * t140 + t181 * t241 + t222;
t191 = qJDD(1) * t126 + t181 * t239 + t192 + (-t212 * t247 + 0.2e1 * t152) * qJD(1);
t166 = rSges(5,2) * t246;
t102 = rSges(6,1) * t144 - rSges(6,2) * t200;
t101 = -rSges(6,1) * t200 - rSges(6,2) * t144;
t95 = -Icges(6,5) * t200 - Icges(6,6) * t144;
t92 = -rSges(6,1) * t134 - rSges(6,2) * t289;
t91 = -Icges(6,1) * t134 - Icges(6,4) * t289;
t90 = -Icges(6,4) * t134 - Icges(6,2) * t289;
t89 = -Icges(6,5) * t134 - Icges(6,6) * t289;
t88 = t102 * t243;
t85 = rSges(6,1) * t119 - rSges(6,2) * t120;
t84 = rSges(6,1) * t117 - rSges(6,2) * t118;
t59 = qJD(1) * t87 - t170;
t58 = qJD(1) * t221 + t169;
t46 = qJD(1) * t56 - t215;
t45 = qJD(1) * t204 + t251;
t40 = -rSges(6,3) * t247 - t218;
t39 = -rSges(6,3) * t246 + t269;
t32 = qJDD(1) * t106 + qJD(1) * (-qJD(1) * t238 + t252) + t192;
t31 = t221 * qJDD(1) + (-qJD(1) * t106 - t116) * qJD(1) + t205;
t30 = -qJD(4) * t186 + qJD(2) + (t181 * t70 - t182 * t71) * qJD(5);
t29 = t144 * t68 - t200 * t65;
t27 = t100 * t120 + t119 * t98 - t181 * t96;
t25 = t27 * qJD(1);
t24 = -t170 + (qJD(5) * t102 + t245) * t181 + (t294 + t264) * qJD(1);
t23 = qJD(1) * t201 + t251 + t88;
t19 = qJDD(1) * t105 + (-t213 * t247 + t166) * qJD(1) + t191;
t18 = t204 * qJDD(1) + (-qJD(1) * t105 + t202) * qJD(1) + t199;
t13 = t100 * t77 + t117 * t90 + t118 * t91 + t182 * t89 - t247 * t96 + t76 * t98;
t12 = t100 * t75 + t119 * t90 + t120 * t91 - t181 * t89 - t246 * t96 + t74 * t98;
t11 = -t131 * t70 - t132 * t71 + (-t181 * t40 - t182 * t39) * qJD(5) + t167;
t10 = -t134 * t68 + t144 * t37 - t200 * t35 - t289 * t65;
t9 = t134 * t67 + t144 * t38 - t200 * t36 - t289 * t63;
t8 = t264 * qJDD(1) + (-qJD(1) * t129 + t39) * qJD(1) + t191 + t92 * t244 - t131 * t102;
t7 = t92 * t243 + t132 * t102 + t201 * qJDD(1) + (qJD(1) * t291 + t202 - t40) * qJD(1) + t199;
t6 = qJD(5) * t210 + t25;
t5 = qJD(5) * t211 + t260;
t1 = [(-t100 * t134 + t144 * t91 - t200 * t90 - t289 * t98) * qJD(1) + (t25 + (t273 * t182 + (t295 + t297) * t181) * qJD(5)) * t229 - m(2) * (-g(1) * t158 + g(2) * t159) + (t29 + t27) * t283 + (t28 + t26) * t282 + (t10 + t12) * t231 + (t5 - t260 + ((t272 + t295) * t182 + t293 * t181) * qJD(5)) * t230 + ((-t139 * t190 - g(2) + t222) * t128 + (-t299 + (-0.2e1 * t141 - t183 + t128) * t190 - g(1)) * t127) * m(3) + (t13 + t9 + t6) * t228 + (-(-t261 - pkin(2) + (-pkin(3) - pkin(4)) * t186) * t274 + t23 * (t170 + t218) + t24 * (t236 + t269) + (t7 * t196 - t23 * t245) * t181 + ((-t24 * t188 - t23 * t189) * pkin(1) + (t196 * t23 + t24 * t279) * t182 + (t23 * (-qJ(3) - t279) + t24 * t196) * t181) * qJD(1) - (-t23 + t88 + t224 + (-t277 + t298) * qJD(1)) * t24 + (-g(2) + t8) * (t181 * t279 + t160 + t253 + t294) + (-g(1) + t7) * (t225 - t275 + t70)) * m(6) + (t45 * t215 + t46 * (t166 + t236) + ((-t46 * t188 - t45 * t189) * pkin(1) + t45 * t193 * t182 + (t45 * (-rSges(5,2) - qJ(3)) + t46 * (t203 - t213)) * t181) * qJD(1) - (-t45 + (-t103 - t277) * qJD(1) + t224) * t46 + (-g(2) + t19) * t56 + (-g(1) + t18) * (t181 * t193 + t176 + t225)) * m(5) + (t58 * t170 + t59 * (t249 + t252) + ((-t59 * t188 - t58 * t189) * pkin(1) + t58 * (rSges(4,2) * t185 + t232) * t182 + (t58 * (-rSges(4,3) - qJ(3)) + t59 * t232) * t181) * qJD(1) - (-t58 + (-t104 - t277) * qJD(1) + t248) * t59 + (t32 - g(2)) * t87 + (t31 - g(1)) * (t232 * t181 + t225 + t250)) * m(4) + (m(2) * (t158 ^ 2 + t159 ^ 2) + m(3) * (t127 ^ 2 + t141 * t128) + t100 * t144 - t200 * t98 + Icges(2,3) + Icges(3,3) + (Icges(5,3) + Icges(4,2)) * t186 ^ 2 + ((Icges(4,1) + Icges(5,1)) * t185 + (2 * Icges(4,4) - 2 * Icges(5,5)) * t186) * t185) * qJDD(1); (m(3) + m(4)) * qJDD(2) + m(5) * t167 + m(6) * t11 + (-m(3) + t242) * g(3); t242 * (-g(2) * t182 + t274) + 0.2e1 * (t280 * t8 + t281 * t7) * m(6) + 0.2e1 * (t18 * t281 + t19 * t280) * m(5) + 0.2e1 * (t280 * t32 + t281 * t31) * m(4); t284 * (-g(3) * t186 + (g(1) * t182 + g(2) * t181) * t185) + m(5) * (-t167 * t186 + t18 * t257 + t19 * t259) + m(6) * (-t11 * t186 + t257 * t7 + t259 * t8); -t6 * t246 / 0.2e1 - t181 * (t12 * qJD(1) + t217 * qJD(5) + t27 * qJDD(1) - t131 * t223 + t16 * t132) / 0.2e1 + t210 * t283 + ((-t16 * t181 + t182 * t223) * qJD(1) + t217) * t231 - t5 * t247 / 0.2e1 + t182 * (t13 * qJD(1) + t216 * qJD(5) + t26 * qJDD(1) + t15 * t131 + t14 * t132) / 0.2e1 + t211 * t282 + ((-t14 * t181 - t15 * t182) * qJD(1) + t216) * t228 + qJDD(1) * (-t29 * t181 + t28 * t182) / 0.2e1 + qJD(1) * (-t10 * t181 + t9 * t182 + (-t28 * t181 - t182 * t29) * qJD(1)) / 0.2e1 + ((t119 * t265 + t120 * t268 - t181 * t95) * qJD(1) + (-t119 * t286 + t120 * t198 + t181 * t209) * qJD(5)) * t230 + ((t117 * t265 + t118 * t268 + t182 * t95) * qJD(1) + (-t117 * t286 + t198 * t118 - t209 * t182) * qJD(5)) * t229 - qJD(1) * ((t268 * t144 - t200 * t265) * qJD(1) + (t144 * t198 + t200 * t286) * qJD(5)) / 0.2e1 + ((t23 * t92 - t11 * t71 + t30 * (qJD(1) * t70 - t39)) * t182 + (t24 * t92 + t11 * t70 + t30 * (qJD(1) * t71 - t40)) * t181 + (t181 * t8 + t7 * t182 + (-t23 * t181 + t24 * t182) * qJD(1)) * t102 - (-t23 * t84 + t24 * t85) * qJD(1) - (t30 * (-t181 * t84 - t182 * t85) + (t24 * t181 + t23 * t182) * t101) * qJD(5) - g(1) * t85 - g(2) * t84 - g(3) * t101) * m(6);];
tau = t1;
