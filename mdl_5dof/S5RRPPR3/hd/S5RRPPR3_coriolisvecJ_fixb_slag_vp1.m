% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:28
% EndTime: 2019-12-31 19:26:33
% DurationCPUTime: 3.26s
% Computational Cost: add. (6758->326), mult. (5034->409), div. (0->0), fcn. (3666->8), ass. (0->199)
t162 = qJ(1) + qJ(2);
t156 = pkin(8) + t162;
t153 = sin(t156);
t158 = cos(t162);
t155 = pkin(2) * t158;
t154 = cos(t156);
t292 = t154 * pkin(3) + t153 * qJ(4);
t229 = t155 + t292;
t303 = -t153 * rSges(5,3) - t229;
t161 = qJD(1) + qJD(2);
t163 = sin(qJ(5));
t165 = cos(qJ(5));
t258 = Icges(6,4) * t163;
t203 = Icges(6,2) * t165 + t258;
t75 = Icges(6,6) * t154 + t153 * t203;
t248 = t153 * t165;
t123 = Icges(6,4) * t248;
t249 = t153 * t163;
t77 = Icges(6,1) * t249 + Icges(6,5) * t154 + t123;
t207 = t163 * t77 + t165 * t75;
t76 = -Icges(6,6) * t153 + t154 * t203;
t257 = Icges(6,4) * t165;
t204 = Icges(6,1) * t163 + t257;
t78 = -Icges(6,5) * t153 + t154 * t204;
t42 = t163 * t76 - t165 * t78;
t188 = t203 * t161;
t129 = -Icges(6,2) * t163 + t257;
t287 = -Icges(6,6) * t161 + qJD(5) * t129;
t49 = t287 * t153 + t154 * t188;
t189 = t204 * t161;
t131 = Icges(6,1) * t165 - t258;
t286 = -Icges(6,5) * t161 + qJD(5) * t131;
t51 = t286 * t153 + t154 * t189;
t302 = -qJD(5) * t207 + t161 * t42 - t163 * t49 + t165 * t51;
t242 = t129 + t204;
t243 = -t203 + t131;
t301 = (t163 * t242 - t165 * t243) * t161;
t300 = -t154 * rSges(4,1) - t155;
t164 = sin(qJ(1));
t262 = pkin(1) * qJD(1);
t233 = t164 * t262;
t157 = sin(t162);
t107 = rSges(3,1) * t157 + rSges(3,2) * t158;
t253 = t107 * t161;
t81 = -t233 - t253;
t211 = rSges(6,1) * t163 + rSges(6,2) * t165;
t299 = 0.2e1 * qJD(5);
t241 = qJD(4) * t161;
t247 = t154 * t161;
t118 = qJ(4) * t247;
t137 = qJD(4) * t153;
t245 = t118 + t137;
t250 = t153 * t161;
t298 = t153 * t241 + t161 * (-pkin(3) * t250 + t245);
t202 = Icges(6,5) * t163 + Icges(6,6) * t165;
t297 = t202 * t161;
t205 = t163 * t78 + t165 * t76;
t296 = t205 * t154;
t294 = -rSges(5,2) * t154 - t303;
t295 = t294 * t161;
t293 = -rSges(4,2) * t153 - t300;
t100 = t153 * rSges(5,2) + t154 * rSges(5,3);
t244 = rSges(5,2) * t250 + rSges(5,3) * t247;
t140 = t154 * qJ(4);
t271 = pkin(3) * t153;
t99 = -t140 + t271;
t91 = t161 * t99;
t259 = t137 - t91;
t291 = -t161 * t100 + t244 - t259;
t150 = t154 * pkin(7);
t213 = -t150 - t155;
t147 = t154 * rSges(6,3);
t79 = rSges(6,1) * t249 + rSges(6,2) * t248 + t147;
t264 = rSges(6,2) * t163;
t265 = rSges(6,1) * t165;
t135 = -t264 + t265;
t239 = qJD(5) * t154;
t97 = t135 * t239;
t290 = (t292 - t213 + t79) * t161 - t97;
t48 = t153 * t188 - t287 * t154;
t50 = t153 * t189 - t286 * t154;
t74 = -Icges(6,3) * t153 + t154 * t202;
t289 = qJD(5) * t42 + t161 * t74 + t163 * t50 + t165 * t48;
t127 = Icges(6,5) * t165 - Icges(6,6) * t163;
t288 = -Icges(6,3) * t161 + qJD(5) * t127;
t110 = t203 * qJD(5);
t111 = t204 * qJD(5);
t285 = qJD(5) * (t129 * t163 - t131 * t165) + t110 * t165 + t111 * t163 + t127 * t161;
t206 = t163 * t75 - t165 * t77;
t73 = Icges(6,3) * t154 + t153 * t202;
t284 = qJD(5) * t206 + t161 * t73 - t163 * t51 - t165 * t49;
t273 = pkin(2) * t157;
t214 = -pkin(7) * t153 - t273;
t232 = qJD(5) * t265;
t230 = t153 * t232 + t211 * t247;
t231 = qJD(5) * t264;
t80 = -t153 * rSges(6,3) + t154 * t211;
t69 = t161 * t80;
t283 = -t153 * t231 - t214 * t161 + t230 + t245 - t69;
t267 = t129 * t154 + t78;
t269 = -t131 * t154 + t76;
t281 = t163 * t269 - t165 * t267;
t268 = -Icges(6,2) * t249 + t123 + t77;
t270 = -t131 * t153 + t75;
t280 = t163 * t270 - t165 * t268;
t160 = t161 ^ 2;
t279 = -pkin(3) - pkin(7);
t278 = t153 / 0.2e1;
t277 = -t154 / 0.2e1;
t275 = -t161 / 0.2e1;
t274 = pkin(1) * t164;
t272 = pkin(2) * t160;
t166 = cos(qJ(1));
t159 = t166 * pkin(1);
t152 = t158 * rSges(3,1);
t201 = t165 * t129 + t163 * t131;
t83 = t153 * t127;
t56 = t154 * t201 - t83;
t260 = t56 * t161;
t251 = t127 * t154;
t246 = t157 * t161;
t240 = qJD(5) * t153;
t238 = -rSges(6,3) + t279;
t26 = t154 * t73 + t75 * t248 + t77 * t249;
t27 = -t154 * t74 - t76 * t248 - t78 * t249;
t167 = qJD(1) ^ 2;
t237 = t167 * t274;
t236 = t167 * t159;
t234 = t166 * t262;
t226 = -t240 / 0.2e1;
t224 = -t239 / 0.2e1;
t101 = rSges(4,1) * t153 + rSges(4,2) * t154;
t184 = -t101 - t273;
t108 = -rSges(3,2) * t157 + t152;
t95 = -rSges(3,2) * t246 + t152 * t161;
t218 = t137 - t233;
t138 = qJD(4) * t154;
t217 = -t138 + t234;
t215 = -t271 - t273;
t96 = t135 * t240;
t199 = t218 + t96;
t30 = (t214 + t80 - t99) * t161 + t199;
t31 = t217 + t290;
t209 = t153 * t30 - t154 * t31;
t208 = -t153 * t79 - t154 * t80;
t198 = (t231 - t232) * t154;
t194 = (t153 * t27 + t154 * t26) * qJD(5);
t63 = t153 * t73;
t28 = -t154 * t207 + t63;
t29 = -t153 * t74 + t296;
t193 = (t153 * t29 + t154 * t28) * qJD(5);
t192 = -t157 * t272 - t237;
t191 = -t158 * t272 - t236;
t190 = t138 - t198;
t186 = -pkin(2) * t246 - t233;
t185 = t150 + t229 + t79;
t181 = t140 + t215 + t100;
t179 = t153 * t297 - t288 * t154 - t161 * t205;
t178 = t288 * t153 + t154 * t297 + t161 * t207;
t177 = -t202 * qJD(5) + t161 * t201;
t52 = (t153 * t211 + t147) * t161 + t198;
t53 = (-rSges(6,3) * t161 - t231) * t153 + t230;
t175 = (-t161 * t79 + t52) * t154 + (-t53 + t69) * t153;
t10 = t193 - t260;
t15 = qJD(5) * t205 - t163 * t48 + t165 * t50;
t20 = t177 * t153 + t285 * t154;
t21 = -t285 * t153 + t177 * t154;
t55 = t153 * t201 + t251;
t54 = t55 * t161;
t9 = t54 + t194;
t174 = (t54 + ((-t28 + t63 + t27) * t153 + (t29 - t296 + (-t207 + t74) * t153 + t26) * t154) * qJD(5)) * t226 + (-qJD(5) * t201 + t110 * t163 - t111 * t165) * t161 + (t260 + (t153 ^ 2 * t74 + (-t63 + t27 + (t207 + t74) * t154) * t154) * qJD(5) + t10) * t224 + (t15 + t20 + t9) * t240 / 0.2e1 + (t21 + t302) * t239 / 0.2e1 + (t154 * t56 + (-t206 + t55) * t153) * qJD(5) * t275;
t171 = t153 * t279 + t140 - t273 + t80;
t66 = t161 * t184 - t233;
t67 = t161 * t293 + t234;
t170 = (t67 * t184 + t66 * t300) * t161;
t44 = (t100 - t99 - t273) * t161 + t218;
t45 = t217 + t295;
t169 = (t45 * t215 + t44 * t303) * t161;
t168 = ((-t157 * t31 - t158 * t30) * pkin(2) + t30 * t238 * t154 + (t30 * (-qJ(4) - t211) + t31 * t238) * t153) * t161;
t122 = rSges(5,2) * t247;
t120 = rSges(4,2) * t250;
t117 = t154 * t241;
t112 = t211 * qJD(5);
t93 = t161 * t101;
t90 = t135 * t154;
t89 = t135 * t153;
t82 = t108 * t161 + t234;
t72 = -t161 * t95 - t236;
t71 = -t161 * t253 - t237;
t70 = t161 * t292 - t138;
t58 = -t161 * (rSges(4,1) * t247 - t120) + t191;
t57 = -t101 * t160 + t192;
t36 = qJD(5) * t208 + qJD(3);
t33 = t117 + (-rSges(5,3) * t250 + t122 - t70) * t161 + t191;
t32 = t161 * t244 + t192 + t298;
t17 = -t236 - t112 * t240 + t117 + t213 * t160 + (-t52 - t70 + t97) * t161;
t16 = -t237 + t161 * t53 + t214 * t160 + (t112 * t154 + t135 * t250) * qJD(5) + t298;
t11 = t175 * qJD(5);
t1 = [m(3) * (t72 * (-t107 - t274) + t71 * (t108 + t159) + (-t95 - t234 + t82) * t81) + t174 + (t17 * (t171 - t274) + t30 * (t190 - t234) + t16 * (t159 + t185) + t168 + (-t199 - t233 + t283 + t30 + t91) * t31) * m(6) + (t33 * (t181 - t274) + t44 * (t122 - t217) + t32 * (t159 + t294) + t169 + (-t186 + t44 + t118 + t218 + t291) * t45) * m(5) + (t58 * (t184 - t274) + t66 * (t120 - t234) + t57 * (t159 + t293) + t170 + (-t186 + t66 + t93 - t233) * t67) * m(4); t174 + (t16 * t185 + t17 * t171 + t168 + (-t259 + t283 - t96) * t31 + (-t138 + t190 + t290) * t30) * m(6) + (t33 * t181 + t32 * t294 + t169 + (t122 + t295) * t44 + (t273 * t161 + t245 + t291) * t45) * m(5) + (t67 * t93 - (-t273 * t67 - t293 * t66) * t161 + t66 * t120 + t58 * t184 + t57 * t293 + t170) * m(4) + (-(-t82 * t107 - t108 * t81) * t161 - t107 * t72 + t108 * t71 - t81 * t95 - t82 * t253) * m(3); m(6) * t11; 0.2e1 * (t16 * t277 + t17 * t278) * m(6) + 0.2e1 * (t277 * t32 + t278 * t33) * m(5); t161 * (t302 * t154 + (t161 * t206 + t15) * t153) / 0.2e1 + ((t83 * t239 - t297) * t154 + (-t301 + (t281 * t153 + (-t280 - t251) * t154) * qJD(5)) * t153) * t224 + ((-t240 * t251 - t297) * t153 + (t301 + (t280 * t154 + (-t281 + t83) * t153) * qJD(5)) * t154) * t226 + ((-t163 * t243 - t165 * t242) * t161 + ((t153 * t269 - t154 * t270) * t165 + (t267 * t153 - t268 * t154) * t163) * qJD(5)) * t275 + (t161 * t20 + ((t178 * t153 + t284 * t154 + t161 * t29) * t154 + (t179 * t153 - t289 * t154 - t161 * t28) * t153) * t299) * t278 + (t161 * t21 + ((-t284 * t153 + t178 * t154 + t161 * t27) * t154 + (t289 * t153 + t179 * t154 - t161 * t26) * t153) * t299) * t154 / 0.2e1 - (t9 + t194) * t250 / 0.2e1 + (t10 + t193) * t247 / 0.2e1 + (t11 * t208 + t36 * t175 - t209 * t112 + ((t161 * t30 - t16) * t154 + (t161 * t31 + t17) * t153) * t135 - (t30 * t90 + t31 * t89) * t161 - (t36 * (-t153 * t89 - t154 * t90) - t209 * t211) * qJD(5)) * m(6);];
tauc = t1(:);
