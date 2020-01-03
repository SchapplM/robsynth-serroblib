% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:21
% EndTime: 2020-01-03 11:57:28
% DurationCPUTime: 4.81s
% Computational Cost: add. (11071->391), mult. (9267->518), div. (0->0), fcn. (8721->10), ass. (0->212)
t193 = sin(pkin(9));
t194 = cos(pkin(9));
t195 = sin(qJ(5));
t197 = cos(qJ(5));
t192 = qJ(1) + qJ(2);
t184 = pkin(8) + t192;
t178 = sin(t184);
t179 = cos(t184);
t279 = t194 * t195;
t113 = -t178 * t197 + t179 * t279;
t278 = t194 * t197;
t114 = t178 * t195 + t179 * t278;
t283 = t179 * t193;
t66 = Icges(6,5) * t114 - Icges(6,6) * t113 + Icges(6,3) * t283;
t106 = Icges(6,4) * t114;
t70 = Icges(6,2) * t113 - Icges(6,6) * t283 - t106;
t105 = Icges(6,4) * t113;
t72 = Icges(6,1) * t114 + Icges(6,5) * t283 - t105;
t31 = -(t195 * t70 + t197 * t72) * t193 + t194 * t66;
t233 = t113 * t70 + t114 * t72;
t27 = t283 * t66 + t233;
t323 = t179 * t27;
t130 = t178 * rSges(4,1) + t179 * rSges(4,2);
t185 = sin(t192);
t180 = pkin(2) * t185;
t309 = t180 + t130;
t191 = qJD(1) + qJD(2);
t186 = cos(t192);
t310 = t185 * rSges(3,1) + t186 * rSges(3,2);
t122 = t310 * t191;
t196 = sin(qJ(1));
t296 = pkin(1) * qJD(1);
t264 = t196 * t296;
t109 = t264 + t122;
t111 = -t178 * t279 - t179 * t197;
t112 = t178 * t278 - t179 * t195;
t322 = -t111 * t70 + t112 * t72;
t167 = -qJD(5) * t194 + t191;
t117 = -Icges(6,3) * t194 + (Icges(6,5) * t197 - Icges(6,6) * t195) * t193;
t288 = Icges(6,4) * t197;
t118 = -Icges(6,6) * t194 + (-Icges(6,2) * t195 + t288) * t193;
t289 = Icges(6,4) * t195;
t119 = -Icges(6,5) * t194 + (Icges(6,1) * t197 - t289) * t193;
t211 = -t113 * t118 + t114 * t119 + t117 * t283;
t318 = t211 * t167;
t131 = t179 * pkin(3) + t178 * qJ(4);
t281 = t186 * t191;
t162 = pkin(2) * t281;
t317 = -t191 * t131 + t162;
t316 = t178 * pkin(3) + t180;
t124 = -rSges(6,3) * t194 + (rSges(6,1) * t197 - rSges(6,2) * t195) * t193;
t266 = qJD(5) * t193;
t257 = t124 * t266;
t221 = (-qJD(4) - t257) * t179;
t75 = rSges(6,1) * t114 - rSges(6,2) * t113 + rSges(6,3) * t283;
t315 = t167 * t75 + t162 + t221;
t282 = t179 * t194;
t116 = pkin(4) * t282 + pkin(7) * t283;
t198 = cos(qJ(1));
t183 = t198 * t296;
t34 = t183 + (t116 + t131) * t191 + t315;
t292 = t191 * t34;
t304 = pkin(4) * t194;
t314 = (-t180 + (-t304 - pkin(3) + (-rSges(6,3) - pkin(7)) * t193) * t178) * t292;
t256 = t178 * t266;
t286 = t178 * t193;
t74 = t112 * rSges(6,1) + t111 * rSges(6,2) + rSges(6,3) * t286;
t313 = -t124 * t256 + t167 * t74;
t312 = t191 * t309;
t251 = -qJ(4) * t179 + t316;
t285 = t178 * t194;
t242 = pkin(4) * t285 + pkin(7) * t286 + t251;
t297 = rSges(4,2) * t178;
t132 = t179 * rSges(4,1) - t297;
t121 = t191 * t132;
t284 = t179 * t191;
t150 = rSges(4,1) * t284;
t311 = -t121 + t150;
t267 = qJD(4) * t179;
t287 = t178 * t191;
t272 = pkin(3) * t284 + qJ(4) * t287;
t236 = -t267 + t272;
t260 = t191 * t282;
t280 = t191 * t193;
t261 = t179 * t280;
t273 = pkin(4) * t260 + pkin(7) * t261;
t81 = -qJD(5) * t112 - t113 * t191;
t82 = qJD(5) * t111 + t114 * t191;
t48 = t82 * rSges(6,1) + t81 * rSges(6,2) + rSges(6,3) * t261;
t308 = -t191 * t116 + t236 + t273 - t315 + t317 + t48;
t270 = rSges(5,1) * t282 + t178 * rSges(5,3);
t102 = -rSges(5,2) * t283 + t270;
t249 = t162 - t267;
t274 = rSges(5,1) * t260 + rSges(5,3) * t287;
t307 = -t191 * t102 - t249 + t272 + t274 + t317;
t104 = Icges(6,4) * t111;
t71 = Icges(6,1) * t112 + Icges(6,5) * t286 + t104;
t213 = t178 * (-Icges(6,2) * t112 + t104 + t71) - t179 * (Icges(6,2) * t114 + t105 - t72);
t290 = Icges(6,4) * t112;
t68 = Icges(6,2) * t111 + Icges(6,6) * t286 + t290;
t306 = t178 * (-Icges(6,1) * t111 + t290 + t68) - t179 * (-Icges(6,1) * t113 - t106 + t70);
t190 = t191 ^ 2;
t305 = pkin(2) * t190;
t188 = t196 * pkin(1);
t189 = t198 * pkin(1);
t303 = -t113 * t68 + t114 * t71;
t298 = rSges(3,2) * t185;
t294 = t178 * t66;
t65 = Icges(6,5) * t112 + Icges(6,6) * t111 + Icges(6,3) * t286;
t293 = t179 * t65;
t61 = t183 + (t102 + t131) * t191 + t249;
t291 = t191 * t61;
t135 = (-Icges(6,1) * t195 - t288) * t193;
t277 = t118 - t135;
t134 = (-Icges(6,2) * t197 - t289) * t193;
t276 = t119 + t134;
t262 = t178 * t280;
t275 = rSges(5,2) * t262 + rSges(5,3) * t284;
t147 = qJ(4) * t284;
t168 = qJD(4) * t178;
t271 = t147 + t168;
t199 = qJD(1) ^ 2;
t182 = t199 * t189;
t269 = t186 * t305 + t182;
t24 = t111 * t68 + t112 * t71 + t65 * t286;
t25 = -t286 * t66 - t322;
t265 = t199 * t188;
t153 = rSges(5,1) * t285;
t263 = t191 * t236 + t269;
t181 = pkin(2) * t186;
t259 = t181 + t131;
t253 = -t266 / 0.2e1;
t252 = t266 / 0.2e1;
t143 = rSges(3,1) * t186 - t298;
t110 = t143 * t191 + t183;
t248 = t178 * t253;
t247 = t178 * t252;
t246 = t179 * t253;
t245 = t179 * t252;
t244 = t191 * t252;
t241 = -rSges(5,2) * t286 + t153;
t243 = -rSges(5,3) * t179 + t241 + t251;
t123 = rSges(3,1) * t281 - t191 * t298;
t240 = -t168 + t264;
t239 = t132 + t181;
t79 = qJD(5) * t114 + t111 * t191;
t80 = qJD(5) * t113 + t112 * t191;
t238 = rSges(6,1) * t80 + rSges(6,2) * t79;
t237 = -rSges(5,2) * t280 - qJD(4);
t33 = t191 * t242 + t240 + t313;
t232 = -t178 * t33 - t179 * t34;
t231 = -t178 * t75 + t179 * t74;
t230 = t178 * (Icges(6,5) * t111 - Icges(6,6) * t112) - t179 * (Icges(6,5) * t113 + Icges(6,6) * t114);
t228 = t259 + t270;
t227 = t178 * t244;
t226 = t179 * t244;
t222 = (-rSges(6,1) * t195 - rSges(6,2) * t197) * t193;
t220 = (t178 * t24 - t179 * t25) * t193;
t26 = -t283 * t65 - t303;
t219 = (t178 * t26 - t323) * t193;
t218 = t241 + t316;
t217 = -t185 * t305 - t265;
t216 = -t238 + t271;
t133 = (-Icges(6,5) * t195 - Icges(6,6) * t197) * t193;
t215 = t191 * t168 + t217;
t210 = t74 + t242;
t207 = t75 + t259 + t116;
t10 = qJD(5) * t219 - t318;
t42 = Icges(6,5) * t82 + Icges(6,6) * t81 + Icges(6,3) * t261;
t44 = Icges(6,4) * t82 + Icges(6,2) * t81 + Icges(6,6) * t261;
t46 = Icges(6,1) * t82 + Icges(6,4) * t81 + Icges(6,5) * t261;
t13 = -t194 * t42 + (-t195 * t44 + t197 * t46 + (-t195 * t71 - t197 * t68) * qJD(5)) * t193;
t41 = Icges(6,5) * t80 + Icges(6,6) * t79 + Icges(6,3) * t262;
t43 = Icges(6,4) * t80 + Icges(6,2) * t79 + Icges(6,6) * t262;
t45 = Icges(6,1) * t80 + Icges(6,4) * t79 + Icges(6,5) * t262;
t14 = -t194 * t41 + (-t195 * t43 + t197 * t45 + (t195 * t72 - t197 * t70) * qJD(5)) * t193;
t125 = qJD(5) * t133;
t126 = qJD(5) * t134;
t127 = qJD(5) * t135;
t19 = t113 * t126 - t114 * t127 + t118 * t79 + t119 * t80 + (t117 * t287 - t125 * t179) * t193;
t20 = t111 * t126 + t112 * t127 + t118 * t81 + t119 * t82 + (t117 * t284 + t125 * t178) * t193;
t30 = -t194 * t65 + (-t195 * t68 + t197 * t71) * t193;
t39 = -t125 * t194 + (-t126 * t195 + t127 * t197 + (-t118 * t197 - t119 * t195) * qJD(5)) * t193;
t36 = t39 * t167;
t49 = t111 * t118 + t112 * t119 + t117 * t286;
t40 = t49 * t167;
t9 = qJD(5) * t220 + t40;
t206 = t36 + (t40 + ((-t233 + t24 + t27) * t178 + (t26 + (t293 - t294) * t193 - t25 + t303) * t179) * t266) * t245 + (t13 + t20) * t247 + (t31 - t211) * t227 + (t30 + t49) * t226 + (t10 + t318 + (t323 + (t303 + t25 + (t293 + t294) * t193 + t322) * t178) * t266) * t248 + (t9 + t14 + t19) * t246;
t205 = ((t191 * t24 - t111 * t43 - t112 * t45 - t70 * t81 + t72 * t82 - (t178 * t41 - t284 * t66) * t193) * t179 + (t191 * t25 + t111 * t44 + t112 * t46 + t68 * t81 + t71 * t82 + (t178 * t42 + t284 * t65) * t193) * t178) * t193;
t204 = ((t191 * t26 - t113 * t43 + t114 * t45 - t70 * t79 + t72 * t80 - (-t179 * t41 - t287 * t66) * t193) * t179 + (t191 * t27 + t113 * t44 - t114 * t46 + t68 * t79 + t71 * t80 + (-t179 * t42 + t287 * t65) * t193) * t178) * t193;
t203 = ((t191 * t30 - t14) * t179 + (t191 * t31 + t13) * t178) * t193;
t93 = t264 + t312;
t94 = t183 + t162 + t121;
t201 = (-t93 * t297 - t94 * t309) * t191;
t97 = pkin(3) * t287 - t271;
t51 = (-t153 * t191 + t275 - t97) * t191 + t215;
t52 = (t179 * t237 + t274) * t191 + t263;
t60 = t191 * t243 + t240;
t200 = (-t51 * rSges(5,2) * t193 + t52 * (-rSges(5,3) - qJ(4)) + t60 * t237) * t179 + (-t180 + (-rSges(5,1) * t194 - pkin(3)) * t178) * t291;
t128 = qJD(5) * t222;
t100 = t123 * t191 + t182;
t99 = -t122 * t191 - t265;
t92 = t191 * (-rSges(4,2) * t287 + t150) + t269;
t91 = -t130 * t190 + t217;
t90 = rSges(6,1) * t113 + rSges(6,2) * t114;
t89 = rSges(6,1) * t111 - rSges(6,2) * t112;
t47 = rSges(6,3) * t262 + t238;
t35 = t231 * t266 + qJD(3);
t22 = -t128 * t256 + t167 * t48 + (t221 + t273) * t191 + t263;
t21 = -t128 * t179 * t266 - t167 * t47 + (-t97 + (t257 + (-pkin(7) * t193 - t304) * t191) * t178) * t191 + t215;
t16 = ((-t191 * t75 + t48) * t179 + (-t191 * t74 + t47) * t178) * t266;
t1 = [m(3) * (t99 * (t143 + t189) + t100 * (t188 + t310) + (-t110 + t123 + t183) * t109) + t206 + (t21 * (t189 + t207) + t34 * (t216 - t264) + t22 * (t188 + t210) + (t34 + t308) * t33 + t314) * m(6) + (t51 * (t189 + t228) + t61 * (t147 - t240 + t275) + t52 * (t188 + t218) + t200 + (t61 + t307) * t60) * m(5) + (t91 * (t189 + t239) - t94 * t264 + t92 * (t188 + t309) + t201 + (t94 + t311) * t93) * m(4); t206 + (t21 * t207 + t22 * t210 + t242 * t292 + (-t168 + t216 + t313) * t34 + t308 * t33 + t314) * m(6) + (t52 * t218 + t51 * t228 + t243 * t291 + t200 + (-t168 + t271 + t275) * t61 + t307 * t60) * m(5) + (t91 * t239 + t309 * t92 + t311 * t93 + t94 * t312 + t201) * m(4) + (t100 * t310 + t109 * t123 - t110 * t122 + t143 * t99 - (t109 * t143 - t110 * t310) * t191) * m(3); m(6) * t16; m(5) * (-t178 * t52 - t179 * t51) + m(6) * (-t22 * t178 - t21 * t179); -t194 * (qJD(5) * t203 + t36) / 0.2e1 + t167 * (-t39 * t194 + t203) / 0.2e1 + (qJD(5) * t205 + t167 * t20) * t286 / 0.2e1 + (-t194 * t49 + t220) * t226 + (-t194 * t20 + t205) * t247 - (qJD(5) * t204 + t167 * t19) * t283 / 0.2e1 + (t194 * t211 + t219) * t227 + (-t19 * t194 + t204) * t246 - t167 * (-t194 * t133 * t167 + ((-t195 * t276 - t197 * t277) * t167 + ((-t213 * t195 - t197 * t306) * t193 - t230 * t194) * qJD(5)) * t193) / 0.2e1 + ((t111 * t276 - t112 * t277 + t133 * t286) * t167 + (t111 * t213 - t112 * t306 + t230 * t286) * t266) * t248 + ((t113 * t276 + t114 * t277 - t133 * t283) * t167 + (t213 * t113 + t306 * t114 - t230 * t283) * t266) * t245 + (t178 * t10 + t179 * t9) * t280 / 0.2e1 + ((-t21 * t75 - t22 * t74 - t33 * t48 + t34 * t47) * t194 + (t16 * t231 + t35 * (t178 * t47 + t179 * t48 - t284 * t75 - t287 * t74) + t232 * t128 + ((-t191 * t33 - t21) * t179 + (-t22 + t292) * t178) * t124) * t193 - (t33 * t89 - t34 * t90) * t167 - (t35 * (t178 * t90 + t179 * t89) + t232 * t222) * t266) * m(6);];
tauc = t1(:);
