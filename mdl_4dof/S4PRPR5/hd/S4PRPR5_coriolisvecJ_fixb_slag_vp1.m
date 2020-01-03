% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:58
% EndTime: 2019-12-31 16:23:11
% DurationCPUTime: 10.09s
% Computational Cost: add. (6804->440), mult. (10500->703), div. (0->0), fcn. (10000->8), ass. (0->243)
t178 = qJ(2) + pkin(7);
t174 = sin(t178);
t175 = cos(t178);
t183 = sin(qJ(2));
t185 = cos(qJ(2));
t352 = -Icges(3,5) * t183 - Icges(4,5) * t174 - Icges(3,6) * t185 - Icges(4,6) * t175;
t353 = t352 * qJD(2);
t308 = Icges(3,4) * t183;
t236 = -Icges(3,2) * t185 - t308;
t307 = Icges(3,4) * t185;
t241 = -Icges(3,1) * t183 - t307;
t285 = qJD(2) * t175;
t286 = qJD(2) * t174;
t305 = Icges(4,4) * t175;
t306 = Icges(4,4) * t174;
t351 = -(-t183 * t236 + t185 * t241) * qJD(2) + (-Icges(4,2) * t175 - t306) * t286 - (-Icges(4,1) * t174 - t305) * t285;
t179 = sin(pkin(6));
t182 = sin(qJ(4));
t184 = cos(qJ(4));
t254 = rSges(5,1) * t184 - rSges(5,2) * t182;
t96 = -rSges(5,3) * t175 + t174 * t254;
t350 = t179 * t96;
t180 = cos(pkin(6));
t349 = t180 * t96;
t176 = t179 ^ 2;
t177 = t180 ^ 2;
t287 = t176 + t177;
t348 = t352 * t179;
t347 = t353 * t179;
t346 = t353 * t180;
t345 = t352 * t180;
t237 = -Icges(3,2) * t183 + t307;
t117 = Icges(3,6) * t179 + t180 * t237;
t242 = Icges(3,1) * t185 - t308;
t119 = Icges(3,5) * t179 + t180 * t242;
t235 = -Icges(4,2) * t174 + t305;
t240 = Icges(4,1) * t175 - t306;
t339 = -(Icges(4,6) * t179 + t180 * t235) * t175 - (Icges(4,5) * t179 + t180 * t240) * t174;
t344 = t351 * t180 + (t117 * t185 + t119 * t183 - t339) * qJD(2);
t116 = -Icges(3,6) * t180 + t179 * t237;
t118 = -Icges(3,5) * t180 + t179 * t242;
t340 = (-Icges(4,6) * t180 + t179 * t235) * t175 + (-Icges(4,5) * t180 + t179 * t240) * t174;
t343 = t351 * t179 + (t116 * t185 + t118 * t183 + t340) * qJD(2);
t279 = qJD(4) * t174;
t284 = qJD(2) * t179;
t160 = t180 * t279 + t284;
t328 = t160 / 0.2e1;
t282 = qJD(2) * t180;
t161 = t179 * t279 - t282;
t326 = t161 / 0.2e1;
t338 = qJD(4) / 0.2e1;
t169 = rSges(3,1) * t183 + rSges(3,2) * t185;
t213 = qJD(2) * t169;
t230 = Icges(5,5) * t184 - Icges(5,6) * t182;
t90 = -Icges(5,3) * t175 + t174 * t230;
t301 = Icges(5,4) * t184;
t233 = -Icges(5,2) * t182 + t301;
t92 = -Icges(5,6) * t175 + t174 * t233;
t302 = Icges(5,4) * t182;
t238 = Icges(5,1) * t184 - t302;
t94 = -Icges(5,5) * t175 + t174 * t238;
t293 = t180 * t182;
t294 = t179 * t184;
t148 = -t175 * t293 + t294;
t296 = t174 * t180;
t292 = t180 * t184;
t295 = t179 * t182;
t149 = t175 * t292 + t295;
t303 = Icges(5,4) * t149;
t56 = Icges(5,2) * t148 + Icges(5,6) * t296 + t303;
t129 = Icges(5,4) * t148;
t58 = Icges(5,1) * t149 + Icges(5,5) * t296 + t129;
t246 = -t182 * t56 + t184 * t58;
t146 = -t175 * t295 - t292;
t297 = t174 * t179;
t147 = t175 * t294 - t293;
t304 = Icges(5,4) * t147;
t55 = Icges(5,2) * t146 + Icges(5,6) * t297 + t304;
t128 = Icges(5,4) * t146;
t57 = Icges(5,1) * t147 + Icges(5,5) * t297 + t128;
t247 = -t182 * t55 + t184 * t57;
t332 = -(-t180 * t90 - t246) * t160 - (-t179 * t90 - t247) * t161;
t331 = -t183 * (t236 * t179 + t118) - t185 * (-t241 * t179 + t116);
t133 = (-Icges(5,2) * t184 - t302) * t174;
t278 = qJD(4) * t175;
t192 = t160 * (-Icges(5,2) * t149 + t129 + t58) + t161 * (-Icges(5,2) * t147 + t128 + t57) - t278 * (t133 + t94);
t186 = qJD(2) ^ 2;
t269 = t175 * t282;
t272 = t182 * t286;
t80 = -qJD(4) * t149 + t180 * t272;
t271 = t184 * t286;
t81 = qJD(4) * t148 - t180 * t271;
t41 = Icges(5,5) * t81 + Icges(5,6) * t80 + Icges(5,3) * t269;
t54 = Icges(5,5) * t149 + Icges(5,6) * t148 + Icges(5,3) * t296;
t218 = t174 * t41 + t285 * t54;
t43 = Icges(5,4) * t81 + Icges(5,2) * t80 + Icges(5,6) * t269;
t45 = Icges(5,1) * t81 + Icges(5,4) * t80 + Icges(5,5) * t269;
t10 = t148 * t43 + t149 * t45 + t180 * t218 + t56 * t80 + t58 * t81;
t132 = (-Icges(5,5) * t182 - Icges(5,6) * t184) * t174;
t91 = Icges(5,3) * t174 + t175 * t230;
t49 = qJD(2) * t91 + qJD(4) * t132;
t217 = t174 * t49 + t285 * t90;
t22 = t148 * t56 + t149 * t58 + t296 * t54;
t53 = Icges(5,5) * t147 + Icges(5,6) * t146 + Icges(5,3) * t297;
t21 = t148 * t55 + t149 * t57 + t296 * t53;
t315 = t179 * t21;
t250 = t180 * t22 + t315;
t36 = t148 * t92 + t149 * t94 + t296 * t90;
t310 = t36 * t174;
t93 = Icges(5,6) * t174 + t175 * t233;
t50 = qJD(2) * t93 + qJD(4) * t133;
t134 = (-Icges(5,1) * t182 - t301) * t174;
t95 = Icges(5,5) * t174 + t175 * t238;
t51 = qJD(2) * t95 + qJD(4) * t134;
t190 = -(t148 * t50 + t149 * t51 + t180 * t217 + t80 * t92 + t81 * t94) * t175 + (t175 * t250 + t310) * qJD(2);
t270 = t175 * t284;
t78 = -qJD(4) * t147 + t179 * t272;
t79 = qJD(4) * t146 - t179 * t271;
t40 = Icges(5,5) * t79 + Icges(5,6) * t78 + Icges(5,3) * t270;
t219 = t174 * t40 + t285 * t53;
t42 = Icges(5,4) * t79 + Icges(5,2) * t78 + Icges(5,6) * t270;
t44 = Icges(5,1) * t79 + Icges(5,4) * t78 + Icges(5,5) * t270;
t9 = t148 * t42 + t149 * t44 + t180 * t219 + t55 * t80 + t57 * t81;
t330 = t10 * t328 + t190 * t338 + t326 * t9;
t329 = -t160 / 0.2e1;
t327 = -t161 / 0.2e1;
t325 = pkin(2) * t183;
t324 = pkin(2) * t185;
t98 = -qJ(3) * t180 + t179 * t324;
t99 = qJ(3) * t179 + t180 * t324;
t319 = t179 * t98 + t180 * t99;
t317 = pkin(2) * qJD(2);
t316 = t175 * t90;
t20 = t146 * t56 + t147 * t58 + t297 * t54;
t313 = t180 * t20;
t35 = t146 * t92 + t147 * t94 + t297 * t90;
t311 = t35 * t174;
t276 = t183 * t317;
t280 = qJD(3) * t180;
t162 = -t179 * t276 - t280;
t173 = qJD(3) * t179;
t163 = -t180 * t276 + t173;
t289 = t162 * t284 + t163 * t282;
t288 = t179 * t162 + t180 * t163;
t277 = t186 * t324;
t275 = t185 * t317;
t274 = t99 * t282 + t98 * t284 + qJD(1);
t267 = t282 / 0.2e1;
t266 = -t278 / 0.2e1;
t265 = t278 / 0.2e1;
t164 = rSges(4,1) * t174 + rSges(4,2) * t175;
t264 = -t164 - t325;
t165 = rSges(4,1) * t175 - rSges(4,2) * t174;
t263 = -t165 - t324;
t166 = pkin(3) * t174 - pkin(5) * t175;
t262 = -t166 - t325;
t167 = pkin(3) * t175 + pkin(5) * t174;
t261 = -t167 - t324;
t260 = qJD(2) * t338;
t259 = t262 - t96;
t258 = t287 * t325;
t257 = t174 * t260;
t256 = t175 * t260;
t152 = t165 * qJD(2);
t255 = -t152 - t275;
t170 = rSges(3,1) * t185 - rSges(3,2) * t183;
t153 = t167 * qJD(2);
t47 = rSges(5,1) * t81 + rSges(5,2) * t80 + rSges(5,3) * t269;
t137 = (-rSges(5,1) * t182 - rSges(5,2) * t184) * t174;
t203 = rSges(5,3) * t174 + t175 * t254;
t52 = qJD(2) * t203 + qJD(4) * t137;
t60 = rSges(5,1) * t149 + rSges(5,2) * t148 + rSges(5,3) * t296;
t16 = -t179 * t277 - t47 * t278 - t160 * t52 + (-t153 * t179 + (t174 * t60 - t175 * t349) * qJD(4)) * qJD(2);
t46 = rSges(5,1) * t79 + rSges(5,2) * t78 + rSges(5,3) * t270;
t59 = rSges(5,1) * t147 + rSges(5,2) * t146 + rSges(5,3) * t297;
t17 = -t180 * t277 + t46 * t278 + t161 * t52 + (-t153 * t180 + (-t174 * t59 + t175 * t350) * qJD(4)) * qJD(2);
t253 = -t16 * t180 + t17 * t179;
t252 = t54 * t160 + t53 * t161;
t19 = t146 * t55 + t147 * t57 + t297 * t53;
t251 = t179 * t19 + t313;
t23 = t174 * t247 - t175 * t53;
t24 = t174 * t246 - t175 * t54;
t249 = t23 * t179 + t24 * t180;
t248 = -t179 * t60 + t180 * t59;
t245 = -t182 * t92 + t184 * t94;
t244 = qJD(2) * t264;
t243 = qJD(2) * t262;
t229 = t287 * t170;
t228 = t287 * t213;
t227 = -t153 - t52 - t275;
t226 = t179 * t256;
t225 = t180 * t256;
t222 = -t245 + t91;
t130 = t167 * t179;
t131 = t167 * t180;
t18 = t160 * t59 - t161 * t60 + (t130 * t179 + t131 * t180) * qJD(2) + t274;
t221 = t18 * t166;
t220 = -qJD(2) * t152 - t277;
t216 = qJD(2) * t166;
t215 = t18 * t248;
t104 = -rSges(4,3) * t180 + t165 * t179;
t105 = rSges(4,3) * t179 + t165 * t180;
t37 = (t104 * t179 + t105 * t180) * qJD(2) + t274;
t214 = t37 * t164;
t212 = qJD(2) * t164;
t205 = (t241 * t180 - t117) * t185 + (-t236 * t180 - t119) * t183;
t204 = t132 * t278 - t160 * (Icges(5,5) * t148 - Icges(5,6) * t149) - t161 * (Icges(5,5) * t146 - Icges(5,6) * t147);
t198 = t174 * t204;
t193 = (Icges(5,1) * t148 - t303 - t56) * t160 + (Icges(5,1) * t146 - t304 - t55) * t161 - (t134 - t92) * t278;
t191 = -(t146 * t50 + t147 * t51 + t179 * t217 + t78 * t92 + t79 * t94) * t175 + (t175 * t251 + t311) * qJD(2);
t38 = t174 * t245 - t316;
t189 = -((qJD(2) * t245 - t49) * t175 + (qJD(2) * t90 - t182 * t50 + t184 * t51 + (-t182 * t94 - t184 * t92) * qJD(4)) * t174) * t175 + (t38 * t174 + t175 * t249) * qJD(2);
t188 = t339 * t179 + t340 * t180;
t187 = (-t222 * t278 - t332) * t174;
t115 = t180 * t216;
t114 = t179 * t216;
t113 = t180 * t212;
t112 = t179 * t212;
t85 = t220 * t180;
t84 = t220 * t179;
t83 = t180 * t244 + t173;
t82 = t179 * t244 - t280;
t75 = t94 * t180;
t74 = t94 * t179;
t73 = t92 * t180;
t72 = t92 * t179;
t69 = rSges(5,1) * t148 - rSges(5,2) * t149;
t68 = rSges(5,1) * t146 - rSges(5,2) * t147;
t61 = t228 * qJD(2);
t48 = qJD(2) * t229 + qJD(1);
t39 = (-t112 * t179 - t113 * t180) * qJD(2) + t289;
t30 = t161 * t96 + t180 * t243 + t278 * t59 + t173;
t29 = -t160 * t96 + t179 * t243 - t278 * t60 - t280;
t12 = t160 * t46 - t161 * t47 + (-t114 * t179 - t115 * t180 + t248 * t278) * qJD(2) + t289;
t11 = t160 * t24 + t161 * t23 - t278 * t38;
t8 = t146 * t43 + t147 * t45 + t179 * t218 + t56 * t78 + t58 * t79;
t7 = t146 * t42 + t147 * t44 + t179 * t219 + t55 * t78 + t57 * t79;
t6 = t160 * t22 + t161 * t21 - t278 * t36;
t5 = t160 * t20 + t161 * t19 - t278 * t35;
t4 = (qJD(2) * t246 - t41) * t175 + (qJD(2) * t54 - t182 * t43 + t184 * t45 + (-t182 * t58 - t184 * t56) * qJD(4)) * t174;
t3 = (qJD(2) * t247 - t40) * t175 + (qJD(2) * t53 - t182 * t42 + t184 * t44 + (-t182 * t57 - t184 * t55) * qJD(4)) * t174;
t1 = qJD(4) * t191 + t160 * t8 + t161 * t7;
t2 = [-m(3) * t61 + m(4) * t39 + m(5) * t12; (t179 * t22 - t180 * t21) * t225 + (t179 * t20 - t180 * t19) * t226 - t11 * t279 / 0.2e1 + (((t182 * t73 - t184 * t75 + t54) * t160 + (t182 * t72 - t184 * t74 + t53) * t161 + t38 * qJD(4)) * t174 + ((t222 * t175 + (t182 * t93 - t184 * t95 - t90) * t174 + t249) * qJD(4) + t332) * t175) * t265 + (t179 * t8 - t180 * t7) * t326 + ((-t146 * t73 - t147 * t75) * t160 + (-t146 * t72 - t147 * t74) * t161 + (t311 + (-t146 * t93 - t147 * t95 + t313) * t175) * qJD(4) + (((t19 - t316) * qJD(4) + t252) * t175 + t187) * t179) * t327 + (t10 * t179 - t180 * t9) * t328 + ((-t148 * t73 - t149 * t75) * t160 + (-t148 * t72 - t149 * t74) * t161 + (t310 + (-t148 * t93 - t149 * t95 + t315) * t175) * qJD(4) + (((t22 - t316) * qJD(4) + t252) * t175 + t187) * t180) * t329 + t179 * t330 - t180 * t1 / 0.2e1 + (t179 * t24 - t180 * t23) * t257 + (t343 * t177 + (t346 * t179 + (-t344 - t347) * t180) * t179) * t284 + (-t347 * t177 + (t344 * t179 + (-t343 + t346) * t180) * t179) * t282 - (t345 * qJD(2) * t176 + (-t331 * t180 + t188 + (t205 - t348) * t179) * t282) * t284 / 0.2e1 + ((t205 * t179 + t188 + (-t331 - t345) * t180) * t284 + t348 * qJD(2) * t177) * t267 + (-(-t160 * t29 + t161 * t30) * t203 - ((t29 * t60 - t30 * t59) * t174 + t215 * t175) * qJD(4) - ((-t180 * t221 + t261 * t30) * t180 + (-t179 * t221 + t261 * t29) * t179) * qJD(2) + t12 * t319 + (t17 * t259 + t30 * t227 + t12 * (t131 + t60)) * t180 + (t16 * t259 + t29 * t227 + t12 * (t130 + t59)) * t179 + (t160 * t350 - t161 * t349 + t258 * qJD(2) + t288 + (-t115 + t47) * t180 + (-t114 + t46) * t179) * t18) * m(5) + (-(-t37 * t258 + (-t180 * t214 + t263 * t83) * t180 + (-t179 * t214 + t263 * t82) * t179) * qJD(2) + t39 * t319 + t37 * t288 + (t39 * t105 - t37 * t113 + t255 * t83 + t264 * t85) * t180 + (t39 * t104 - t37 * t112 + t255 * t82 + t264 * t84) * t179) * m(4) + ((-t3 + t6) * t180 + (t4 + t5) * t179) * t266 + (-t48 * t228 - t61 * t229 + (t169 * t170 * t186 + t213 * t48) * t287) * m(3); m(4) * (t179 * t85 - t180 * t84) + m(5) * t253; t175 * t6 * t267 + t296 * t330 + (t174 * t250 - t175 * t36) * t225 + ((t10 * t180 + t179 * t9) * t174 + t190) * t328 + t5 * t270 / 0.2e1 + t1 * t297 / 0.2e1 + (t174 * t251 - t175 * t35) * t226 + ((t179 * t7 + t180 * t8) * t174 + t191) * t326 + t11 * t286 / 0.2e1 - t175 * (qJD(4) * t189 + t160 * t4 + t161 * t3) / 0.2e1 + (t174 * t249 - t175 * t38) * t257 + ((t179 * t3 + t180 * t4) * t174 + t189) * t266 + (t148 * t192 + t149 * t193 - t180 * t198) * t329 + (t146 * t192 + t147 * t193 - t179 * t198) * t327 + (t204 * t175 + (-t182 * t192 + t193 * t184) * t174) * t265 + ((-t16 * t60 + t17 * t59 - t29 * t47 + t30 * t46 + (t215 + (t179 * t30 - t180 * t29) * t96) * qJD(2)) * t175 + (t30 * (-qJD(2) * t59 + t179 * t52) + t29 * (qJD(2) * t60 - t180 * t52) + t12 * t248 + t18 * (-t179 * t47 + t180 * t46) + t253 * t96) * t174 - t30 * (t137 * t161 + t278 * t68) - t29 * (-t137 * t160 - t278 * t69) - t18 * (t160 * t68 - t161 * t69)) * m(5);];
tauc = t2(:);
