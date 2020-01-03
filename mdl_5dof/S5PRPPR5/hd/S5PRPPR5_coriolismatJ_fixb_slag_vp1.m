% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPPR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:04
% EndTime: 2019-12-31 17:38:13
% DurationCPUTime: 7.04s
% Computational Cost: add. (9527->413), mult. (25480->639), div. (0->0), fcn. (30502->8), ass. (0->221)
t219 = sin(qJ(2));
t221 = cos(qJ(2));
t279 = t219 * t221;
t339 = Icges(4,1) - Icges(4,3);
t338 = -0.2e1 * t219;
t337 = t219 * t339;
t336 = 0.2e1 * (Icges(3,1) - Icges(3,2)) * t279 + (t338 * t219 + 0.2e1 * t221 ^ 2) * Icges(3,4);
t216 = sin(pkin(7));
t241 = -Icges(3,5) * t219 - Icges(3,6) * t221;
t180 = t241 * t216;
t217 = cos(pkin(7));
t181 = t241 * t217;
t334 = Icges(5,1) - Icges(5,2);
t243 = -Icges(4,4) * t219 + Icges(4,6) * t221;
t333 = t243 * t216 + t180;
t332 = t243 * t217 + t181;
t214 = t216 ^ 2;
t215 = t217 ^ 2;
t266 = t214 + t215;
t280 = t217 * t221;
t329 = Icges(4,5) * t338 - t339 * t221;
t331 = t180 - t336 * t217 + (0.2e1 * Icges(4,5) * t280 + Icges(4,6) * t216 - t217 * t337) * t221 + (-Icges(4,4) * t216 + t329 * t217) * t219;
t282 = t216 * t221;
t330 = t181 + t336 * t216 + (-0.2e1 * Icges(4,5) * t282 + Icges(4,6) * t217 + t216 * t337) * t221 + (-Icges(4,4) * t217 - t329 * t216) * t219;
t326 = -2 * Icges(5,4);
t325 = 2 * Icges(5,4);
t257 = m(4) / 0.4e1 + m(5) / 0.4e1 + m(6) / 0.4e1;
t267 = t266 * t279;
t324 = t257 * (t267 - t279);
t323 = t266 * t219;
t294 = sin(pkin(8));
t295 = cos(pkin(8));
t191 = t219 * t294 + t221 * t295;
t192 = t219 * t295 - t221 * t294;
t321 = 2 * qJD(2);
t320 = m(4) / 0.2e1;
t319 = m(5) / 0.2e1;
t318 = m(6) / 0.2e1;
t193 = t219 * pkin(2) - qJ(3) * t221;
t194 = t219 * rSges(4,1) - rSges(4,3) * t221;
t269 = -t193 - t194;
t150 = t269 * t216;
t152 = t269 * t217;
t272 = t150 * t282 + t152 * t280;
t197 = rSges(4,1) * t221 + t219 * rSges(4,3);
t196 = pkin(2) * t221 + t219 * qJ(3);
t270 = t266 * t196;
t80 = t197 * t266 + t270;
t317 = m(4) * (t323 * t80 + t272);
t264 = -pkin(3) * t219 - t193;
t253 = -rSges(5,1) * t192 + rSges(5,2) * t191 + t264;
t110 = t253 * t216;
t112 = t253 * t217;
t275 = t110 * t282 + t112 * t280;
t165 = t191 * t216;
t166 = t192 * t216;
t167 = t191 * t217;
t168 = t192 * t217;
t255 = t270 + (t216 * t282 + t217 * t280) * pkin(3);
t63 = t216 * (rSges(5,1) * t165 + rSges(5,2) * t166) + t217 * (rSges(5,1) * t167 + rSges(5,2) * t168) + t255;
t316 = m(5) * (t323 * t63 + t275);
t218 = sin(qJ(5));
t220 = cos(qJ(5));
t252 = rSges(6,1) * t220 - rSges(6,2) * t218;
t121 = rSges(6,3) * t191 + t252 * t192;
t146 = -t165 * t218 + t217 * t220;
t147 = t165 * t220 + t217 * t218;
t87 = rSges(6,1) * t147 + rSges(6,2) * t146 - rSges(6,3) * t166;
t64 = -t121 * t166 - t191 * t87;
t148 = -t167 * t218 - t216 * t220;
t149 = t167 * t220 - t216 * t218;
t88 = rSges(6,1) * t149 + rSges(6,2) * t148 - rSges(6,3) * t168;
t65 = t121 * t168 + t191 * t88;
t302 = t64 * t280 + t65 * t282;
t96 = -rSges(6,3) * t165 - t252 * t166;
t97 = -rSges(6,3) * t167 - t252 * t168;
t43 = t165 * t88 + t166 * t97 - t167 * t87 - t168 * t96;
t120 = -rSges(6,3) * t192 + t252 * t191;
t49 = -t120 * t166 - t121 * t165 - t191 * t96 + t192 * t87;
t50 = t120 * t168 + t121 * t167 + t191 * t97 - t192 * t88;
t59 = t166 * t88 - t168 * t87;
t315 = m(6) * (-t221 * t43 + (t216 * t50 + t217 * t49 + t59) * t219 + t302);
t235 = -pkin(4) * t192 - pkin(6) * t191 - t121 + t264;
t75 = t235 * t216;
t77 = t235 * t217;
t301 = t77 * t280 + t75 * t282;
t44 = (pkin(4) * t167 - pkin(6) * t168 + t88) * t217 + (pkin(4) * t165 - pkin(6) * t166 + t87) * t216 + t255;
t314 = m(6) * (t323 * t44 + t301);
t313 = m(6) * (t323 * t59 + t302);
t136 = (-rSges(6,1) * t218 - rSges(6,2) * t220) * t192;
t106 = rSges(6,1) * t146 - rSges(6,2) * t147;
t107 = rSges(6,1) * t148 - rSges(6,2) * t149;
t70 = t106 * t216 + t107 * t217;
t312 = m(6) * (-t136 * t323 - t221 * t70);
t311 = -t165 / 0.2e1;
t310 = -t166 / 0.2e1;
t309 = -t167 / 0.2e1;
t308 = -t168 / 0.2e1;
t307 = t191 / 0.2e1;
t306 = -t192 / 0.2e1;
t305 = t216 / 0.2e1;
t304 = -t217 / 0.2e1;
t300 = m(6) * qJD(5);
t140 = Icges(6,4) * t146;
t85 = Icges(6,1) * t147 - Icges(6,5) * t166 + t140;
t299 = -Icges(6,2) * t147 + t140 + t85;
t141 = Icges(6,4) * t148;
t86 = Icges(6,1) * t149 - Icges(6,5) * t168 + t141;
t298 = -Icges(6,2) * t149 + t141 + t86;
t289 = Icges(6,4) * t147;
t83 = Icges(6,2) * t146 - Icges(6,6) * t166 + t289;
t297 = Icges(6,1) * t146 - t289 - t83;
t288 = Icges(6,4) * t149;
t84 = Icges(6,2) * t148 - Icges(6,6) * t168 + t288;
t296 = Icges(6,1) * t148 - t288 - t84;
t287 = Icges(6,4) * t218;
t286 = Icges(6,4) * t220;
t24 = -t216 * t49 + t217 * t50;
t278 = t24 * qJD(4);
t34 = 0.2e1 * (t43 / 0.4e1 - t70 / 0.4e1) * m(6);
t277 = t34 * qJD(1);
t232 = 0.2e1 * t257 * t323;
t256 = t318 + t319 + t320;
t99 = -t256 * t219 + t232;
t276 = t99 * qJD(1);
t242 = -Icges(6,2) * t218 + t286;
t117 = Icges(6,6) * t191 + t242 * t192;
t274 = -t117 + (-Icges(6,1) * t218 - t286) * t192;
t246 = Icges(6,1) * t220 - t287;
t119 = Icges(6,5) * t191 + t246 * t192;
t273 = t119 + (-Icges(6,2) * t220 - t287) * t192;
t271 = t266 * t193;
t268 = -t196 - t197;
t265 = t24 * t318;
t263 = -pkin(3) * t221 - t196;
t254 = -rSges(5,1) * t191 - rSges(5,2) * t192 + t263;
t195 = t219 * rSges(3,1) + rSges(3,2) * t221;
t251 = -t218 * t83 + t220 * t85;
t250 = -t218 * t84 + t220 * t86;
t239 = Icges(6,5) * t220 - Icges(6,6) * t218;
t238 = -t117 * t218 + t119 * t220;
t237 = (-Icges(5,5) * t166 + Icges(5,6) * t165) * t217 - (-Icges(5,5) * t168 + Icges(5,6) * t167) * t216;
t236 = -pkin(4) * t191 + pkin(6) * t192 - t120 + t263;
t234 = -Icges(6,3) * t165 - t239 * t166 + t251;
t233 = -Icges(6,3) * t167 - t239 * t168 + t250;
t100 = Icges(6,5) * t146 - Icges(6,6) * t147;
t36 = -t100 * t166 + t299 * t146 + t297 * t147;
t101 = Icges(6,5) * t148 - Icges(6,6) * t149;
t37 = -t101 * t166 + t298 * t146 + t296 * t147;
t17 = t216 * t37 - t217 * t36;
t38 = -t100 * t168 + t299 * t148 + t297 * t149;
t39 = -t101 * t168 + t298 * t148 + t296 * t149;
t18 = t216 * t39 - t217 * t38;
t231 = t17 * t304 + t18 * t305;
t230 = -pkin(3) * t323 - t271;
t229 = (-Icges(6,3) * t192 + t239 * t191 + t238) * t191;
t228 = (Icges(5,6) * t216 - t167 * t325 + t334 * t168) * t216 - (-Icges(5,6) * t217 - t165 * t325 + t334 * t166) * t217;
t227 = (-Icges(5,5) * t216 + t334 * t167 - t168 * t326) * t216 - (Icges(5,5) * t217 + t334 * t165 - t166 * t326) * t217;
t115 = Icges(6,3) * t191 + t239 * t192;
t116 = -Icges(6,6) * t192 + t242 * t191;
t118 = -Icges(6,5) * t192 + t246 * t191;
t81 = Icges(6,5) * t147 + Icges(6,6) * t146 - Icges(6,3) * t166;
t92 = -Icges(6,6) * t165 - t242 * t166;
t94 = -Icges(6,5) * t165 - t246 * t166;
t25 = t146 * t92 + t147 * t94 - t165 * t81 - t234 * t166;
t82 = Icges(6,5) * t149 + Icges(6,6) * t148 - Icges(6,3) * t168;
t93 = -Icges(6,6) * t167 - t242 * t168;
t95 = -Icges(6,5) * t167 - t246 * t168;
t26 = t146 * t93 + t147 * t95 - t165 * t82 - t233 * t166;
t45 = t146 * t83 + t147 * t85 - t166 * t81;
t46 = t146 * t84 + t147 * t86 - t166 * t82;
t55 = -t115 * t166 + t117 * t146 + t119 * t147;
t3 = -t26 * t168 - t46 * t167 - t45 * t165 + (-t115 * t165 + t116 * t146 + t118 * t147) * t191 - t55 * t192 - (t25 + t229) * t166;
t32 = (-t218 * t92 + t220 * t94 - t81) * t192 + t234 * t191;
t33 = (-t218 * t93 + t220 * t95 - t82) * t192 + t233 * t191;
t27 = t148 * t92 + t149 * t94 - t167 * t81 - t234 * t168;
t28 = t148 * t93 + t149 * t95 - t167 * t82 - t233 * t168;
t47 = t148 * t83 + t149 * t85 - t168 * t81;
t48 = t148 * t84 + t149 * t86 - t168 * t82;
t56 = -t115 * t168 + t117 * t148 + t119 * t149;
t4 = -t48 * t167 - t27 * t166 - t47 * t165 + (-t115 * t167 + t116 * t148 + t118 * t149) * t191 - t56 * t192 - (t28 + t229) * t168;
t52 = t191 * t81 + t251 * t192;
t53 = t191 * t82 + t250 * t192;
t60 = t115 * t191 + t238 * t192;
t222 = (-t166 * t45 - t168 * t46 + t191 * t55) * t311 + (-t166 * t47 - t168 * t48 + t191 * t56) * t309 + (-t166 * t52 - t168 * t53 + t191 * t60) * t306 + t3 * t310 + (-t52 * t165 - t32 * t166 - t53 * t167 - t33 * t168 - t60 * t192 + (t229 + (-t116 * t218 + t118 * t220 - t115) * t192) * t191) * t307 + t4 * t308;
t153 = t268 * t217;
t151 = t268 * t216;
t133 = (-Icges(6,5) * t218 - Icges(6,6) * t220) * t192;
t132 = t266 * t195;
t113 = t254 * t217;
t111 = t254 * t216;
t98 = t232 + (m(4) + m(5) + m(6)) * t219 / 0.2e1;
t89 = -t194 * t266 - t271;
t79 = 0.4e1 * t324;
t78 = t236 * t217;
t76 = t236 * t216;
t72 = t107 * t191 + t136 * t168;
t71 = -t106 * t191 - t136 * t166;
t68 = t216 * (-rSges(5,1) * t166 + rSges(5,2) * t165) + t217 * (-rSges(5,1) * t168 + rSges(5,2) * t167) + t230;
t67 = -t106 * t168 + t107 * t166;
t57 = t312 / 0.2e1;
t54 = (-pkin(4) * t168 - pkin(6) * t167 + t97) * t217 + (-pkin(4) * t166 - pkin(6) * t165 + t96) * t216 + t230;
t41 = t101 * t191 + (-t298 * t218 + t296 * t220) * t192;
t40 = t100 * t191 + (-t299 * t218 + t297 * t220) * t192;
t35 = (t43 + t70) * t318;
t30 = t313 / 0.2e1;
t23 = qJD(2) * t265;
t22 = t44 * t70 + (-t216 * t75 - t217 * t77) * t136;
t16 = t314 + t316 + t317;
t15 = t216 * t28 - t217 * t27;
t14 = t216 * t26 - t217 * t25;
t13 = -t39 * t168 - t38 * t166 + (-t133 * t168 + t273 * t148 + t274 * t149) * t191;
t12 = -t37 * t168 - t36 * t166 + (-t133 * t166 + t273 * t146 + t274 * t147) * t191;
t11 = t43 * t59 + t49 * t64 + t50 * t65;
t9 = t315 / 0.2e1;
t8 = t30 + t9 - t312 / 0.2e1;
t7 = t57 + t30 - t315 / 0.2e1;
t6 = t57 + t9 - t313 / 0.2e1;
t2 = m(6) * t22 + t231;
t1 = m(6) * t11 + t222;
t5 = [0, t98 * qJD(3) + t35 * qJD(5) + (-m(3) * t132 / 0.2e1 + t89 * t320 + t68 * t319 + t54 * t318) * t321, t98 * qJD(2), 0, t35 * qJD(2) + t67 * t300; qJD(3) * t99 - qJD(5) * t34, t16 * qJD(3) + t2 * qJD(5) + (m(5) * (t110 * t111 + t112 * t113 + t63 * t68) + m(4) * (t150 * t151 + t152 * t153 + t80 * t89) + m(6) * (t44 * t54 + t75 * t76 + t77 * t78) + m(3) * (-t132 + t195) * t266 * (rSges(3,1) * t221 - t219 * rSges(3,2)) + (-t228 * t167 - t227 * t168 + t237 * t216 + t15 + t332 * t214 + (t330 * t217 + (t331 - t333) * t216) * t217) * t305 + (-t228 * t165 - t227 * t166 - t237 * t217 + t14 + t333 * t215 + (t331 * t216 + (t330 - t332) * t217) * t216) * t304) * qJD(2), t16 * qJD(2) + t7 * qJD(5) + t276 + (-0.4e1 * t324 + 0.2e1 * t256 * (-t221 * t323 + t267)) * qJD(3), -t24 * t300 / 0.2e1, -t277 + t2 * qJD(2) + t7 * qJD(3) + (t13 * t305 + t12 * t304 + t18 * t308 + t17 * t310 + (t216 * t41 - t217 * t40) * t307 - t222) * qJD(5) + (-t278 / 0.2e1 + (t44 * t67 + t59 * t70 + t71 * t77 + t72 * t75 + (-t216 * t65 - t217 * t64) * t136 - t11) * qJD(5)) * m(6); -t99 * qJD(2), -t276 + t79 * qJD(3) + t6 * qJD(5) + 0.4e1 * (-t314 / 0.4e1 - t316 / 0.4e1 - t317 / 0.4e1) * qJD(2) + ((-t221 * t54 + t301) * t318 + (-t221 * t68 + t275) * t319 + (-t221 * t89 + t272) * t320 + ((t216 * t76 + t217 * t78 + t44) * t318 + (t111 * t216 + t113 * t217 + t63) * t319 + (t151 * t216 + t153 * t217 + t80) * t320) * t219) * t321, t79 * qJD(2), 0, t6 * qJD(2) + (-t221 * t67 + (t216 * t72 + t217 * t71) * t219) * t300; 0, ((-t216 * t78 + t217 * t76) * t318 + (t111 * t217 - t113 * t216) * t319) * t321 + qJD(5) * t265, 0, 0, t23 + (-t216 * t71 + t217 * t72) * t300; t34 * qJD(2), t277 + (t4 * t305 + t3 * t304 + (t216 * t48 - t217 * t47) * t309 + t15 * t308 + (t216 * t46 - t217 * t45) * t311 + t14 * t310 + (t216 * t53 - t217 * t52) * t306 + (t216 * t33 - t217 * t32) * t307 - t231) * qJD(2) + t8 * qJD(3) + t1 * qJD(5) + ((t43 * t44 + t49 * t77 + t50 * t75 + t54 * t59 + t64 * t78 + t65 * t76 - t22) * qJD(2) + t278 / 0.2e1) * m(6), t8 * qJD(2), t23, t1 * qJD(2) + (m(6) * (t59 * t67 + t64 * t71 + t65 * t72) + t13 * t308 + t12 * t310 + (-t40 * t166 - t41 * t168 + ((-t273 * t218 + t274 * t220) * t192 + t133 * t191) * t191) * t307) * qJD(5);];
Cq = t5;
