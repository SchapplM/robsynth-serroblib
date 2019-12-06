% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:23:59
% EndTime: 2019-12-05 15:24:13
% DurationCPUTime: 7.27s
% Computational Cost: add. (14482->399), mult. (14130->635), div. (0->0), fcn. (14970->10), ass. (0->240)
t206 = sin(pkin(7));
t201 = t206 ^ 2;
t208 = cos(pkin(7));
t202 = t208 ^ 2;
t272 = t201 + t202;
t204 = qJ(2) + pkin(8);
t198 = sin(t204);
t200 = cos(t204);
t203 = pkin(9) + qJ(5);
t197 = sin(t203);
t199 = cos(t203);
t255 = rSges(6,1) * t199 - rSges(6,2) * t197;
t135 = -rSges(6,3) * t200 + t255 * t198;
t117 = t135 * t206;
t118 = t135 * t208;
t239 = -Icges(4,5) * t198 - Icges(4,6) * t200;
t211 = sin(qJ(2));
t212 = cos(qJ(2));
t351 = 0.2e1 * t212 * (Icges(3,1) - Icges(3,2)) * t211 + (-0.2e1 * t211 ^ 2 + 0.2e1 * t212 ^ 2) * Icges(3,4);
t289 = t198 * t200;
t240 = -Icges(3,5) * t211 - Icges(3,6) * t212;
t178 = t240 * t206;
t179 = t240 * t208;
t350 = t239 * t206 + t178;
t349 = t239 * t208 + t179;
t207 = cos(pkin(9));
t205 = sin(pkin(9));
t285 = t205 * t206;
t173 = -t200 * t285 - t207 * t208;
t284 = t206 * t207;
t174 = t200 * t284 - t205 * t208;
t283 = t208 * t200;
t175 = -t205 * t283 + t284;
t176 = t207 * t283 + t285;
t287 = t198 * t208;
t288 = t198 * t206;
t344 = ((Icges(5,5) * t176 + Icges(5,6) * t175 + Icges(5,3) * t287) * t206 - t208 * (Icges(5,5) * t174 + Icges(5,6) * t173 + Icges(5,3) * t288)) * t200 + (((Icges(5,1) * t174 + Icges(5,4) * t173 + Icges(5,5) * t288) * t207 - t205 * (Icges(5,4) * t174 + Icges(5,2) * t173 + Icges(5,6) * t288)) * t208 + (-(Icges(5,1) * t176 + Icges(5,4) * t175 + Icges(5,5) * t287) * t207 + t205 * (Icges(5,4) * t176 + Icges(5,2) * t175 + Icges(5,6) * t287)) * t206) * t198 + t272 * t239;
t291 = (-Icges(6,5) * t197 - Icges(6,6) * t199) * t289;
t270 = m(5) / 0.4e1 + m(6) / 0.4e1;
t273 = t272 * t289;
t343 = t270 * (t273 - t289);
t177 = t272 * t198;
t296 = Icges(6,4) * t199;
t241 = -Icges(6,2) * t197 + t296;
t131 = -Icges(6,6) * t200 + t241 * t198;
t278 = -t131 + (-Icges(6,1) * t197 - t296) * t198;
t297 = Icges(6,4) * t197;
t246 = Icges(6,1) * t199 - t297;
t133 = -Icges(6,5) * t200 + t246 * t198;
t277 = t133 + (-Icges(6,2) * t199 - t297) * t198;
t238 = Icges(6,5) * t199 - Icges(6,6) * t197;
t129 = -Icges(6,3) * t200 + t238 * t198;
t256 = rSges(5,1) * t207 - rSges(5,2) * t205;
t341 = rSges(5,3) * t200 - t256 * t198;
t336 = 2 * qJD(2);
t335 = m(4) / 0.2e1;
t334 = m(5) / 0.2e1;
t333 = m(6) / 0.2e1;
t286 = t200 * t206;
t184 = pkin(3) * t198 - qJ(4) * t200;
t322 = pkin(2) * t211;
t267 = -t184 - t322;
t259 = t267 + t341;
t92 = t259 * t206;
t94 = t259 * t208;
t315 = t94 * t283 + t92 * t286;
t186 = pkin(3) * t200 + qJ(4) * t198;
t321 = pkin(2) * t212;
t276 = t272 * t321;
t261 = t272 * t186 + t276;
t51 = t206 * (rSges(5,1) * t174 + rSges(5,2) * t173 + rSges(5,3) * t288) + t208 * (rSges(5,1) * t176 + rSges(5,2) * t175 + rSges(5,3) * t287) + t261;
t332 = m(5) * (t177 * t51 + t315);
t157 = -t197 * t286 - t199 * t208;
t158 = -t208 * t197 + t199 * t286;
t88 = rSges(6,1) * t158 + rSges(6,2) * t157 + rSges(6,3) * t288;
t70 = t135 * t288 + t200 * t88;
t159 = -t197 * t283 + t199 * t206;
t160 = t197 * t206 + t199 * t283;
t89 = rSges(6,1) * t160 + rSges(6,2) * t159 + rSges(6,3) * t287;
t71 = -t135 * t287 - t200 * t89;
t317 = t70 * t283 + t71 * t286;
t251 = -t206 * t89 + t208 * t88;
t50 = t251 * t200 + (-t117 * t208 + t118 * t206) * t198;
t136 = rSges(6,3) * t198 + t255 * t200;
t54 = (t206 * t136 - t88) * t198;
t55 = (-t208 * t136 + t89) * t198;
t61 = t251 * t198;
t331 = m(6) * (-t200 * t50 + (t206 * t55 + t208 * t54 + t61) * t198 + t317);
t210 = -pkin(6) - qJ(4);
t279 = qJ(4) + t210;
t195 = pkin(4) * t207 + pkin(3);
t319 = -pkin(3) + t195;
t236 = -t319 * t198 - t279 * t200 - t135 + t267;
t72 = t236 * t206;
t74 = t236 * t208;
t316 = t74 * t283 + t72 * t286;
t128 = -t279 * t198 + t319 * t200;
t41 = (t208 * t128 + t89) * t208 + (t128 * t206 + t88) * t206 + t261;
t330 = m(6) * (t177 * t41 + t316);
t329 = m(6) * (t177 * t61 + t317);
t152 = (-rSges(6,1) * t197 - rSges(6,2) * t199) * t198;
t109 = rSges(6,1) * t157 - rSges(6,2) * t158;
t110 = rSges(6,1) * t159 - rSges(6,2) * t160;
t67 = t109 * t206 + t110 * t208;
t328 = m(6) * (-t152 * t177 - t200 * t67);
t327 = t198 / 0.2e1;
t326 = -t200 / 0.2e1;
t325 = t206 / 0.2e1;
t324 = -t208 / 0.2e1;
t323 = t208 / 0.2e1;
t314 = m(6) * qJD(5);
t82 = Icges(6,5) * t158 + Icges(6,6) * t157 + Icges(6,3) * t288;
t311 = t200 * t82;
t83 = Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t287;
t310 = t200 * t83;
t298 = Icges(6,4) * t160;
t85 = Icges(6,2) * t159 + Icges(6,6) * t287 + t298;
t148 = Icges(6,4) * t159;
t87 = Icges(6,1) * t160 + Icges(6,5) * t287 + t148;
t47 = t157 * t85 + t158 * t87 + t83 * t288;
t309 = t47 * t208;
t299 = Icges(6,4) * t158;
t84 = Icges(6,2) * t157 + Icges(6,6) * t288 + t299;
t147 = Icges(6,4) * t157;
t86 = Icges(6,1) * t158 + Icges(6,5) * t288 + t147;
t48 = t159 * t84 + t160 * t86 + t82 * t287;
t308 = t48 * t206;
t307 = -Icges(6,2) * t158 + t147 + t86;
t306 = -Icges(6,2) * t160 + t148 + t87;
t305 = Icges(6,1) * t157 - t299 - t84;
t304 = Icges(6,1) * t159 - t298 - t85;
t292 = t129 * t200;
t30 = t206 * t54 - t208 * t55;
t282 = t30 * qJD(3);
t39 = 0.2e1 * (t50 / 0.4e1 - t67 / 0.4e1) * m(6);
t281 = t39 * qJD(1);
t257 = 0.2e1 * t270 * t177;
t269 = t333 + t334;
t91 = -t269 * t198 + t257;
t280 = t91 * qJD(1);
t271 = qJD(5) * t198;
t268 = t30 * t333;
t185 = rSges(4,1) * t198 + rSges(4,2) * t200;
t266 = -t185 - t322;
t265 = -t186 - t321;
t187 = rSges(4,1) * t200 - rSges(4,2) * t198;
t264 = -t187 - t321;
t260 = t272 * t322;
t258 = -rSges(5,3) * t198 - t256 * t200 + t265;
t190 = t211 * rSges(3,1) + rSges(3,2) * t212;
t254 = -t197 * t84 + t199 * t86;
t253 = -t197 * t85 + t199 * t87;
t52 = t254 * t198 - t311;
t53 = t253 * t198 - t310;
t252 = t52 * t206 + t53 * t208;
t237 = -t131 * t197 + t133 * t199;
t235 = -t128 - t136 + t265;
t234 = -t129 * t206 - t254;
t233 = -t129 * t208 - t253;
t103 = Icges(6,5) * t157 - Icges(6,6) * t158;
t32 = t103 * t288 + t307 * t157 + t305 * t158;
t104 = Icges(6,5) * t159 - Icges(6,6) * t160;
t33 = t104 * t288 + t306 * t157 + t304 * t158;
t17 = t206 * t33 - t208 * t32;
t34 = t103 * t287 + t307 * t159 + t305 * t160;
t35 = t104 * t287 + t306 * t159 + t304 * t160;
t18 = t206 * t35 - t208 * t34;
t232 = t17 * t324 + t18 * t325;
t81 = -t272 * t185 - t260;
t230 = t81 * t187;
t228 = -t272 * t184 - t260;
t226 = (Icges(6,3) * t198 + t238 * t200 - t237) * t200;
t225 = t351 * t206 + t179;
t224 = -t351 * t208 + t178;
t221 = Icges(5,5) * t200 + (-Icges(5,1) * t207 + Icges(5,4) * t205) * t198;
t219 = Icges(5,6) * t200 + (-Icges(5,4) * t207 + Icges(5,2) * t205) * t198;
t215 = -t195 * t198 - t200 * t210 + t184;
t214 = t234 * t198 + t311;
t213 = t233 * t198 + t310;
t142 = t264 * t208;
t141 = t264 * t206;
t134 = Icges(6,5) * t198 + t246 * t200;
t132 = Icges(6,6) * t198 + t241 * t200;
t126 = t221 * t208;
t125 = t221 * t206;
t124 = t219 * t208;
t123 = t219 * t206;
t119 = t272 * t190;
t116 = t133 * t208;
t115 = t133 * t206;
t114 = t131 * t208;
t113 = t131 * t206;
t95 = t258 * t208;
t93 = t258 * t206;
t90 = t257 + (m(5) + m(6)) * t327;
t78 = 0.4e1 * t343;
t77 = -t110 * t200 - t152 * t287;
t76 = t109 * t200 + t152 * t288;
t75 = t235 * t208;
t73 = t235 * t206;
t63 = (t109 * t208 - t110 * t206) * t198;
t62 = t237 * t198 - t292;
t60 = t272 * t341 + t228;
t59 = t129 * t287 + t131 * t159 + t133 * t160;
t58 = t129 * t288 + t131 * t157 + t133 * t158;
t56 = t328 / 0.2e1;
t49 = t159 * t85 + t160 * t87 + t83 * t287;
t46 = t157 * t84 + t158 * t86 + t82 * t288;
t44 = (t215 * t208 - t118) * t208 + (t206 * t215 - t117) * t206 + t228;
t43 = -t104 * t200 + (-t306 * t197 + t304 * t199) * t198;
t42 = -t103 * t200 + (-t307 * t197 + t305 * t199) * t198;
t40 = (t50 + t67) * t333;
t38 = -t233 * t200 + (t114 * t197 - t116 * t199 + t83) * t198;
t37 = -t234 * t200 + (t113 * t197 - t115 * t199 + t82) * t198;
t29 = t329 / 0.2e1;
t28 = -t114 * t159 - t116 * t160 + t213 * t208;
t27 = -t113 * t159 - t115 * t160 + t214 * t208;
t26 = -t114 * t157 - t116 * t158 + t206 * t213;
t25 = -t113 * t157 - t115 * t158 + t206 * t214;
t24 = qJD(2) * t268;
t23 = t252 * t198 - t200 * t62;
t21 = -t200 * t59 + (t208 * t49 + t308) * t198;
t20 = -t200 * t58 + (t206 * t46 + t309) * t198;
t19 = t41 * t67 + (-t206 * t72 - t208 * t74) * t152;
t16 = t206 * t28 - t208 * t27;
t15 = t206 * t26 - t208 * t25;
t14 = t330 + t332;
t13 = t50 * t61 + t54 * t70 + t55 * t71;
t11 = t331 / 0.2e1;
t10 = -(t277 * t159 + t278 * t160) * t200 + (t34 * t206 + (t35 - t291) * t208) * t198;
t9 = -(t277 * t157 + t278 * t158) * t200 + (t33 * t208 + (t32 - t291) * t206) * t198;
t8 = (t226 + t252) * t200 + (t38 * t208 + t37 * t206 - (-t132 * t197 + t134 * t199 + t129) * t200 + t62) * t198;
t7 = t29 + t11 - t328 / 0.2e1;
t6 = t56 + t29 - t331 / 0.2e1;
t5 = t56 + t11 - t329 / 0.2e1;
t4 = (-t132 * t159 - t134 * t160 + t308 + (t49 - t292) * t208) * t200 + (t27 * t206 + t59 + (t28 - t226) * t208) * t198;
t3 = (-t132 * t157 - t134 * t158 + t309 + (t46 - t292) * t206) * t200 + (t26 * t208 + t58 + (t25 - t226) * t206) * t198;
t2 = m(6) * t19 + t232;
t1 = m(6) * t13 + (t21 * t323 + t20 * t325 - t8 / 0.2e1) * t200 + (t4 * t323 + t3 * t325 + t23 / 0.2e1) * t198;
t12 = [0, t90 * qJD(4) + t40 * qJD(5) + (-m(3) * t119 / 0.2e1 + t81 * t335 + t60 * t334 + t44 * t333) * t336, 0, t90 * qJD(2), t40 * qJD(2) + t63 * t314; qJD(4) * t91 - qJD(5) * t39, t14 * qJD(4) + t2 * qJD(5) + (m(6) * (t41 * t44 + t72 * t73 + t74 * t75) + m(5) * (t51 * t60 + t92 * t93 + t94 * t95) + m(4) * (t276 * t81 + (t266 * t142 + t208 * t230) * t208 + (t266 * t141 + t206 * t230) * t206) + m(3) * (-t119 + t190) * t272 * (rSges(3,1) * t212 - t211 * rSges(3,2)) + (t16 + (t124 * t175 + t126 * t176) * t206 + t349 * t201 + (-t175 * t123 - t176 * t125 + t225 * t208 + (t224 - t350) * t206 + t344) * t208) * t325 + (t15 - (t123 * t173 + t125 * t174) * t208 + t350 * t202 + (t173 * t124 + t174 * t126 + t224 * t206 + (t225 - t349) * t208 + t344) * t206) * t324) * qJD(2), -t30 * t314 / 0.2e1, t14 * qJD(2) + t6 * qJD(5) + t280 + (-0.4e1 * t343 + 0.2e1 * t269 * (-t177 * t200 + t273)) * qJD(4), -t281 + t2 * qJD(2) + t6 * qJD(4) - t282 * m(6) / 0.2e1 + (-t23 / 0.2e1 + (t18 / 0.2e1 - t4 / 0.2e1) * t208 + (t17 / 0.2e1 - t3 / 0.2e1) * t206) * t271 + (t10 * t325 + t9 * t324 + (t8 / 0.2e1 + (t42 / 0.2e1 - t21 / 0.2e1) * t208 + (-t43 / 0.2e1 - t20 / 0.2e1) * t206) * t200 + (t41 * t63 + t61 * t67 + t72 * t77 + t74 * t76 + (-t206 * t71 - t208 * t70) * t152 - t13) * m(6)) * qJD(5); 0, ((t206 * t75 - t208 * t73) * t333 + (t206 * t95 - t208 * t93) * t334 + (-t141 * t208 + t142 * t206) * t335) * t336 + qJD(5) * t268, 0, 0, t24 + (t206 * t76 - t208 * t77) * t314; -t91 * qJD(2), -t280 + t78 * qJD(4) + t5 * qJD(5) + 0.4e1 * (-t330 / 0.4e1 - t332 / 0.4e1) * qJD(2) + ((-t200 * t44 + t316) * t333 + (-t200 * t60 + t315) * t334 + ((t206 * t73 + t208 * t75 + t41) * t333 + (t206 * t93 + t208 * t95 + t51) * t334) * t198) * t336, 0, t78 * qJD(2), t5 * qJD(2) + (-t200 * t63 + (t206 * t77 + t208 * t76) * t198) * t314; t39 * qJD(2), t281 + (t4 * t325 + t3 * t324 + (t206 * t53 - t208 * t52) * t327 + (t206 * t38 - t208 * t37) * t326 + (t206 * t49 - t208 * t48) * t283 / 0.2e1 + t16 * t287 / 0.2e1 + (t206 * t47 - t208 * t46) * t286 / 0.2e1 + t15 * t288 / 0.2e1 - t232) * qJD(2) + t7 * qJD(4) + t1 * qJD(5) + ((t41 * t50 + t44 * t61 + t54 * t74 + t55 * t72 + t70 * t75 + t71 * t73 - t19) * qJD(2) + t282 / 0.2e1) * m(6), t24, t7 * qJD(2), t1 * qJD(2) + (m(6) * (t61 * t63 + t70 * t76 + t71 * t77) - t200 ^ 2 * t291 / 0.2e1) * qJD(5) + (t10 * t323 + t9 * t325 + (t43 * t208 + t42 * t206 - (-t277 * t197 + t278 * t199) * t200) * t326) * t271;];
Cq = t12;
