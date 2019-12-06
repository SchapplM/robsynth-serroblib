% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:54
% EndTime: 2019-12-05 15:01:04
% DurationCPUTime: 5.56s
% Computational Cost: add. (14147->370), mult. (13578->594), div. (0->0), fcn. (14460->8), ass. (0->220)
t188 = sin(pkin(7));
t183 = t188 ^ 2;
t190 = cos(pkin(7));
t184 = t190 ^ 2;
t234 = t183 + t184;
t186 = pkin(8) + qJ(3);
t180 = sin(t186);
t182 = cos(t186);
t311 = Icges(4,5) * t180 + Icges(4,6) * t182;
t185 = pkin(9) + qJ(5);
t179 = sin(t185);
t181 = cos(t185);
t223 = rSges(6,1) * t181 - rSges(6,2) * t179;
t134 = -rSges(6,3) * t182 + t180 * t223;
t117 = t134 * t188;
t118 = t134 * t190;
t257 = t180 * t182;
t189 = cos(pkin(9));
t249 = t189 * t190;
t187 = sin(pkin(9));
t252 = t187 * t188;
t164 = -t182 * t252 - t249;
t250 = t188 * t189;
t251 = t187 * t190;
t165 = t182 * t250 - t251;
t166 = -t182 * t251 + t250;
t167 = t182 * t249 + t252;
t255 = t180 * t190;
t256 = t180 * t188;
t305 = (t188 * (Icges(5,5) * t167 + Icges(5,6) * t166 + Icges(5,3) * t255) - t190 * (Icges(5,5) * t165 + Icges(5,6) * t164 + Icges(5,3) * t256)) * t182 + (-t190 * (t187 * (Icges(5,4) * t165 + Icges(5,2) * t164 + Icges(5,6) * t256) - t189 * (Icges(5,1) * t165 + Icges(5,4) * t164 + Icges(5,5) * t256)) + (t187 * (Icges(5,4) * t167 + Icges(5,2) * t166 + Icges(5,6) * t255) - t189 * (Icges(5,1) * t167 + Icges(5,4) * t166 + Icges(5,5) * t255)) * t188) * t180 - t311 * t234;
t259 = (-Icges(6,5) * t179 - Icges(6,6) * t181) * t257;
t232 = m(5) / 0.4e1 + m(6) / 0.4e1;
t235 = t234 * t257;
t304 = t232 * (t235 - t257);
t168 = t234 * t180;
t264 = Icges(6,4) * t181;
t213 = -Icges(6,2) * t179 + t264;
t130 = -Icges(6,6) * t182 + t180 * t213;
t243 = -t130 + (-Icges(6,1) * t179 - t264) * t180;
t265 = Icges(6,4) * t179;
t216 = Icges(6,1) * t181 - t265;
t132 = -Icges(6,5) * t182 + t180 * t216;
t242 = t132 + (-Icges(6,2) * t181 - t265) * t180;
t211 = Icges(6,5) * t181 - Icges(6,6) * t179;
t128 = -Icges(6,3) * t182 + t180 * t211;
t224 = rSges(5,1) * t189 - rSges(5,2) * t187;
t302 = rSges(5,3) * t182 - t180 * t224;
t297 = 2 * qJD(3);
t296 = m(5) / 0.2e1;
t295 = m(6) / 0.2e1;
t169 = t180 * pkin(3) - t182 * qJ(4);
t241 = -t169 + t302;
t106 = t241 * t188;
t108 = t241 * t190;
t253 = t182 * t190;
t254 = t182 * t188;
t281 = t106 * t254 + t108 * t253;
t171 = t182 * pkin(3) + t180 * qJ(4);
t236 = t234 * t171;
t57 = t188 * (rSges(5,1) * t165 + rSges(5,2) * t164 + rSges(5,3) * t256) + t190 * (rSges(5,1) * t167 + rSges(5,2) * t166 + rSges(5,3) * t255) + t236;
t294 = m(5) * (t168 * t57 + t281);
t152 = -t179 * t254 - t181 * t190;
t248 = t190 * t179;
t153 = t181 * t254 - t248;
t85 = rSges(6,1) * t153 + rSges(6,2) * t152 + rSges(6,3) * t256;
t68 = t134 * t256 + t182 * t85;
t154 = t181 * t188 - t182 * t248;
t155 = t179 * t188 + t181 * t253;
t86 = rSges(6,1) * t155 + rSges(6,2) * t154 + rSges(6,3) * t255;
t69 = -t134 * t255 - t182 * t86;
t283 = t68 * t253 + t69 * t254;
t219 = -t188 * t86 + t190 * t85;
t50 = t219 * t182 + (-t117 * t190 + t118 * t188) * t180;
t135 = rSges(6,3) * t180 + t182 * t223;
t53 = (t135 * t188 - t85) * t180;
t54 = (-t135 * t190 + t86) * t180;
t60 = t219 * t180;
t293 = m(6) * (-t182 * t50 + (t188 * t54 + t190 * t53 + t60) * t180 + t283);
t191 = -pkin(6) - qJ(4);
t244 = qJ(4) + t191;
t178 = pkin(4) * t189 + pkin(3);
t284 = -pkin(3) + t178;
t229 = -t180 * t284 - t182 * t244 - t134 - t169;
t72 = t229 * t188;
t74 = t229 * t190;
t282 = t74 * t253 + t72 * t254;
t127 = -t180 * t244 + t182 * t284;
t43 = (t127 * t190 + t86) * t190 + (t127 * t188 + t85) * t188 + t236;
t292 = m(6) * (t168 * t43 + t282);
t291 = m(6) * (t168 * t60 + t283);
t147 = (-rSges(6,1) * t179 - rSges(6,2) * t181) * t180;
t104 = rSges(6,1) * t152 - rSges(6,2) * t153;
t105 = rSges(6,1) * t154 - rSges(6,2) * t155;
t67 = t104 * t188 + t105 * t190;
t290 = m(6) * (-t147 * t168 - t182 * t67);
t289 = t180 / 0.2e1;
t288 = -t182 / 0.2e1;
t287 = t188 / 0.2e1;
t286 = -t190 / 0.2e1;
t285 = t190 / 0.2e1;
t280 = m(6) * qJD(5);
t79 = Icges(6,5) * t153 + Icges(6,6) * t152 + Icges(6,3) * t256;
t277 = t182 * t79;
t80 = Icges(6,5) * t155 + Icges(6,6) * t154 + Icges(6,3) * t255;
t276 = t182 * t80;
t266 = Icges(6,4) * t155;
t82 = Icges(6,2) * t154 + Icges(6,6) * t255 + t266;
t143 = Icges(6,4) * t154;
t84 = Icges(6,1) * t155 + Icges(6,5) * t255 + t143;
t46 = t152 * t82 + t153 * t84 + t256 * t80;
t275 = t46 * t190;
t267 = Icges(6,4) * t153;
t81 = Icges(6,2) * t152 + Icges(6,6) * t256 + t267;
t142 = Icges(6,4) * t152;
t83 = Icges(6,1) * t153 + Icges(6,5) * t256 + t142;
t47 = t154 * t81 + t155 * t83 + t255 * t79;
t274 = t47 * t188;
t273 = -Icges(6,2) * t153 + t142 + t83;
t272 = -Icges(6,2) * t155 + t143 + t84;
t271 = Icges(6,1) * t152 - t267 - t81;
t270 = Icges(6,1) * t154 - t266 - t82;
t260 = t128 * t182;
t30 = t188 * t53 - t190 * t54;
t247 = t30 * qJD(2);
t38 = 0.2e1 * (t50 / 0.4e1 - t67 / 0.4e1) * m(6);
t246 = t38 * qJD(1);
t225 = 0.2e1 * t232 * t168;
t231 = t295 + t296;
t88 = -t180 * t231 + t225;
t245 = t88 * qJD(1);
t240 = -rSges(5,3) * t180 - t182 * t224 - t171;
t237 = t234 * t169;
t233 = qJD(5) * t180;
t230 = t30 * t295;
t228 = -t127 - t135 - t171;
t170 = rSges(4,1) * t180 + rSges(4,2) * t182;
t222 = -t179 * t81 + t181 * t83;
t221 = -t179 * t82 + t181 * t84;
t51 = t180 * t222 - t277;
t52 = t180 * t221 - t276;
t220 = t51 * t188 + t52 * t190;
t210 = -t130 * t179 + t132 * t181;
t209 = -t128 * t188 - t222;
t208 = -t128 * t190 - t221;
t98 = Icges(6,5) * t152 - Icges(6,6) * t153;
t32 = t152 * t273 + t153 * t271 + t256 * t98;
t99 = Icges(6,5) * t154 - Icges(6,6) * t155;
t33 = t152 * t272 + t153 * t270 + t256 * t99;
t17 = t188 * t33 - t190 * t32;
t34 = t154 * t273 + t155 * t271 + t255 * t98;
t35 = t154 * t272 + t155 * t270 + t255 * t99;
t18 = t188 * t35 - t190 * t34;
t206 = t17 * t286 + t18 * t287;
t203 = (Icges(6,3) * t180 + t182 * t211 - t210) * t182;
t200 = Icges(5,5) * t182 + (-Icges(5,1) * t189 + Icges(5,4) * t187) * t180;
t198 = Icges(5,6) * t182 + (-Icges(5,4) * t189 + Icges(5,2) * t187) * t180;
t194 = -t180 * t178 - t182 * t191 + t169;
t193 = t180 * t209 + t277;
t192 = t180 * t208 + t276;
t159 = t311 * t190;
t158 = t311 * t188;
t133 = Icges(6,5) * t180 + t182 * t216;
t131 = Icges(6,6) * t180 + t182 * t213;
t125 = t200 * t190;
t124 = t200 * t188;
t123 = t198 * t190;
t122 = t198 * t188;
t116 = t132 * t190;
t115 = t132 * t188;
t114 = t130 * t190;
t113 = t130 * t188;
t110 = t234 * t170;
t109 = t240 * t190;
t107 = t240 * t188;
t87 = t225 + (m(5) + m(6)) * t289;
t78 = 0.4e1 * t304;
t77 = -t105 * t182 - t147 * t255;
t76 = t104 * t182 + t147 * t256;
t75 = t228 * t190;
t73 = t228 * t188;
t63 = (t104 * t190 - t105 * t188) * t180;
t62 = t180 * t210 - t260;
t61 = t234 * t302 - t237;
t59 = t128 * t255 + t130 * t154 + t132 * t155;
t58 = t128 * t256 + t130 * t152 + t132 * t153;
t55 = t290 / 0.2e1;
t49 = (t190 * t194 - t118) * t190 + (t188 * t194 - t117) * t188 - t237;
t48 = t154 * t82 + t155 * t84 + t255 * t80;
t45 = t152 * t81 + t153 * t83 + t256 * t79;
t41 = -t182 * t99 + (-t179 * t272 + t181 * t270) * t180;
t40 = -t182 * t98 + (-t179 * t273 + t181 * t271) * t180;
t39 = (t50 + t67) * t295;
t37 = -t208 * t182 + (t114 * t179 - t116 * t181 + t80) * t180;
t36 = -t209 * t182 + (t113 * t179 - t115 * t181 + t79) * t180;
t29 = t291 / 0.2e1;
t28 = -t114 * t154 - t116 * t155 + t190 * t192;
t27 = -t113 * t154 - t115 * t155 + t190 * t193;
t26 = -t114 * t152 - t116 * t153 + t188 * t192;
t25 = -t113 * t152 - t115 * t153 + t188 * t193;
t24 = qJD(3) * t230;
t22 = t180 * t220 - t182 * t62;
t21 = t43 * t67 + (-t188 * t72 - t190 * t74) * t147;
t20 = -t182 * t59 + (t190 * t48 + t274) * t180;
t19 = -t182 * t58 + (t188 * t45 + t275) * t180;
t16 = t188 * t28 - t190 * t27;
t15 = t188 * t26 - t190 * t25;
t14 = t292 + t294;
t13 = t50 * t60 + t53 * t68 + t54 * t69;
t11 = t293 / 0.2e1;
t10 = -(t154 * t242 + t155 * t243) * t182 + (t34 * t188 + (t35 - t259) * t190) * t180;
t9 = -(t152 * t242 + t153 * t243) * t182 + (t33 * t190 + (t32 - t259) * t188) * t180;
t8 = (t203 + t220) * t182 + (t37 * t190 + t36 * t188 - (-t131 * t179 + t133 * t181 + t128) * t182 + t62) * t180;
t7 = t29 + t11 - t290 / 0.2e1;
t6 = t55 + t29 - t293 / 0.2e1;
t5 = t55 + t11 - t291 / 0.2e1;
t4 = (-t131 * t154 - t133 * t155 + t274 + (t48 - t260) * t190) * t182 + (t27 * t188 + t59 + (t28 - t203) * t190) * t180;
t3 = (-t131 * t152 - t133 * t153 + t275 + (t45 - t260) * t188) * t182 + (t26 * t190 + t58 + (t25 - t203) * t188) * t180;
t2 = m(6) * t21 + t206;
t1 = m(6) * t13 + (t20 * t285 + t19 * t287 - t8 / 0.2e1) * t182 + (t4 * t285 + t3 * t287 + t22 / 0.2e1) * t180;
t12 = [0, 0, t87 * qJD(4) + t39 * qJD(5) + (-m(4) * t110 / 0.2e1 + t61 * t296 + t49 * t295) * t297, t87 * qJD(3), t39 * qJD(3) + t280 * t63; 0, 0, ((-t107 * t190 + t109 * t188) * t296 + (t188 * t75 - t190 * t73) * t295) * t297 + qJD(5) * t230, 0, t24 + (t188 * t76 - t190 * t77) * t280; qJD(4) * t88 - qJD(5) * t38, -t30 * t280 / 0.2e1, t14 * qJD(4) + t2 * qJD(5) + (m(6) * (t43 * t49 + t72 * t73 + t74 * t75) + m(5) * (t106 * t107 + t108 * t109 + t57 * t61) + m(4) * (-t110 + t170) * t234 * (rSges(4,1) * t182 - rSges(4,2) * t180) + (t16 - t183 * t159 + (t166 * t123 + t167 * t125) * t188 + (-t166 * t122 - t167 * t124 + t188 * t158 + t305) * t190) * t287 + (t15 - t184 * t158 - (t164 * t122 + t165 * t124) * t190 + (t164 * t123 + t165 * t125 + t190 * t159 + t305) * t188) * t286) * qJD(3), t14 * qJD(3) + t6 * qJD(5) + t245 + (-0.4e1 * t304 + 0.2e1 * t231 * (-t168 * t182 + t235)) * qJD(4), -t246 + t2 * qJD(3) + t6 * qJD(4) - t247 * m(6) / 0.2e1 + (-t22 / 0.2e1 + (t18 / 0.2e1 - t4 / 0.2e1) * t190 + (t17 / 0.2e1 - t3 / 0.2e1) * t188) * t233 + (t10 * t287 + t9 * t286 + (t8 / 0.2e1 + (t40 / 0.2e1 - t20 / 0.2e1) * t190 + (-t41 / 0.2e1 - t19 / 0.2e1) * t188) * t182 + (t43 * t63 + t60 * t67 + t72 * t77 + t74 * t76 + (-t188 * t69 - t190 * t68) * t147 - t13) * m(6)) * qJD(5); -t88 * qJD(3), 0, -t245 + t78 * qJD(4) + t5 * qJD(5) + 0.4e1 * (-t292 / 0.4e1 - t294 / 0.4e1) * qJD(3) + ((-t182 * t49 + t282) * t295 + (-t182 * t61 + t281) * t296 + ((t188 * t73 + t190 * t75 + t43) * t295 + (t107 * t188 + t109 * t190 + t57) * t296) * t180) * t297, t78 * qJD(3), t5 * qJD(3) + (-t182 * t63 + (t188 * t77 + t190 * t76) * t180) * t280; t38 * qJD(3), t24, t246 + (t3 * t286 + (t188 * t52 - t190 * t51) * t289 + (t188 * t37 - t190 * t36) * t288 + (t188 * t48 - t190 * t47) * t253 / 0.2e1 + t16 * t255 / 0.2e1 + (t188 * t46 - t190 * t45) * t254 / 0.2e1 + t15 * t256 / 0.2e1 + t4 * t287 - t206) * qJD(3) + t7 * qJD(4) + t1 * qJD(5) + (t247 / 0.2e1 + (t43 * t50 + t49 * t60 + t53 * t74 + t54 * t72 + t68 * t75 + t69 * t73 - t21) * qJD(3)) * m(6), t7 * qJD(3), t1 * qJD(3) + (m(6) * (t60 * t63 + t68 * t76 + t69 * t77) - t182 ^ 2 * t259 / 0.2e1) * qJD(5) + (t10 * t285 + t9 * t287 + (t41 * t190 + t40 * t188 - (-t179 * t242 + t181 * t243) * t182) * t288) * t233;];
Cq = t12;
