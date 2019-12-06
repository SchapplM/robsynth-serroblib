% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:00
% EndTime: 2019-12-05 15:03:09
% DurationCPUTime: 5.03s
% Computational Cost: add. (9563->328), mult. (12964->532), div. (0->0), fcn. (13786->6), ass. (0->203)
t176 = pkin(8) + qJ(3);
t172 = sin(t176);
t173 = cos(t176);
t301 = -Icges(4,5) * t172 + (Icges(5,5) - Icges(4,6)) * t173;
t177 = sin(pkin(7));
t174 = t177 ^ 2;
t178 = cos(pkin(7));
t175 = t178 ^ 2;
t221 = t174 + t175;
t300 = Icges(5,4) * t172 + t301;
t179 = sin(qJ(5));
t180 = cos(qJ(5));
t208 = rSges(6,1) * t179 + rSges(6,2) * t180;
t186 = -t172 * rSges(6,3) + t173 * t208;
t112 = t186 * t177;
t113 = t186 * t178;
t244 = t172 * t173;
t298 = -Icges(5,2) + Icges(5,3);
t297 = t300 * t177;
t296 = t300 * t178;
t240 = t173 * t178;
t241 = t173 * t177;
t242 = t172 * t178;
t243 = t172 * t177;
t291 = 2 * Icges(5,6);
t292 = (-(-Icges(5,4) * t178 + t298 * t241 + t243 * t291) * t178 + (Icges(5,4) * t177 + t298 * t240 + t242 * t291) * t177) * t172 + t301 * t221;
t245 = (-Icges(6,5) * t180 + Icges(6,6) * t179) * t244;
t218 = m(5) / 0.4e1 + m(6) / 0.4e1;
t222 = t221 * t244;
t290 = t218 * (t222 - t244);
t289 = t221 * t172;
t253 = Icges(6,4) * t179;
t196 = Icges(6,2) * t180 + t253;
t184 = -Icges(6,6) * t172 + t173 * t196;
t232 = -t184 - (-Icges(6,1) * t180 + t253) * t173;
t252 = Icges(6,4) * t180;
t200 = Icges(6,1) * t179 + t252;
t185 = -Icges(6,5) * t172 + t173 * t200;
t231 = -t185 + (Icges(6,2) * t179 - t252) * t173;
t194 = Icges(6,5) * t179 + Icges(6,6) * t180;
t183 = -Icges(6,3) * t172 + t173 * t194;
t283 = 2 * qJD(3);
t282 = m(5) / 0.2e1;
t281 = m(6) / 0.2e1;
t159 = pkin(3) * t172 - qJ(4) * t173;
t207 = rSges(5,2) * t172 + rSges(5,3) * t173;
t224 = -t159 + t207;
t102 = t224 * t177;
t104 = t224 * t178;
t264 = t102 * t241 + t104 * t240;
t163 = -rSges(5,2) * t173 + rSges(5,3) * t172;
t162 = pkin(3) * t173 + qJ(4) * t172;
t225 = t221 * t162;
t60 = t163 * t221 + t225;
t280 = m(5) * (t289 * t60 + t264);
t237 = t178 * t179;
t238 = t177 * t180;
t156 = t172 * t238 + t237;
t236 = t178 * t180;
t239 = t177 * t179;
t157 = -t172 * t239 + t236;
t89 = -rSges(6,1) * t157 + rSges(6,2) * t156 + rSges(6,3) * t241;
t67 = -t172 * t89 - t186 * t241;
t154 = t172 * t236 - t239;
t155 = t172 * t237 + t238;
t88 = rSges(6,1) * t155 + rSges(6,2) * t154 + rSges(6,3) * t240;
t68 = t172 * t88 + t186 * t240;
t270 = t67 * t240 + t68 * t241;
t205 = t177 * t88 - t178 * t89;
t44 = (t112 * t178 - t113 * t177) * t173 + t205 * t172;
t120 = t173 * rSges(6,3) + t172 * t208;
t52 = (t120 * t177 - t89) * t173;
t53 = (-t120 * t178 + t88) * t173;
t59 = t205 * t173;
t279 = m(6) * (-t173 * t44 + (t177 * t53 + t178 * t52 - t59) * t172 + t270);
t210 = -pkin(6) * t172 - t159 + t186;
t77 = t210 * t177;
t79 = t210 * t178;
t269 = t79 * t240 + t77 * t241;
t271 = pkin(6) * t173;
t49 = t177 * t89 + t178 * t88 + t221 * t271 + t225;
t278 = m(6) * (t289 * t49 + t269);
t277 = m(6) * (-t289 * t59 + t270);
t153 = (-rSges(6,1) * t180 + rSges(6,2) * t179) * t173;
t100 = rSges(6,1) * t156 + rSges(6,2) * t157;
t99 = rSges(6,1) * t154 - rSges(6,2) * t155;
t69 = t100 * t177 + t178 * t99;
t276 = m(6) * (-t153 * t289 - t173 * t69);
t275 = t172 / 0.2e1;
t274 = t177 / 0.2e1;
t273 = -t178 / 0.2e1;
t272 = t178 / 0.2e1;
t255 = Icges(6,4) * t155;
t84 = Icges(6,2) * t154 + Icges(6,6) * t240 + t255;
t268 = -Icges(6,1) * t154 + t255 + t84;
t254 = Icges(6,4) * t157;
t85 = Icges(6,2) * t156 + Icges(6,6) * t241 - t254;
t267 = -Icges(6,1) * t156 - t254 + t85;
t148 = Icges(6,4) * t154;
t86 = Icges(6,1) * t155 + Icges(6,5) * t240 + t148;
t266 = -Icges(6,2) * t155 + t148 + t86;
t149 = Icges(6,4) * t156;
t87 = -Icges(6,1) * t157 + Icges(6,5) * t241 + t149;
t265 = Icges(6,2) * t157 + t149 + t87;
t263 = m(6) * qJD(5);
t82 = Icges(6,5) * t155 + Icges(6,6) * t154 + Icges(6,3) * t240;
t261 = t172 * t82;
t83 = -Icges(6,5) * t157 + Icges(6,6) * t156 + Icges(6,3) * t241;
t260 = t172 * t83;
t46 = t154 * t85 + t155 * t87 + t83 * t240;
t259 = t46 * t177;
t47 = t156 * t84 - t157 * t86 + t82 * t241;
t258 = t47 * t178;
t246 = t172 * t183;
t26 = t177 * t52 - t178 * t53;
t235 = t26 * qJD(2);
t36 = 0.2e1 * (t44 / 0.4e1 - t69 / 0.4e1) * m(6);
t234 = t36 * qJD(1);
t209 = 0.2e1 * t218 * t289;
t219 = t282 + t281;
t76 = -t172 * t219 + t209;
t233 = t76 * qJD(1);
t226 = t221 * t159;
t223 = -t162 - t163;
t220 = qJD(5) * t173;
t217 = t26 * t281;
t211 = -t120 - t162 - t271;
t161 = rSges(4,1) * t172 + rSges(4,2) * t173;
t204 = -t179 * t86 - t180 * t84;
t50 = t173 * t204 + t261;
t203 = -t179 * t87 - t180 * t85;
t51 = t173 * t203 + t260;
t206 = t51 * t177 + t50 * t178;
t191 = t179 * t185 + t180 * t184;
t190 = t183 * t177 - t203;
t189 = t183 * t178 - t204;
t93 = Icges(6,5) * t154 - Icges(6,6) * t155;
t32 = t266 * t154 - t268 * t155 + t93 * t240;
t94 = Icges(6,5) * t156 + Icges(6,6) * t157;
t33 = t265 * t154 - t267 * t155 + t94 * t240;
t16 = t177 * t32 - t178 * t33;
t34 = t266 * t156 + t268 * t157 + t93 * t241;
t35 = t265 * t156 + t267 * t157 + t94 * t241;
t17 = t177 * t34 - t178 * t35;
t188 = t16 * t274 + t17 * t273;
t187 = (Icges(6,3) * t173 + t172 * t194 - t191) * t172;
t182 = t173 * t189 - t261;
t181 = t173 * t190 - t260;
t118 = Icges(6,5) * t173 + t172 * t200;
t116 = Icges(6,6) * t173 + t172 * t196;
t111 = t185 * t178;
t110 = t185 * t177;
t109 = t184 * t178;
t108 = t184 * t177;
t105 = t223 * t178;
t103 = t223 * t177;
t90 = t221 * t161;
t80 = t211 * t178;
t78 = t211 * t177;
t75 = t209 + (m(5) + m(6)) * t275;
t72 = -t153 * t240 + t172 * t99;
t71 = -t100 * t172 + t153 * t241;
t70 = 0.4e1 * t290;
t65 = t207 * t221 - t226;
t64 = (t100 * t178 - t177 * t99) * t173;
t61 = t173 * t191 - t246;
t58 = -t156 * t184 + t157 * t185 - t183 * t241;
t57 = -t154 * t184 - t155 * t185 - t183 * t240;
t56 = -pkin(6) * t289 + t112 * t177 + t113 * t178 - t226;
t54 = t276 / 0.2e1;
t48 = t156 * t85 - t157 * t87 + t83 * t241;
t45 = t154 * t84 + t155 * t86 + t82 * t240;
t41 = t172 * t94 + (t267 * t179 - t265 * t180) * t173;
t40 = t172 * t93 + (t268 * t179 - t266 * t180) * t173;
t39 = (-t108 * t180 - t110 * t179 + t83) * t173 + t190 * t172;
t38 = (-t109 * t180 - t111 * t179 + t82) * t173 + t189 * t172;
t37 = (t44 + t69) * t281;
t31 = t154 * t108 + t155 * t110 + t178 * t181;
t30 = t154 * t109 + t155 * t111 + t178 * t182;
t29 = t156 * t108 - t157 * t110 + t177 * t181;
t28 = t156 * t109 - t157 * t111 + t177 * t182;
t25 = t277 / 0.2e1;
t23 = qJD(3) * t217;
t22 = t49 * t69 + (-t177 * t77 - t178 * t79) * t153;
t21 = t172 * t61 + t173 * t206;
t20 = t172 * t58 + (t177 * t48 + t258) * t173;
t19 = t172 * t57 + (t178 * t45 + t259) * t173;
t18 = t278 + t280;
t15 = t177 * t30 - t178 * t31;
t14 = t177 * t28 - t178 * t29;
t13 = (t231 * t156 + t232 * t157) * t172 + (t34 * t178 + (t35 + t245) * t177) * t173;
t12 = (t231 * t154 - t232 * t155) * t172 + (t33 * t177 + (t32 + t245) * t178) * t173;
t11 = -t44 * t59 + t52 * t67 + t53 * t68;
t9 = t279 / 0.2e1;
t8 = (t39 * t177 + t38 * t178 + t61) * t173 + (t187 + (-t116 * t180 - t118 * t179 - t183) * t173 - t206) * t172;
t7 = t25 + t9 - t276 / 0.2e1;
t6 = t54 + t25 - t279 / 0.2e1;
t5 = t54 + t9 - t277 / 0.2e1;
t4 = (t154 * t116 + t155 * t118 - t259 + (-t45 + t246) * t178) * t172 + (t31 * t177 + t57 + (t30 + t187) * t178) * t173;
t3 = (t156 * t116 - t157 * t118 - t258 + (-t48 + t246) * t177) * t172 + (t28 * t178 + t58 + (t29 + t187) * t177) * t173;
t2 = m(6) * t22 + t188;
t1 = m(6) * t11 + (t4 * t272 + t3 * t274 + t21 / 0.2e1) * t173 + (t19 * t273 - t177 * t20 / 0.2e1 + t8 / 0.2e1) * t172;
t10 = [0, 0, t75 * qJD(4) + t37 * qJD(5) + (-m(4) * t90 / 0.2e1 + t65 * t282 + t56 * t281) * t283, t75 * qJD(3), t37 * qJD(3) + t64 * t263; 0, 0, ((-t103 * t178 + t105 * t177) * t282 + (t177 * t80 - t178 * t78) * t281) * t283 + qJD(5) * t217, 0, t23 + (t177 * t71 - t178 * t72) * t263; qJD(4) * t76 - qJD(5) * t36, -t26 * t263 / 0.2e1, t18 * qJD(4) + t2 * qJD(5) + (m(5) * (t102 * t103 + t104 * t105 + t60 * t65) + m(6) * (t49 * t56 + t77 * t78 + t79 * t80) + m(4) * (t161 - t90) * t221 * (rSges(4,1) * t173 - rSges(4,2) * t172) + (t15 + (-t297 * t177 + t292) * t178 + t296 * t174) * t274 + (t14 + (-t296 * t178 + t292) * t177 + t297 * t175) * t273) * qJD(3), t18 * qJD(3) + t6 * qJD(5) + t233 + (-0.4e1 * t290 + 0.2e1 * t219 * (-t173 * t289 + t222)) * qJD(4), -t234 + t2 * qJD(3) + t6 * qJD(4) - t235 * m(6) / 0.2e1 + (-t21 / 0.2e1 + (t16 / 0.2e1 - t4 / 0.2e1) * t178 + (t17 / 0.2e1 - t3 / 0.2e1) * t177) * t220 + (t12 * t274 + t13 * t273 + (-t8 / 0.2e1 + (-t41 / 0.2e1 + t19 / 0.2e1) * t178 + (t40 / 0.2e1 + t20 / 0.2e1) * t177) * t172 + (t49 * t64 - t59 * t69 + t71 * t79 + t72 * t77 + (-t177 * t68 - t178 * t67) * t153 - t11) * m(6)) * qJD(5); -t76 * qJD(3), 0, -t233 + t70 * qJD(4) + t5 * qJD(5) + 0.4e1 * (-t280 / 0.4e1 - t278 / 0.4e1) * qJD(3) + ((-t173 * t65 + t264) * t282 + (-t173 * t56 + t269) * t281 + ((t103 * t177 + t105 * t178 + t60) * t282 + (t177 * t78 + t178 * t80 + t49) * t281) * t172) * t283, t70 * qJD(3), t5 * qJD(3) + (-t173 * t64 + (t177 * t72 + t178 * t71) * t172) * t263; t36 * qJD(3), t23, t234 + (t173 * (t177 * t50 - t178 * t51) / 0.2e1 + (t177 * t38 - t178 * t39) * t275 - (t177 * t45 - t178 * t46) * t242 / 0.2e1 + t15 * t240 / 0.2e1 - (t177 * t47 - t178 * t48) * t243 / 0.2e1 + t14 * t241 / 0.2e1 + t4 * t274 + t3 * t273 - t188) * qJD(3) + t7 * qJD(4) + t1 * qJD(5) + (t235 / 0.2e1 + (t44 * t49 + t52 * t79 + t53 * t77 - t56 * t59 + t67 * t80 + t68 * t78 - t22) * qJD(3)) * m(6), t7 * qJD(3), t1 * qJD(3) + (m(6) * (-t59 * t64 + t67 * t71 + t68 * t72) + t172 ^ 2 * t245 / 0.2e1) * qJD(5) + (t12 * t272 + t13 * t274 + (t40 * t178 + t41 * t177 + (t179 * t232 - t180 * t231) * t172) * t275) * t220;];
Cq = t10;
