% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:57
% EndTime: 2020-01-03 11:36:12
% DurationCPUTime: 5.51s
% Computational Cost: add. (10877->372), mult. (9167->504), div. (0->0), fcn. (8673->10), ass. (0->194)
t186 = qJ(1) + pkin(8);
t181 = qJ(3) + t186;
t176 = cos(t181);
t188 = cos(pkin(9));
t275 = t176 * t188;
t187 = sin(pkin(9));
t276 = t176 * t187;
t116 = pkin(4) * t275 + pkin(7) * t276;
t175 = sin(t181);
t128 = pkin(3) * t176 + qJ(4) * t175;
t185 = qJD(1) + qJD(3);
t120 = t185 * t128;
t254 = qJD(4) * t176;
t277 = t176 * t185;
t280 = t175 * t185;
t264 = pkin(3) * t277 + qJ(4) * t280;
t230 = -t254 + t264;
t249 = t185 * t275;
t274 = t185 * t187;
t250 = t176 * t274;
t265 = pkin(4) * t249 + pkin(7) * t250;
t160 = -qJD(5) * t188 + t185;
t189 = sin(qJ(5));
t191 = cos(qJ(5));
t121 = -rSges(6,3) * t188 + (rSges(6,1) * t191 - rSges(6,2) * t189) * t187;
t253 = qJD(5) * t187;
t248 = t121 * t253;
t213 = (-qJD(4) - t248) * t176;
t273 = t188 * t189;
t111 = -t175 * t191 + t176 * t273;
t272 = t188 * t191;
t112 = t175 * t189 + t176 * t272;
t75 = rSges(6,1) * t112 - rSges(6,2) * t111 + rSges(6,3) * t276;
t304 = t160 * t75 + t213;
t110 = t175 * t272 - t176 * t189;
t81 = -qJD(5) * t110 - t111 * t185;
t109 = -t175 * t273 - t176 * t191;
t82 = qJD(5) * t109 + t112 * t185;
t48 = rSges(6,1) * t82 + rSges(6,2) * t81 + rSges(6,3) * t250;
t313 = -t116 * t185 - t120 + t230 + t265 - t304 + t48;
t262 = rSges(5,1) * t275 + rSges(5,3) * t175;
t100 = -rSges(5,2) * t276 + t262;
t266 = rSges(5,1) * t249 + rSges(5,3) * t280;
t312 = -t100 * t185 - t120 + t264 + t266;
t66 = Icges(6,5) * t112 - Icges(6,6) * t111 + Icges(6,3) * t276;
t104 = Icges(6,4) * t112;
t70 = Icges(6,2) * t111 - Icges(6,6) * t276 - t104;
t103 = Icges(6,4) * t111;
t72 = Icges(6,1) * t112 + Icges(6,5) * t276 - t103;
t31 = -(t189 * t70 + t191 * t72) * t187 + t188 * t66;
t226 = t111 * t70 + t112 * t72;
t27 = t276 * t66 + t226;
t311 = t176 * t27;
t301 = rSges(4,1) * t175 + rSges(4,2) * t176;
t107 = t301 * t185;
t179 = sin(t186);
t190 = sin(qJ(1));
t183 = t190 * pkin(1);
t300 = pkin(2) * t179 + t183;
t297 = t300 * qJD(1);
t94 = t297 + t107;
t310 = -t109 * t70 + t110 * t72;
t117 = -Icges(6,3) * t188 + (Icges(6,5) * t191 - Icges(6,6) * t189) * t187;
t281 = Icges(6,4) * t191;
t118 = -Icges(6,6) * t188 + (-Icges(6,2) * t189 + t281) * t187;
t282 = Icges(6,4) * t189;
t119 = -Icges(6,5) * t188 + (Icges(6,1) * t191 - t282) * t187;
t204 = -t111 * t118 + t112 * t119 + t117 * t276;
t306 = t204 * t160;
t169 = t175 * pkin(3);
t126 = -qJ(4) * t176 + t169;
t278 = t175 * t188;
t149 = rSges(5,1) * t278;
t279 = t175 * t187;
t235 = -rSges(5,2) * t279 + t149;
t302 = (-rSges(5,3) * t176 + t126 + t235) * t185;
t247 = t175 * t253;
t271 = pkin(4) * t278 + pkin(7) * t279 + t126;
t74 = rSges(6,1) * t110 + rSges(6,2) * t109 + rSges(6,3) * t279;
t296 = -t121 * t247 + t160 * t74 + t185 * t271;
t102 = Icges(6,4) * t109;
t71 = Icges(6,1) * t110 + Icges(6,5) * t279 + t102;
t207 = (-Icges(6,2) * t110 + t102 + t71) * t175 - (Icges(6,2) * t112 + t103 - t72) * t176;
t283 = Icges(6,4) * t110;
t68 = Icges(6,2) * t109 + Icges(6,6) * t279 + t283;
t295 = (-Icges(6,1) * t109 + t283 + t68) * t175 - (-Icges(6,1) * t111 - t104 + t70) * t176;
t193 = qJD(1) ^ 2;
t180 = cos(t186);
t174 = pkin(2) * t180;
t294 = pkin(4) * t188;
t192 = cos(qJ(1));
t184 = t192 * pkin(1);
t293 = -t111 * t68 + t112 * t71;
t288 = pkin(1) * qJD(1);
t286 = t175 * t66;
t65 = Icges(6,5) * t110 + Icges(6,6) * t109 + Icges(6,3) * t279;
t285 = t176 * t65;
t270 = t116 + t128;
t132 = (-Icges(6,1) * t189 - t281) * t187;
t269 = t118 - t132;
t131 = (-Icges(6,2) * t191 - t282) * t187;
t268 = t119 + t131;
t251 = t175 * t274;
t267 = rSges(5,2) * t251 + rSges(5,3) * t277;
t161 = qJD(4) * t175;
t263 = qJ(4) * t277 + t161;
t177 = t193 * t184;
t261 = t174 * t193 + t177;
t255 = qJD(1) * t180;
t260 = pkin(2) * t255 + t192 * t288;
t257 = t174 + t184;
t256 = qJD(1) * t179;
t24 = t109 * t68 + t110 * t71 + t279 * t65;
t25 = -t279 * t66 - t310;
t252 = t185 * t230 + t261;
t244 = -t253 / 0.2e1;
t243 = t253 / 0.2e1;
t241 = t128 + t262;
t240 = t175 * t244;
t239 = t175 * t243;
t238 = t176 * t244;
t237 = t176 * t243;
t236 = t185 * t243;
t108 = rSges(4,1) * t277 - rSges(4,2) * t280;
t79 = qJD(5) * t112 + t109 * t185;
t80 = qJD(5) * t111 + t110 * t185;
t234 = rSges(6,1) * t80 + rSges(6,2) * t79;
t232 = -rSges(5,2) * t274 - qJD(4);
t129 = rSges(4,1) * t176 - rSges(4,2) * t175;
t95 = t129 * t185 + t260;
t229 = -t254 + t260;
t139 = rSges(3,1) * t180 - rSges(3,2) * t179;
t206 = t297 - t161;
t33 = t206 + t296;
t34 = t185 * t270 + t260 + t304;
t225 = -t175 * t33 - t176 * t34;
t224 = -t175 * t75 + t176 * t74;
t223 = t175 * (Icges(6,5) * t109 - Icges(6,6) * t110) - t176 * (Icges(6,5) * t111 + Icges(6,6) * t112);
t222 = t175 * t236;
t221 = t176 * t236;
t220 = t169 + t235;
t219 = t300 * t193;
t214 = (-rSges(6,1) * t189 - rSges(6,2) * t191) * t187;
t212 = (t175 * t24 - t176 * t25) * t187;
t26 = -t276 * t65 - t293;
t211 = (t175 * t26 - t311) * t187;
t210 = -t234 + t263;
t130 = (-Icges(6,5) * t189 - Icges(6,6) * t191) * t187;
t209 = t161 * t185 - t219;
t205 = t74 + t271;
t201 = t75 + t270;
t200 = t34 * (-t294 - pkin(3) + (-rSges(6,3) - pkin(7)) * t187) * t280;
t199 = (-rSges(5,1) * t188 - pkin(3)) * t280 + t263 + t267;
t10 = qJD(5) * t211 - t306;
t42 = Icges(6,5) * t82 + Icges(6,6) * t81 + Icges(6,3) * t250;
t44 = Icges(6,4) * t82 + Icges(6,2) * t81 + Icges(6,6) * t250;
t46 = Icges(6,1) * t82 + Icges(6,4) * t81 + Icges(6,5) * t250;
t13 = -t188 * t42 + (-t189 * t44 + t191 * t46 + (-t189 * t71 - t191 * t68) * qJD(5)) * t187;
t41 = Icges(6,5) * t80 + Icges(6,6) * t79 + Icges(6,3) * t251;
t43 = Icges(6,4) * t80 + Icges(6,2) * t79 + Icges(6,6) * t251;
t45 = Icges(6,1) * t80 + Icges(6,4) * t79 + Icges(6,5) * t251;
t14 = -t188 * t41 + (-t189 * t43 + t191 * t45 + (t189 * t72 - t191 * t70) * qJD(5)) * t187;
t122 = qJD(5) * t130;
t123 = qJD(5) * t131;
t124 = qJD(5) * t132;
t19 = t111 * t123 - t112 * t124 + t118 * t79 + t119 * t80 + (t117 * t280 - t122 * t176) * t187;
t20 = t109 * t123 + t110 * t124 + t118 * t81 + t119 * t82 + (t117 * t277 + t122 * t175) * t187;
t30 = -t188 * t65 + (-t189 * t68 + t191 * t71) * t187;
t39 = -t122 * t188 + (-t123 * t189 + t124 * t191 + (-t118 * t191 - t119 * t189) * qJD(5)) * t187;
t36 = t39 * t160;
t49 = t109 * t118 + t110 * t119 + t117 * t279;
t40 = t49 * t160;
t9 = qJD(5) * t212 + t40;
t198 = t36 + (t40 + ((-t226 + t24 + t27) * t175 + (t26 + (t285 - t286) * t187 - t25 + t293) * t176) * t253) * t237 + (t13 + t20) * t239 + (t31 - t204) * t222 + (t30 + t49) * t221 + (t10 + t306 + (t311 + (t293 + t25 + (t285 + t286) * t187 + t310) * t175) * t253) * t240 + (t9 + t14 + t19) * t238;
t197 = ((t185 * t24 - t109 * t43 - t110 * t45 - t70 * t81 + t72 * t82 - (t175 * t41 - t277 * t66) * t187) * t176 + (t185 * t25 + t109 * t44 + t110 * t46 + t68 * t81 + t71 * t82 + (t175 * t42 + t277 * t65) * t187) * t175) * t187;
t196 = ((t185 * t26 - t111 * t43 + t112 * t45 - t70 * t79 + t72 * t80 - (-t176 * t41 - t280 * t66) * t187) * t176 + (t185 * t27 + t111 * t44 - t112 * t46 + t68 * t79 + t71 * t80 + (-t176 * t42 + t280 * t65) * t187) * t175) * t187;
t195 = ((t185 * t30 - t14) * t176 + (t185 * t31 + t13) * t175) * t187;
t97 = pkin(3) * t280 - t263;
t51 = (-t149 * t185 + t267 - t97) * t185 + t209;
t52 = (t176 * t232 + t266) * t185 + t252;
t60 = t206 + t302;
t194 = (-t51 * rSges(5,2) * t187 + t52 * (-rSges(5,3) - qJ(4)) + t60 * t232) * t176;
t125 = qJD(5) * t214;
t92 = t108 * t185 + t261;
t91 = -t107 * t185 - t219;
t90 = rSges(6,1) * t111 + rSges(6,2) * t112;
t89 = rSges(6,1) * t109 - rSges(6,2) * t110;
t61 = (t100 + t128) * t185 + t229;
t47 = rSges(6,3) * t251 + t234;
t35 = t224 * t253 + qJD(2);
t22 = -t125 * t247 + t160 * t48 + (t213 + t265) * t185 + t252;
t21 = -t125 * t176 * t253 - t160 * t47 + (-t97 + (t248 + (-pkin(7) * t187 - t294) * t185) * t175) * t185 + t209;
t16 = ((-t185 * t75 + t48) * t176 + (-t185 * t74 + t47) * t175) * t253;
t1 = [m(3) * (-t177 + t193 * (t139 + t184) + (-0.2e1 * rSges(3,1) * t255 + 0.2e1 * rSges(3,2) * t256 + qJD(1) * t139) * qJD(1)) * (-rSges(3,1) * t179 - rSges(3,2) * t180 - t183) + m(4) * (t91 * (t129 + t257) + t92 * (t300 + t301) + (-t95 + t108 + t260) * t94) + t198 + (t21 * (t201 + t257) + t34 * (-pkin(2) * t256 - t190 * t288 + t210) + t22 * (t205 + t300) + t200 + (t34 + t313) * t33) * m(6) + (t51 * (t241 + t257) + t61 * (-t297 + t199) + t52 * (t220 + t300) + t194 + (t260 - t229 + t61 + t312) * t60) * m(5); m(6) * t16; t198 + (t21 * t201 + t22 * t205 + t200 + (-t161 + t210 + t296) * t34 + t313 * t33) * m(6) + (t52 * t220 + t51 * t241 + t194 + (-t161 + t199 + t302) * t61 + (t254 + t312) * t60) * m(5) + (-t107 * t95 + t108 * t94 + t301 * t92 + t129 * t91 - (t129 * t94 - t301 * t95) * t185) * m(4); m(5) * (-t175 * t52 - t176 * t51) + m(6) * (-t22 * t175 - t21 * t176); -t188 * (qJD(5) * t195 + t36) / 0.2e1 + t160 * (-t188 * t39 + t195) / 0.2e1 + (qJD(5) * t197 + t160 * t20) * t279 / 0.2e1 + (-t188 * t49 + t212) * t221 + (-t188 * t20 + t197) * t239 - (qJD(5) * t196 + t160 * t19) * t276 / 0.2e1 + (t188 * t204 + t211) * t222 + (-t188 * t19 + t196) * t238 - t160 * (-t188 * t130 * t160 + ((-t189 * t268 - t191 * t269) * t160 + ((-t189 * t207 - t191 * t295) * t187 - t223 * t188) * qJD(5)) * t187) / 0.2e1 + ((t109 * t268 - t110 * t269 + t130 * t279) * t160 + (t109 * t207 - t110 * t295 + t223 * t279) * t253) * t240 + ((t111 * t268 + t112 * t269 - t130 * t276) * t160 + (t207 * t111 + t112 * t295 - t223 * t276) * t253) * t237 + (t10 * t175 + t176 * t9) * t274 / 0.2e1 + ((-t21 * t75 - t22 * t74 - t33 * t48 + t34 * t47) * t188 + (t16 * t224 + t35 * (t175 * t47 + t176 * t48 - t277 * t75 - t280 * t74) + t225 * t125 + ((-t185 * t33 - t21) * t176 + (t185 * t34 - t22) * t175) * t121) * t187 - (t33 * t89 - t34 * t90) * t160 - (t35 * (t175 * t90 + t176 * t89) + t225 * t214) * t253) * m(6);];
tauc = t1(:);
