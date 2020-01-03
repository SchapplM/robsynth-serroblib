% Calculate vector of inverse dynamics joint torques for
% S5RPRPR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:42
% EndTime: 2019-12-31 18:17:47
% DurationCPUTime: 4.00s
% Computational Cost: add. (7331->366), mult. (5409->470), div. (0->0), fcn. (3956->8), ass. (0->203)
t180 = qJ(1) + pkin(8);
t175 = qJ(3) + t180;
t168 = sin(t175);
t169 = cos(t175);
t179 = qJD(1) + qJD(3);
t236 = qJD(5) * t179;
t103 = qJDD(5) * t168 + t169 * t236;
t181 = sin(qJ(5));
t183 = cos(qJ(5));
t221 = rSges(6,1) * t181 + rSges(6,2) * t183;
t126 = t221 * qJD(5);
t263 = rSges(6,2) * t181;
t265 = rSges(6,1) * t183;
t149 = -t263 + t265;
t164 = t169 * pkin(7);
t177 = t179 ^ 2;
t178 = qJDD(1) + qJDD(3);
t173 = sin(t180);
t174 = cos(t180);
t185 = qJD(1) ^ 2;
t182 = sin(qJ(1));
t184 = cos(qJ(1));
t202 = (-qJDD(1) * t182 - t184 * t185) * pkin(1);
t189 = (-qJDD(1) * t173 - t174 * t185) * pkin(2) + t202;
t239 = qJD(4) * t179;
t188 = qJDD(4) * t168 + t169 * t239 + t189;
t156 = t169 * qJ(4);
t112 = pkin(3) * t168 - t156;
t271 = pkin(7) * t168;
t87 = -t168 * rSges(6,3) + t169 * t221;
t226 = -t112 + t87 - t271;
t238 = qJD(5) * t168;
t162 = t169 * rSges(6,3);
t232 = qJD(5) * t263;
t233 = qJD(5) * t265;
t208 = (t232 - t233) * t169;
t52 = (t168 * t221 + t162) * t179 + t208;
t165 = t169 * pkin(3);
t115 = t168 * qJ(4) + t165;
t153 = qJD(4) * t169;
t78 = t115 * t179 - t153;
t11 = -t177 * t164 - t126 * t238 + t103 * t149 + (-t52 - t78) * t179 + t226 * t178 + t188;
t299 = -g(1) + t11;
t248 = t169 * t179;
t137 = rSges(5,2) * t248;
t261 = rSges(5,3) * t169;
t113 = rSges(5,2) * t168 + t261;
t246 = -t112 + t113;
t251 = t168 * t179;
t24 = (-rSges(5,3) * t251 + t137 - t78) * t179 + t246 * t178 + t188;
t298 = -g(1) + t24;
t104 = qJDD(5) * t169 - t168 * t236;
t167 = pkin(2) * t174;
t176 = t184 * pkin(1);
t170 = qJDD(1) * t176;
t272 = pkin(1) * t182;
t224 = -pkin(2) * t173 - t272;
t193 = qJDD(1) * t167 + t185 * t224 + t170;
t152 = qJD(4) * t168;
t245 = qJ(4) * t248 + t152;
t191 = t178 * t115 + t168 * t239 + t179 * (-pkin(3) * t251 + t245) + t193;
t231 = t168 * t233 + t221 * t248;
t53 = (-rSges(6,3) * t179 - t232) * t168 + t231;
t249 = t168 * t183;
t250 = t168 * t181;
t86 = rSges(6,1) * t250 + rSges(6,2) * t249 + t162;
t12 = -t177 * t271 - t104 * t149 + t178 * t86 + t179 * t53 + (pkin(7) * t178 + qJD(5) * t126 - qJDD(4)) * t169 + t191;
t297 = -g(2) + t12;
t116 = -rSges(5,2) * t169 + t168 * rSges(5,3);
t244 = rSges(5,2) * t251 + rSges(5,3) * t248;
t25 = -qJDD(4) * t169 + t178 * t116 + t179 * t244 + t191;
t296 = -g(2) + t25;
t114 = rSges(4,1) * t168 + rSges(4,2) * t169;
t135 = rSges(4,2) * t251;
t91 = rSges(4,1) * t248 - t135;
t295 = -t114 * t178 - t179 * t91 - g(1) + t189;
t163 = t169 * rSges(4,1);
t117 = -rSges(4,2) * t168 + t163;
t247 = t179 * t114;
t294 = t117 * t178 - t179 * t247 - g(2) + t193;
t257 = Icges(6,4) * t183;
t144 = -Icges(6,2) * t181 + t257;
t213 = Icges(6,1) * t181 + t257;
t242 = t144 + t213;
t258 = Icges(6,4) * t181;
t146 = Icges(6,1) * t183 - t258;
t212 = Icges(6,2) * t183 + t258;
t243 = -t212 + t146;
t293 = (t181 * t242 - t183 * t243) * t179;
t109 = t149 * t238;
t75 = t179 * t87;
t292 = t109 + t75;
t211 = Icges(6,5) * t181 + Icges(6,6) * t183;
t291 = t211 * t179;
t83 = -Icges(6,6) * t168 + t169 * t212;
t85 = -Icges(6,5) * t168 + t169 * t213;
t214 = t181 * t85 + t183 * t83;
t290 = t214 * t169;
t120 = t174 * rSges(3,1) - rSges(3,2) * t173;
t108 = t120 + t176;
t289 = t167 + t176;
t287 = t164 + t86 + t115;
t44 = t181 * t83 - t183 * t85;
t204 = t212 * t179;
t284 = -Icges(6,6) * t179 + qJD(5) * t144;
t48 = t168 * t204 - t169 * t284;
t205 = t213 * t179;
t283 = -Icges(6,5) * t179 + qJD(5) * t146;
t50 = t168 * t205 - t169 * t283;
t81 = -Icges(6,3) * t168 + t169 * t211;
t286 = qJD(5) * t44 + t179 * t81 + t181 * t50 + t183 * t48;
t142 = Icges(6,5) * t183 - Icges(6,6) * t181;
t285 = -Icges(6,3) * t179 + qJD(5) * t142;
t124 = t212 * qJD(5);
t125 = t213 * qJD(5);
t209 = t181 * t144 - t146 * t183;
t282 = qJD(5) * t209 + t124 * t183 + t125 * t181 + t142 * t179;
t82 = Icges(6,6) * t169 + t168 * t212;
t138 = Icges(6,4) * t249;
t84 = Icges(6,1) * t250 + Icges(6,5) * t169 + t138;
t215 = t181 * t82 - t183 * t84;
t49 = t168 * t284 + t169 * t204;
t51 = t168 * t283 + t169 * t205;
t80 = Icges(6,3) * t169 + t168 * t211;
t281 = qJD(5) * t215 + t179 * t80 - t181 * t51 - t183 * t49;
t267 = t144 * t169 + t85;
t269 = -t146 * t169 + t83;
t279 = t181 * t269 - t183 * t267;
t268 = -Icges(6,2) * t250 + t138 + t84;
t270 = -t146 * t168 + t82;
t278 = t181 * t270 - t183 * t268;
t277 = -pkin(3) - pkin(7);
t276 = t103 / 0.2e1;
t275 = t104 / 0.2e1;
t274 = t168 / 0.2e1;
t273 = -t169 / 0.2e1;
t206 = t289 * qJD(1);
t74 = t117 * t179 + t206;
t260 = t114 * t74;
t210 = t183 * t144 + t181 * t146;
t94 = t168 * t142;
t60 = t169 * t210 - t94;
t259 = t60 * t179;
t252 = t142 * t169;
t77 = t115 + t116;
t105 = t179 * t112;
t241 = t152 - t105;
t237 = qJD(5) * t169;
t235 = m(3) + m(4) + m(5);
t234 = -rSges(6,3) + t277;
t26 = t169 * t80 + t82 * t249 + t84 * t250;
t27 = -t169 * t81 - t83 * t249 - t85 * t250;
t230 = -t238 / 0.2e1;
t229 = t238 / 0.2e1;
t228 = -t237 / 0.2e1;
t227 = t237 / 0.2e1;
t150 = rSges(2,1) * t184 - rSges(2,2) * t182;
t148 = rSges(2,1) * t182 + rSges(2,2) * t184;
t119 = rSges(3,1) * t173 + rSges(3,2) * t174;
t220 = t168 * t27 + t169 * t26;
t216 = t181 * t84 + t183 * t82;
t69 = t168 * t80;
t28 = -t169 * t216 + t69;
t29 = -t168 * t81 + t290;
t219 = t168 * t29 + t169 * t28;
t207 = t224 * qJD(1);
t198 = t152 + t207;
t32 = t179 * t226 + t109 + t198;
t110 = t149 * t237;
t195 = t206 - t153;
t33 = t179 * t287 - t110 + t195;
t218 = t168 * t32 - t169 * t33;
t217 = -t168 * t86 - t169 * t87;
t76 = t261 + t156 + (rSges(5,2) - pkin(3)) * t168;
t197 = t168 * t291 - t169 * t285 - t179 * t214;
t196 = t168 * t285 + t169 * t291 + t179 * t216;
t194 = -t211 * qJD(5) + t179 * t210;
t192 = -t105 + t198;
t73 = t207 - t247;
t61 = t168 * t277 + t156 + t87;
t10 = qJD(5) * t219 - t259;
t16 = -qJD(5) * t216 - t181 * t49 + t183 * t51;
t17 = qJD(5) * t214 - t181 * t48 + t183 * t50;
t20 = t194 * t168 + t169 * t282;
t21 = -t168 * t282 + t194 * t169;
t59 = t168 * t210 + t252;
t58 = t59 * t179;
t9 = qJD(5) * t220 + t58;
t190 = (t58 + ((-t28 + t69 + t27) * t168 + (t29 - t290 + (-t216 + t81) * t168 + t26) * t169) * qJD(5)) * t230 + t44 * t276 - t103 * t60 / 0.2e1 + (-qJD(5) * t210 + t124 * t181 - t125 * t183) * t179 + (-t215 + t59) * t275 + (t259 + (t168 ^ 2 * t81 + (-t69 + t27 + (t216 + t81) * t169) * t169) * qJD(5) + t10) * t228 + (t16 + t21) * t227 + (t17 + t20 + t9) * t229 + (Icges(4,3) + Icges(5,1) - t209) * t178;
t56 = t179 * t246 + t198;
t57 = t179 * t77 + t195;
t187 = t56 * (t137 + t153) + t57 * (t244 + t245) + (-t56 * t165 + (t56 * (-rSges(5,3) - qJ(4)) - t57 * pkin(3)) * t168) * t179;
t186 = t32 * (t153 - t208) + t33 * (-t168 * t232 + t231 + t245) + (t32 * t234 * t169 + (t32 * (-qJ(4) - t221) + t33 * t234) * t168) * t179;
t106 = t179 * t113;
t101 = t149 * t169;
t100 = t149 * t168;
t38 = qJD(5) * t217 + qJD(2);
t13 = -t103 * t86 - t104 * t87 + qJDD(2) + (-t168 * t53 + t169 * t52) * qJD(5);
t6 = t168 * t286 + t197 * t169;
t5 = -t168 * t281 + t196 * t169;
t4 = t197 * t168 - t169 * t286;
t3 = t196 * t168 + t169 * t281;
t1 = [t190 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t73 * t135 + (-t163 * t73 - t260) * t179 + (t224 * t74 - t289 * t73) * qJD(1) + t294 * (t117 + t289) + t295 * (-t114 + t224)) * m(4) + ((qJDD(1) * t120 - g(2) + t170) * t108 + (-qJDD(1) * t119 + t202 + (-0.2e1 * t120 + 0.2e1 * t108 - t176) * t185 - g(1)) * (-t119 - t272)) * m(3) + ((t148 ^ 2 + t150 ^ 2) * qJDD(1) + g(1) * t148 - g(2) * t150) * m(2) + ((t224 * t33 - t289 * t32) * qJD(1) + t186 - (-pkin(7) * t251 + t192 + t292 - t32) * t33 + t297 * (t287 + t289) + t299 * (t224 + t61)) * m(6) + ((t224 * t57 - t289 * t56) * qJD(1) + t187 - (t106 - t56 + t192) * t57 + t296 * (t77 + t289) + t298 * (t224 + t76)) * m(5); m(6) * t13 + t235 * qJDD(2) + (-m(6) - t235) * g(3); t190 + (-t32 * (t110 + t153) - t33 * (t241 + t292) - (-t271 * t33 - t287 * t32) * t179 + t186 + t297 * t287 + t299 * t61) * m(6) + (-t56 * t153 - t57 * (t106 + t241) + t187 + (t56 * t179 + t296) * t77 + t298 * t76) * m(5) + (t260 * t179 - t73 * t91 - t74 * t247 + (t73 * t179 + t294) * t117 - t295 * t114) * m(4); (-m(5) - m(6)) * (g(1) * t168 - g(2) * t169) + 0.2e1 * (t11 * t274 + t12 * t273) * m(6) + 0.2e1 * (t24 * t274 + t25 * t273) * m(5); -t9 * t251 / 0.2e1 + t169 * (t103 * t27 + t104 * t26 + t178 * t59 + t179 * t21 + (t168 * t6 + t169 * t5) * qJD(5)) / 0.2e1 + t220 * t275 + ((t179 * t27 + t5) * t169 + (-t179 * t26 + t6) * t168) * t227 + t10 * t248 / 0.2e1 + (t103 * t29 + t104 * t28 - t178 * t60 + t179 * t20 + (t168 * t4 + t169 * t3) * qJD(5)) * t274 + t219 * t276 + ((t179 * t29 + t3) * t169 + (-t179 * t28 + t4) * t168) * t229 + t178 * (t168 * t44 - t169 * t215) / 0.2e1 + t179 * ((t179 * t44 + t16) * t169 + (t179 * t215 + t17) * t168) / 0.2e1 + ((t94 * t237 - t291) * t169 + (-t293 + (t279 * t168 + (-t278 - t252) * t169) * qJD(5)) * t168) * t228 + ((-t238 * t252 - t291) * t168 + (t293 + (t278 * t169 + (-t279 + t94) * t168) * qJD(5)) * t169) * t230 - t179 * ((-t181 * t243 - t183 * t242) * t179 + ((t168 * t269 - t169 * t270) * t183 + (t168 * t267 - t169 * t268) * t181) * qJD(5)) / 0.2e1 + (t13 * t217 + t38 * ((-t179 * t86 + t52) * t169 + (-t53 + t75) * t168) - t218 * t126 + ((t179 * t32 - t12) * t169 + (t179 * t33 + t11) * t168) * t149 - (t100 * t33 + t101 * t32) * t179 - (t38 * (-t100 * t168 - t101 * t169) - t218 * t221) * qJD(5) - g(1) * t100 + g(2) * t101 + g(3) * t221) * m(6);];
tau = t1;
