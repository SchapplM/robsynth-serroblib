% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR1
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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:55:43
% EndTime: 2020-01-03 11:55:50
% DurationCPUTime: 4.11s
% Computational Cost: add. (8624->348), mult. (5274->431), div. (0->0), fcn. (3876->10), ass. (0->217)
t171 = qJ(1) + qJ(2);
t162 = pkin(8) + t171;
t153 = sin(t162);
t154 = cos(t162);
t169 = pkin(9) + qJ(5);
t160 = sin(t169);
t161 = cos(t169);
t170 = qJD(1) + qJD(2);
t150 = Icges(6,4) * t161;
t207 = -Icges(6,2) * t160 + t150;
t196 = t207 * t170;
t284 = Icges(6,4) * t160;
t112 = Icges(6,1) * t161 - t284;
t197 = t112 * t170;
t71 = -Icges(6,5) * t154 + t112 * t153;
t289 = t161 * t71;
t69 = -Icges(6,6) * t154 + t153 * t207;
t291 = t160 * t69;
t213 = t289 - t291;
t318 = Icges(6,1) * t160 + t150;
t308 = -Icges(6,5) * t170 + qJD(5) * t318;
t109 = Icges(6,2) * t161 + t284;
t311 = -Icges(6,6) * t170 + qJD(5) * t109;
t273 = t154 * t161;
t274 = t154 * t160;
t70 = Icges(6,4) * t273 - Icges(6,2) * t274 + Icges(6,6) * t153;
t121 = Icges(6,4) * t274;
t72 = Icges(6,1) * t273 + Icges(6,5) * t153 - t121;
t326 = -qJD(5) * t213 - t160 * (-t153 * t308 + t154 * t197) - t161 * (-t153 * t311 + t154 * t196) + t170 * (t160 * t72 + t161 * t70);
t101 = t153 * rSges(4,1) + t154 * rSges(4,2);
t163 = sin(t171);
t156 = pkin(2) * t163;
t316 = t156 + t101;
t279 = qJ(4) * t153;
t302 = pkin(3) * t154;
t102 = t279 + t302;
t164 = cos(t171);
t269 = t164 * t170;
t140 = pkin(2) * t269;
t325 = -t170 * t102 + t140;
t175 = sin(qJ(1));
t292 = pkin(1) * qJD(1);
t245 = t175 * t292;
t317 = t163 * rSges(3,1) + t164 * rSges(3,2);
t94 = t317 * t170;
t324 = -t245 - t94;
t264 = t318 + t207;
t265 = t109 - t112;
t322 = (t160 * t264 + t161 * t265) * t170;
t321 = 0.2e1 * qJD(5);
t272 = t154 * t170;
t134 = rSges(4,1) * t272;
t295 = rSges(4,2) * t153;
t103 = t154 * rSges(4,1) - t295;
t93 = t170 * t103;
t320 = t134 - t93;
t319 = t170 * t316;
t173 = cos(pkin(9));
t155 = pkin(4) * t173 + pkin(3);
t174 = -pkin(7) - qJ(4);
t261 = t153 * t155 + t154 * t174;
t104 = t155 * t272;
t115 = rSges(6,1) * t160 + rSges(6,2) * t161;
t255 = qJD(4) * t154;
t227 = t140 - t255;
t253 = qJD(5) * t153;
t192 = -t115 * t253 + t227;
t248 = rSges(6,1) * t273;
t275 = t153 * t170;
t266 = rSges(6,3) * t275 + t170 * t248;
t123 = t154 * t155;
t62 = t302 - t123 + (qJ(4) + t174) * t153;
t125 = rSges(6,2) * t274;
t75 = rSges(6,3) * t153 - t125 + t248;
t63 = t170 * t75;
t315 = t170 * t62 + t104 - t192 + t266 + t325 - t63;
t260 = pkin(3) * t272 + qJ(4) * t275;
t297 = rSges(5,1) * t173;
t137 = t154 * t297;
t262 = rSges(5,3) * t275 + t170 * t137;
t258 = t153 * rSges(5,3) + t137;
t294 = rSges(5,2) * sin(pkin(9));
t80 = -t154 * t294 + t258;
t314 = -t170 * t80 - t227 + t260 + t262 + t325;
t107 = Icges(6,5) * t160 + Icges(6,6) * t161;
t312 = -Icges(6,3) * t170 + qJD(5) * t107;
t97 = t207 * qJD(5);
t98 = t112 * qJD(5);
t310 = qJD(5) * (t109 * t161 + t160 * t318) - t107 * t170 + t160 * t97 - t161 * t98;
t168 = t170 ^ 2;
t304 = t170 / 0.2e1;
t303 = pkin(2) * t168;
t149 = t153 * pkin(3);
t166 = t175 * pkin(1);
t176 = cos(qJ(1));
t167 = t176 * pkin(1);
t301 = t153 * t318 + t69;
t300 = t154 * t318 + t70;
t299 = -t109 * t153 + t71;
t298 = -Icges(6,2) * t273 - t121 + t72;
t296 = rSges(6,1) * t161;
t293 = rSges(6,2) * t160;
t290 = t160 * t70;
t159 = t176 * t292;
t26 = t159 + (t102 - t62 + t75) * t170 + t192;
t288 = t170 * t26;
t43 = t159 + (t102 + t80) * t170 + t227;
t287 = t170 * t43;
t278 = t107 * t153;
t46 = -t109 * t274 + t273 * t318 + t278;
t286 = t46 * t170;
t82 = t107 * t154;
t108 = Icges(6,5) * t161 - Icges(6,6) * t160;
t195 = t108 * t170;
t277 = t153 * t160;
t276 = t153 * t161;
t271 = t160 * t170;
t270 = t163 * t170;
t268 = t170 * t174;
t246 = rSges(6,2) * t271;
t267 = -rSges(6,3) * t272 - t153 * t246;
t247 = t170 * t294;
t263 = rSges(5,3) * t272 + t153 * t247;
t129 = qJ(4) * t272;
t144 = qJD(4) * t153;
t259 = t129 + t144;
t177 = qJD(1) ^ 2;
t158 = t177 * t167;
t257 = t164 * t303 + t158;
t254 = qJD(5) * t115;
t252 = qJD(5) * t154;
t251 = qJD(5) * t160;
t250 = qJD(5) * t161;
t249 = t177 * t166;
t136 = t153 * t297;
t244 = t170 * (-t255 + t260) + t257;
t243 = rSges(6,1) * t251;
t242 = rSges(6,2) * t250;
t157 = pkin(2) * t164;
t241 = t123 - t125 + t157;
t239 = t115 * t252;
t235 = t253 / 0.2e1;
t233 = t252 / 0.2e1;
t100 = -qJ(4) * t154 + t149;
t231 = t100 + t156;
t53 = t154 * t242 + (t154 * t251 + t161 * t275) * rSges(6,1) + t267;
t124 = rSges(6,1) * t276;
t74 = -rSges(6,2) * t277 - rSges(6,3) * t154 + t124;
t229 = t170 * t74 - t53;
t54 = -t153 * t243 + (-t153 * t250 - t154 * t271) * rSges(6,2) + t266;
t228 = t54 - t63;
t118 = rSges(3,1) * t164 - rSges(3,2) * t163;
t91 = t118 * t170 + t159;
t225 = t124 + t156 + t261;
t221 = -t153 * t294 + t136;
t224 = -rSges(5,3) * t154 + t221 + t231;
t95 = rSges(3,1) * t269 - rSges(3,2) * t270;
t220 = -t144 + t245;
t219 = t103 + t157;
t218 = -qJD(4) - t247;
t217 = t144 - t239;
t214 = -t293 + t296;
t40 = t160 * t71 + t161 * t69;
t211 = -t161 * t72 + t290;
t209 = t231 - t100 + t261 + t74;
t205 = -t109 * t160 + t161 * t318;
t203 = t157 + t258 + t279;
t202 = -pkin(2) * t270 + t144 - t267;
t201 = t149 + t156 + t221;
t55 = t71 * t276;
t67 = -Icges(6,3) * t154 + t108 * t153;
t27 = -t154 * t67 - t277 * t69 + t55;
t56 = t72 * t276;
t68 = Icges(6,5) * t273 - Icges(6,6) * t274 + Icges(6,3) * t153;
t28 = t154 * t68 + t277 * t70 - t56;
t200 = (-t153 * t28 - t154 * t27) * qJD(5);
t57 = t69 * t274;
t29 = -t153 * t67 - t273 * t71 + t57;
t30 = t153 * t68 - t211 * t154;
t199 = (-t153 * t30 - t154 * t29) * qJD(5);
t198 = -t163 * t303 - t249;
t190 = t170 * t144 + t198;
t189 = t160 * t299 + t161 * t301;
t188 = t160 * t298 + t161 * t300;
t187 = -t153 * t195 - t154 * t312 + t170 * t211;
t186 = t153 * t312 - t154 * t195 + t170 * t213;
t185 = -t108 * qJD(5) + t170 * t205;
t10 = t199 - t286;
t17 = qJD(5) * t211 + t160 * (t153 * t197 + t154 * t308) + t161 * (t153 * t196 + t154 * t311);
t20 = t185 * t153 + t154 * t310;
t21 = -t153 * t310 + t185 * t154;
t45 = t153 * t205 - t82;
t44 = t45 * t170;
t9 = t44 + t200;
t183 = (qJD(5) * t205 + t160 * t98 + t161 * t97) * t170 + (t44 + ((t29 + t56 - t57 + (t67 - t290) * t153) * t153 + (-t55 - t30 + (-t211 + t67) * t154 + (t289 + t291) * t153) * t154) * qJD(5)) * t235 + (t10 + t286 + ((-t28 + t57 + (t68 - t289) * t154) * t154 + (t27 - t55 + (t68 + t291) * t153) * t153) * qJD(5)) * t233 - (t17 + t20 + t9) * t253 / 0.2e1 - (t21 - t326) * t252 / 0.2e1 + (t154 * t46 + (t40 + t45) * t153) * qJD(5) * t304;
t64 = t245 + t319;
t65 = t159 + t140 + t93;
t180 = (-t64 * t295 - t65 * t316) * t170;
t76 = pkin(3) * t275 - t259;
t31 = (-t136 * t170 + t263 - t76) * t170 + t190;
t32 = (t154 * t218 + t262) * t170 + t244;
t42 = t170 * t224 + t220;
t179 = (t31 * (pkin(3) - t294) + t32 * (-rSges(5,3) - qJ(4)) + t42 * t218) * t154 + (-t156 + (-pkin(3) - t297) * t153) * t287;
t99 = t214 * qJD(5);
t14 = -t99 * t253 + (-t239 - t129 - t53 - t76 + (t149 - t261) * t170) * t170 + t190;
t15 = t99 * t252 + (-t255 + t104 + t54 + (-t254 - t268) * t153 - t260) * t170 + t244;
t25 = t170 * t209 - t217 + t245;
t178 = (t14 * (rSges(6,3) - t174) - t15 * t293 - t25 * t254 + (t26 * (-t155 - t296) - t25 * t174) * t170) * t153 + (t14 * t296 + t26 * (-t242 - t243 - t268) - t15 * rSges(6,3) + t25 * (-qJD(4) - t246)) * t154;
t88 = t115 * t154;
t87 = t115 * t153;
t78 = t170 * t95 + t158;
t77 = -t170 * t94 - t249;
t59 = t170 * (-rSges(4,2) * t275 + t134) + t257;
t58 = -t101 * t168 + t198;
t39 = qJD(3) + (t153 * t74 + t154 * t75) * qJD(5);
t11 = (t228 * t153 + t229 * t154) * qJD(5);
t1 = [t183 + m(3) * (t77 * (t118 + t167) + t78 * (t166 + t317) + (t91 - t159 - t95) * t324) + (t14 * (t167 + t241) + t26 * (t202 - t245) + t15 * (t166 + t225) + t178 + (t26 + t315) * t25) * m(6) + (t31 * (t167 + t203) + t43 * (t129 - t220 + t263) + t32 * (t166 + t201) + t179 + (t43 + t314) * t42) * m(5) + (t58 * (t167 + t219) - t65 * t245 + t59 * (t166 + t316) + t180 + (t65 + t320) * t64) * m(4); t183 + (t14 * t241 + t15 * t225 + t209 * t288 + t178 + (t202 - t217) * t26 + t315 * t25) * m(6) + (t32 * t201 + t31 * t203 + t224 * t287 + t179 + (-t144 + t259 + t263) * t43 + t314 * t42) * m(5) + (t58 * t219 + t316 * t59 + t319 * t65 + t320 * t64 + t180) * m(4) + (-(-t118 * t324 - t317 * t91) * t170 + t317 * t78 + t118 * t77 - t324 * t95 - t91 * t94) * m(3); m(6) * t11; m(5) * (-t153 * t32 - t31 * t154) + m(6) * (-t14 * t154 - t15 * t153); (t326 * t154 + (t170 * t40 - t17) * t153) * t304 - t170 * ((-t160 * t265 + t161 * t264) * t170 + ((t153 * t298 - t154 * t299) * t161 + (-t153 * t300 + t154 * t301) * t160) * qJD(5)) / 0.2e1 + ((-t252 * t278 - t195) * t154 + (-t322 + (-t188 * t153 + (t82 + t189) * t154) * qJD(5)) * t153) * t233 + ((t82 * t253 - t195) * t153 + (t322 + (-t189 * t154 + (-t278 + t188) * t153) * qJD(5)) * t154) * t235 - (t170 * t20 + ((-t186 * t153 - t170 * t30) * t154 + (-t187 * t153 + t170 * t29) * t153) * t321) * t153 / 0.2e1 - (t170 * t21 + ((-t186 * t154 - t170 * t28) * t154 + (-t187 * t154 + t170 * t27) * t153) * t321) * t154 / 0.2e1 + (t9 + t200) * t275 / 0.2e1 - (t10 + t199) * t272 / 0.2e1 + ((t11 * t75 + t39 * t229 + t25 * t99) * t154 + (t11 * t74 + t39 * t228 - t26 * t99) * t153 + ((t15 - t288) * t154 + (-t170 * t25 - t14) * t153) * t115 - (-t25 * t87 - t26 * t88) * t170 - (t39 * (-t153 * t87 - t154 * t88) + (-t153 * t26 + t154 * t25) * t214) * qJD(5)) * m(6);];
tauc = t1(:);
