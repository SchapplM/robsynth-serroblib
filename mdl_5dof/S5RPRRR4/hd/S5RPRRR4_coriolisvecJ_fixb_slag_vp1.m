% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:19
% EndTime: 2019-12-05 18:14:26
% DurationCPUTime: 3.50s
% Computational Cost: add. (9883->297), mult. (5619->383), div. (0->0), fcn. (4190->10), ass. (0->195)
t155 = qJD(1) ^ 2;
t151 = sin(qJ(5));
t153 = cos(qJ(5));
t126 = rSges(6,1) * t151 + rSges(6,2) * t153;
t301 = qJD(5) * t126;
t150 = qJ(1) + pkin(9);
t145 = qJ(3) + t150;
t139 = cos(t145);
t149 = qJD(1) + qJD(3);
t236 = t139 * t149;
t227 = pkin(3) * t236;
t143 = cos(t150);
t154 = cos(qJ(1));
t278 = pkin(1) * t154;
t205 = pkin(2) * t143 + t278;
t292 = t205 * qJD(1);
t163 = t292 + t227;
t144 = qJD(4) + t149;
t140 = qJ(4) + t145;
t135 = sin(t140);
t136 = cos(t140);
t91 = rSges(5,1) * t136 - t135 * rSges(5,2);
t81 = t144 * t91;
t52 = -t163 - t81;
t146 = Icges(6,4) * t153;
t194 = -Icges(6,2) * t151 + t146;
t295 = Icges(6,1) * t151 + t146;
t233 = t295 + t194;
t251 = Icges(6,4) * t151;
t122 = Icges(6,2) * t153 + t251;
t125 = Icges(6,1) * t153 - t251;
t234 = t122 - t125;
t302 = (t233 * t151 + t234 * t153) * t144;
t142 = sin(t150);
t276 = pkin(2) * t142;
t152 = sin(qJ(1));
t279 = pkin(1) * t152;
t206 = t276 + t279;
t186 = t206 * qJD(1);
t300 = 0.2e1 * qJD(5);
t231 = qJD(5) * t135;
t92 = t126 * t231;
t272 = pkin(4) * t136;
t93 = pkin(8) * t135 + t272;
t299 = -t144 * t93 + t92;
t138 = sin(t145);
t200 = rSges(4,1) * t138 + rSges(4,2) * t139;
t263 = rSges(4,1) * t139;
t96 = -t138 * rSges(4,2) + t263;
t256 = t149 * t96;
t64 = -t292 - t256;
t298 = t64 * t200;
t238 = t136 * t153;
t239 = t136 * t151;
t68 = Icges(6,4) * t238 - Icges(6,2) * t239 + Icges(6,6) * t135;
t115 = Icges(6,4) * t239;
t70 = Icges(6,1) * t238 + Icges(6,5) * t135 - t115;
t196 = t151 * t68 - t153 * t70;
t297 = t196 * t136;
t86 = t200 * t149;
t242 = t135 * t151;
t235 = rSges(6,2) * t242 + t136 * rSges(6,3);
t137 = t142 * rSges(3,2);
t264 = rSges(3,1) * t143;
t204 = -t264 - t278;
t296 = t137 + t204;
t237 = t138 * t149;
t228 = pkin(3) * t237;
t165 = t186 + t228;
t230 = qJD(5) * t136;
t271 = pkin(8) * t136;
t241 = t135 * t153;
t226 = rSges(6,1) * t241;
t191 = t226 - t235;
t61 = t144 * t191;
t190 = t126 * t230 + t61 - t144 * (-pkin(4) * t135 + t271);
t28 = t165 + t190;
t117 = rSges(6,2) * t239;
t225 = rSges(6,1) * t238;
t188 = -rSges(6,3) * t135 - t225;
t71 = -t117 - t188;
t29 = t92 + (-t71 - t93) * t144 - t163;
t293 = t135 * t29 + t136 * t28;
t120 = Icges(6,5) * t151 + Icges(6,6) * t153;
t289 = -Icges(6,3) * t144 + t120 * qJD(5);
t288 = -Icges(6,6) * t144 + t122 * qJD(5);
t108 = t194 * qJD(5);
t109 = t125 * qJD(5);
t287 = (t122 * t153 + t151 * t295) * qJD(5) + t108 * t151 - t109 * t153 - t120 * t144;
t286 = -Icges(6,5) * t144 + qJD(5) * t295;
t265 = -Icges(6,2) * t238 - t115 + t70;
t267 = t136 * t295 + t68;
t284 = t265 * t151 + t267 * t153;
t114 = Icges(6,4) * t242;
t69 = -Icges(6,1) * t241 + Icges(6,5) * t136 + t114;
t266 = Icges(6,2) * t241 + t114 + t69;
t67 = Icges(6,6) * t136 - t194 * t135;
t268 = -t135 * t295 + t67;
t283 = -t266 * t151 - t268 * t153;
t280 = -rSges(6,3) - pkin(8);
t277 = pkin(1) * t155;
t275 = pkin(3) * t138;
t274 = pkin(3) * t139;
t273 = pkin(3) * t149 ^ 2;
t121 = Icges(6,5) * t153 - Icges(6,6) * t151;
t65 = Icges(6,3) * t136 - t121 * t135;
t270 = t136 * t65 + t67 * t242;
t269 = t135 * t65 + t69 * t238;
t261 = rSges(6,1) * t153;
t260 = rSges(6,2) * t151;
t62 = t144 * t71;
t255 = t151 * t67;
t254 = t153 * t69;
t192 = t122 * t151 - t153 * t295;
t74 = t120 * t135;
t50 = -t192 * t136 + t74;
t253 = t50 * t144;
t246 = t120 * t136;
t245 = t121 * t144;
t82 = t126 * t135;
t244 = t126 * t136;
t243 = t135 * t144;
t240 = t136 * t144;
t72 = rSges(5,1) * t243 + rSges(5,2) * t240;
t141 = t152 * t277;
t232 = t155 * t276 + t141;
t229 = -t62 + t299;
t224 = t144 * t117 + t135 * t301;
t223 = t136 * t301 + t144 * t226;
t220 = t138 * t273 + t232;
t217 = -pkin(4) - t261;
t215 = -t231 / 0.2e1;
t213 = -t230 / 0.2e1;
t212 = t230 / 0.2e1;
t66 = Icges(6,5) * t238 - Icges(6,6) * t239 + Icges(6,3) * t135;
t211 = -t66 - t254;
t102 = pkin(4) * t243;
t210 = t102 + t223;
t73 = -rSges(5,1) * t240 + rSges(5,2) * t243;
t207 = t235 + t271;
t201 = rSges(3,1) * t142 + rSges(3,2) * t143;
t199 = -rSges(5,1) * t135 - rSges(5,2) * t136;
t198 = -t260 + t261;
t41 = t151 * t69 + t153 * t67;
t197 = -t254 + t255;
t42 = t151 * t70 + t153 * t68;
t25 = t136 * t66 - t70 * t241 + t68 * t242;
t189 = -t91 - t274;
t187 = t205 * t155;
t24 = -t69 * t241 + t270;
t183 = (t135 * t25 + t136 * t24) * qJD(5);
t26 = -t67 * t239 + t269;
t27 = t135 * t66 - t297;
t182 = (t135 * t27 + t136 * t26) * qJD(5);
t181 = t207 - t275;
t180 = t125 * t144;
t179 = t194 * t144;
t177 = t117 - t225 - t272;
t175 = t199 - t275;
t171 = -t245 * t135 - t136 * t289 + t196 * t144;
t170 = t135 * t289 - t245 * t136 + t197 * t144;
t169 = t121 * qJD(5) + t192 * t144;
t168 = t177 - t274;
t167 = t135 * t191 + t136 * t71;
t166 = -t139 * t273 - t187;
t10 = t182 + t253;
t14 = -t197 * qJD(5) + t151 * (t135 * t286 - t136 * t180) + t153 * (t135 * t288 - t136 * t179);
t15 = -t196 * qJD(5) + t151 * (-t135 * t180 - t136 * t286) + t153 * (-t135 * t179 - t136 * t288);
t18 = t169 * t135 - t136 * t287;
t19 = t135 * t287 + t169 * t136;
t49 = t192 * t135 + t246;
t46 = t49 * t144;
t9 = t46 + t183;
t162 = (t46 + ((t27 + t270 + t297) * t136 + (-t26 + (t211 - t255) * t136 + t25 + t269) * t135) * qJD(5)) * t215 + (-t253 + ((t25 + (-t66 + t255) * t136 - t269) * t136 + (t211 * t135 - t24 + t270) * t135) * qJD(5) + t10) * t213 + (t14 + t19) * t212 + (t15 + t18 + t9) * t231 / 0.2e1 + (-t192 * qJD(5) + t108 * t153 + t109 * t151 + (t41 + t49) * t215 + (t42 + t50) * t212) * t144;
t43 = t144 * t235 - t223;
t44 = t188 * t144 + t224;
t159 = (-t44 - t62) * t135 + (t43 + t61) * t136;
t80 = t144 * t199;
t51 = -t80 + t165;
t158 = t52 * (t228 + t72) - t51 * (t73 - t227);
t110 = t198 * qJD(5);
t20 = t110 * t231 + (t102 - t43 + (-pkin(8) * t144 + t301) * t136) * t144 + t220;
t21 = -t110 * t230 + (t44 + t299) * t144 + t166;
t157 = (t20 * t280 + t21 * t217) * t135 + ((-t29 * t260 - t28 * t280) * t135 + (-t28 * t217 + t29 * t280) * t136) * t144;
t156 = t29 * (t210 + t228) - t28 * (t224 - t227) + t157;
t118 = rSges(4,2) * t237;
t87 = -rSges(4,1) * t236 + t118;
t63 = t186 + t86;
t57 = t149 * t87 - t187;
t56 = t149 * t86 + t232;
t48 = t144 * t73 + t166;
t47 = t144 * t72 + t220;
t32 = t167 * qJD(5) + qJD(2);
t11 = t159 * qJD(5);
t1 = [m(3) * ((t201 * t155 + t141) * t296 + (t154 * t277 + (-0.2e1 * t137 - t204 + t264 + t296) * t155) * (t201 + t279)) + t162 + (t20 * (t168 - t205) + t21 * (t181 - t206) + (t28 * t205 + t29 * t206) * qJD(1) + t156 - (t29 + t163 - t229) * t28) * m(6) + m(5) * (t47 * (t189 - t205) + t48 * (t175 - t206) + (t51 * t205 + t52 * t206) * qJD(1) + t158) + (t56 * (-t205 - t96) + t57 * (-t200 - t206) + t298 * t149 + t64 * t186 + (t149 * t263 - t118 - t256 - t64) * t63) * m(4); m(6) * t11; t162 + (t20 * t168 + t21 * t181 + t156 - t29 * (t190 + t228) + t28 * (-t227 + t229)) * m(6) + (-t52 * (-t80 + t228) + t51 * (-t81 - t227) + t48 * t175 + t47 * t189 + t158) * m(5) + (-t57 * t200 - t56 * t96 - t63 * t87 + t64 * t86 - (t63 * t96 + t298) * t149) * m(4); t162 + (t20 * t177 + t21 * t207 + t157 + (-t190 + t210) * t29 + (-t224 + t229) * t28) * m(6) + (t48 * t199 - t47 * t91 - t51 * t73 + t52 * t72 - (-t52 * t199 + t51 * t91) * t144) * m(5); t144 * ((t144 * t42 + t14) * t136 + (-t144 * t41 + t15) * t135) / 0.2e1 - t144 * ((-t234 * t151 + t233 * t153) * t144 + ((t265 * t135 + t266 * t136) * t153 + (-t267 * t135 - t268 * t136) * t151) * qJD(5)) / 0.2e1 + ((t74 * t230 + t245) * t136 + (t302 + (t284 * t135 + (-t283 - t246) * t136) * qJD(5)) * t135) * t213 + ((-t231 * t246 + t245) * t135 + (-t302 + (t283 * t136 + (-t284 + t74) * t135) * qJD(5)) * t136) * t215 + (t144 * t18 + ((t170 * t135 + t144 * t27) * t136 + (t171 * t135 - t144 * t26) * t135) * t300) * t135 / 0.2e1 + (t144 * t19 + ((t170 * t136 + t144 * t25) * t136 + (t171 * t136 - t144 * t24) * t135) * t300) * t136 / 0.2e1 - (t9 + t183) * t243 / 0.2e1 + (t10 + t182) * t240 / 0.2e1 + (t11 * t167 + t32 * t159 + t20 * t82 - t21 * t244 - (t244 * t29 - t28 * t82) * t144 - (t32 * (-t135 * t82 - t136 * t244) + t293 * t198) * qJD(5) + (t240 * t29 - t243 * t28) * t126 + t293 * t110) * m(6);];
tauc = t1(:);
