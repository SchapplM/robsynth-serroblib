% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR6_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR6_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:42
% EndTime: 2019-12-31 18:17:47
% DurationCPUTime: 3.30s
% Computational Cost: add. (6562->288), mult. (4934->390), div. (0->0), fcn. (3618->8), ass. (0->171)
t162 = qJD(1) ^ 2;
t156 = qJD(1) + qJD(3);
t158 = sin(qJ(5));
t160 = cos(qJ(5));
t157 = qJ(1) + pkin(8);
t153 = qJ(3) + t157;
t149 = sin(t153);
t150 = cos(t153);
t243 = Icges(6,4) * t158;
t199 = Icges(6,2) * t160 + t243;
t73 = Icges(6,6) * t150 + t149 * t199;
t234 = t149 * t160;
t119 = Icges(6,4) * t234;
t235 = t149 * t158;
t75 = Icges(6,1) * t235 + Icges(6,5) * t150 + t119;
t203 = t158 * t75 + t160 * t73;
t74 = -Icges(6,6) * t149 + t150 * t199;
t242 = Icges(6,4) * t160;
t200 = Icges(6,1) * t158 + t242;
t76 = -Icges(6,5) * t149 + t150 * t200;
t42 = t158 * t74 - t160 * t76;
t183 = t199 * t156;
t125 = -Icges(6,2) * t158 + t242;
t271 = -Icges(6,6) * t156 + qJD(5) * t125;
t47 = t149 * t271 + t150 * t183;
t184 = t200 * t156;
t127 = Icges(6,1) * t160 - t243;
t270 = -Icges(6,5) * t156 + qJD(5) * t127;
t49 = t149 * t270 + t150 * t184;
t286 = -qJD(5) * t203 + t156 * t42 - t158 * t47 + t160 * t49;
t229 = t125 + t200;
t230 = -t199 + t127;
t285 = (t158 * t229 - t160 * t230) * t156;
t206 = rSges(6,1) * t158 + rSges(6,2) * t160;
t282 = 0.2e1 * qJD(5);
t78 = -t149 * rSges(6,3) + t150 * t206;
t69 = t156 * t78;
t249 = rSges(6,2) * t158;
t251 = rSges(6,1) * t160;
t130 = -t249 + t251;
t226 = qJD(5) * t149;
t93 = t130 * t226;
t281 = t69 + t93;
t146 = t150 * pkin(3);
t99 = t149 * qJ(4) + t146;
t194 = -rSges(5,2) * t150 + t149 * rSges(5,3) + t99;
t280 = t156 * t194;
t198 = Icges(6,5) * t158 + Icges(6,6) * t160;
t279 = t198 * t156;
t201 = t158 * t76 + t160 * t74;
t278 = t201 * t150;
t277 = pkin(2) * cos(t157) + cos(qJ(1)) * pkin(1);
t145 = t150 * pkin(7);
t142 = t150 * rSges(6,3);
t77 = rSges(6,1) * t235 + rSges(6,2) * t234 + t142;
t275 = t145 + t77 + t99;
t46 = t149 * t183 - t150 * t271;
t48 = t149 * t184 - t150 * t270;
t72 = -Icges(6,3) * t149 + t150 * t198;
t273 = qJD(5) * t42 + t156 * t72 + t158 * t48 + t160 * t46;
t123 = Icges(6,5) * t160 - Icges(6,6) * t158;
t272 = -Icges(6,3) * t156 + qJD(5) * t123;
t106 = t199 * qJD(5);
t107 = t200 * qJD(5);
t269 = (t125 * t158 - t127 * t160) * qJD(5) + t106 * t160 + t107 * t158 + t123 * t156;
t202 = t158 * t73 - t160 * t75;
t71 = Icges(6,3) * t150 + t149 * t198;
t268 = qJD(5) * t202 + t156 * t71 - t158 * t49 - t160 * t47;
t253 = t125 * t150 + t76;
t255 = -t127 * t150 + t74;
t266 = t158 * t255 - t160 * t253;
t254 = -Icges(6,2) * t235 + t119 + t75;
t256 = -t127 * t149 + t73;
t265 = t158 * t256 - t160 * t254;
t264 = -pkin(3) - pkin(7);
t263 = t149 / 0.2e1;
t262 = -t150 / 0.2e1;
t260 = -t156 / 0.2e1;
t258 = pkin(7) * t149;
t144 = t150 * rSges(4,1);
t101 = -rSges(4,2) * t149 + t144;
t188 = t277 * qJD(1);
t68 = t101 * t156 + t188;
t98 = rSges(4,1) * t149 + rSges(4,2) * t150;
t257 = t68 * t98;
t135 = t150 * qJ(4);
t96 = pkin(3) * t149 - t135;
t97 = t149 * rSges(5,2) + t150 * rSges(5,3);
t195 = -t96 + t97;
t246 = t156 * t98;
t197 = t160 * t125 + t158 * t127;
t83 = t149 * t123;
t56 = t150 * t197 - t83;
t245 = t56 * t156;
t132 = qJD(4) * t149;
t91 = t156 * t96;
t244 = t132 - t91;
t237 = t123 * t150;
t236 = t149 * t156;
t233 = t150 * t156;
t232 = qJ(4) * t233 + t132;
t231 = rSges(5,2) * t236 + rSges(5,3) * t233;
t227 = qJD(4) * t156;
t225 = qJD(5) * t150;
t224 = -rSges(6,3) + t264;
t24 = t150 * t71 + t73 * t234 + t75 * t235;
t25 = -t150 * t72 - t74 * t234 - t76 * t235;
t223 = qJD(5) * t251;
t222 = qJD(5) * t249;
t221 = t149 * t223 + t206 * t233;
t94 = t130 * t225;
t218 = -t226 / 0.2e1;
t216 = -t225 / 0.2e1;
t210 = -pkin(2) * sin(t157) - sin(qJ(1)) * pkin(1);
t189 = t210 * qJD(1);
t176 = t132 + t189;
t30 = t93 + (t78 - t96 - t258) * t156 + t176;
t133 = qJD(4) * t150;
t173 = t188 - t133;
t31 = t156 * t275 + t173 - t94;
t205 = t149 * t30 - t150 * t31;
t204 = -t149 * t77 - t150 * t78;
t193 = (t222 - t223) * t150;
t191 = t210 * t162;
t190 = t277 * t162;
t186 = (t149 * t25 + t150 * t24) * qJD(5);
t63 = t149 * t71;
t26 = -t150 * t203 + t63;
t27 = -t149 * t72 + t278;
t185 = (t149 * t27 + t150 * t26) * qJD(5);
t178 = t150 * t227 - t190;
t175 = t149 * t279 - t150 * t272 - t156 * t201;
t174 = t149 * t272 + t150 * t279 + t156 * t203;
t172 = -t198 * qJD(5) + t156 * t197;
t171 = t149 * t227 + t156 * (-pkin(3) * t236 + t232) + t191;
t170 = -t91 + t176;
t67 = t189 - t246;
t50 = (t149 * t206 + t142) * t156 + t193;
t51 = (-rSges(6,3) * t156 - t222) * t149 + t221;
t169 = (-t156 * t77 + t50) * t150 + (-t51 + t69) * t149;
t10 = t185 - t245;
t15 = qJD(5) * t201 - t158 * t46 + t160 * t48;
t20 = t172 * t149 + t150 * t269;
t21 = -t149 * t269 + t172 * t150;
t55 = t149 * t197 + t237;
t54 = t55 * t156;
t9 = t54 + t186;
t168 = (t54 + ((-t26 + t63 + t25) * t149 + (t27 - t278 + (-t203 + t72) * t149 + t24) * t150) * qJD(5)) * t218 + (-qJD(5) * t197 + t106 * t158 - t107 * t160) * t156 + (t245 + (t149 ^ 2 * t72 + (-t63 + t25 + (t203 + t72) * t150) * t150) * qJD(5) + t10) * t216 + (t15 + t20 + t9) * t226 / 0.2e1 + (t21 + t286) * t225 / 0.2e1 + (t150 * t56 + (-t202 + t55) * t149) * qJD(5) * t260;
t165 = t149 * t264 + t135 + t78;
t118 = rSges(5,2) * t233;
t52 = t156 * t195 + t176;
t53 = t173 + t280;
t164 = t52 * (t118 + t133) + t53 * (t231 + t232) + (-t52 * t146 + (t52 * (-rSges(5,3) - qJ(4)) - t53 * pkin(3)) * t149) * t156;
t163 = t30 * (t133 - t193) + t31 * (-t149 * t222 + t221 + t232) + (t30 * t224 * t150 + (t30 * (-qJ(4) - t206) + t31 * t224) * t149) * t156;
t155 = t156 ^ 2;
t116 = rSges(4,2) * t236;
t108 = t206 * qJD(5);
t92 = t156 * t97;
t90 = t130 * t150;
t89 = t130 * t149;
t80 = rSges(4,1) * t233 - t116;
t70 = t156 * t99 - t133;
t58 = -t156 * t80 - t190;
t57 = -t156 * t246 + t191;
t36 = qJD(5) * t204 + qJD(2);
t33 = (-rSges(5,3) * t236 + t118 - t70) * t156 + t178;
t32 = t156 * t231 + t171;
t17 = -t155 * t145 - t108 * t226 + (-t50 - t70 + t94) * t156 + t178;
t16 = -t155 * t258 + t156 * t51 + (t108 * t150 + t130 * t236) * qJD(5) + t171;
t11 = t169 * qJD(5);
t1 = [m(4) * (t58 * (t210 - t98) + t67 * t116 + t57 * (t101 + t277) + (-t144 * t67 - t257) * t156 + (t210 * t68 - t277 * t67) * qJD(1)) + t168 + (-(-pkin(7) * t236 + t170 + t281 - t30) * t31 + t17 * (t165 + t210) + t16 * (t275 + t277) + (t210 * t31 - t277 * t30) * qJD(1) + t163) * m(6) + (-(-t52 + t92 + t170) * t53 + t33 * (t195 + t210) + t32 * (t194 + t277) + (t210 * t53 - t277 * t52) * qJD(1) + t164) * m(5); m(6) * t11; t168 + (-t30 * (t133 + t94) - t31 * (t244 + t281) - (-t258 * t31 - t275 * t30) * t156 + t16 * t275 + t165 * t17 + t163) * m(6) + (t194 * t32 + t195 * t33 + t164 - t52 * (t133 - t280) - t53 * (t92 + t244)) * m(5) + (-(-t101 * t67 - t257) * t156 + t101 * t57 - t58 * t98 - t67 * t80 - t68 * t246) * m(4); 0.2e1 * (t16 * t262 + t17 * t263) * m(6) + 0.2e1 * (t262 * t32 + t263 * t33) * m(5); t156 * (t286 * t150 + (t156 * t202 + t15) * t149) / 0.2e1 + ((t83 * t225 - t279) * t150 + (-t285 + (t266 * t149 + (-t265 - t237) * t150) * qJD(5)) * t149) * t216 + ((-t226 * t237 - t279) * t149 + (t285 + (t265 * t150 + (-t266 + t83) * t149) * qJD(5)) * t150) * t218 + ((-t158 * t230 - t160 * t229) * t156 + ((t149 * t255 - t150 * t256) * t160 + (t149 * t253 - t150 * t254) * t158) * qJD(5)) * t260 + (t156 * t20 + ((t174 * t149 + t150 * t268 + t156 * t27) * t150 + (t175 * t149 - t150 * t273 - t156 * t26) * t149) * t282) * t263 + (t156 * t21 + ((-t149 * t268 + t174 * t150 + t156 * t25) * t150 + (t149 * t273 + t175 * t150 - t156 * t24) * t149) * t282) * t150 / 0.2e1 - (t9 + t186) * t236 / 0.2e1 + (t10 + t185) * t233 / 0.2e1 + (t11 * t204 + t36 * t169 - t205 * t108 + ((t156 * t30 - t16) * t150 + (t156 * t31 + t17) * t149) * t130 - (t30 * t90 + t31 * t89) * t156 - (t36 * (-t149 * t89 - t150 * t90) - t205 * t206) * qJD(5)) * m(6);];
tauc = t1(:);
