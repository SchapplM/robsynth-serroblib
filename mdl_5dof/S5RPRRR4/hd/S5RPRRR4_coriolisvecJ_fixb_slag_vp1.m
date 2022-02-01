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
% m [6x1]
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:34:28
% EndTime: 2022-01-23 09:34:33
% DurationCPUTime: 3.23s
% Computational Cost: add. (9883->279), mult. (5619->369), div. (0->0), fcn. (4190->10), ass. (0->183)
t154 = qJD(1) ^ 2;
t149 = qJ(1) + pkin(9);
t144 = qJ(3) + t149;
t138 = sin(t144);
t148 = qJD(1) + qJD(3);
t236 = t138 * t148;
t226 = pkin(3) * t236;
t278 = pkin(2) * cos(t149) + cos(qJ(1)) * pkin(1);
t189 = t278 * qJD(1);
t139 = cos(t144);
t235 = t139 * t148;
t227 = pkin(3) * t235;
t168 = t189 + t227;
t143 = qJD(4) + t148;
t140 = qJ(4) + t144;
t134 = sin(t140);
t135 = cos(t140);
t88 = rSges(5,1) * t135 - rSges(5,2) * t134;
t253 = t143 * t88;
t52 = t168 + t253;
t87 = rSges(5,1) * t134 + rSges(5,2) * t135;
t79 = t143 * t87;
t284 = t52 * (-t79 - t226);
t282 = 2 * qJD(5);
t150 = sin(qJ(5));
t241 = t134 * t150;
t110 = rSges(6,2) * t241;
t234 = rSges(6,3) * t135 + t110;
t152 = cos(qJ(5));
t240 = t134 * t152;
t69 = rSges(6,1) * t240 - t234;
t60 = t143 * t69;
t130 = t135 * pkin(8);
t90 = pkin(4) * t134 - t130;
t281 = -t143 * t90 - t60;
t206 = -pkin(2) * sin(t149) - sin(qJ(1)) * pkin(1);
t190 = t206 * qJD(1);
t169 = t190 - t226;
t51 = t169 - t79;
t237 = t135 * t152;
t280 = rSges(6,1) * t237 + t134 * rSges(6,3);
t279 = -pkin(4) * t135 - t134 * pkin(8);
t145 = Icges(6,4) * t152;
t195 = -Icges(6,2) * t150 + t145;
t117 = Icges(6,1) * t150 + t145;
t238 = t135 * t150;
t225 = rSges(6,2) * t238;
t70 = -t225 + t280;
t179 = t70 - t279;
t256 = rSges(6,2) * t152;
t119 = rSges(6,1) * t150 + t256;
t230 = qJD(5) * t134;
t89 = t119 * t230;
t204 = -t143 * t179 + t89;
t114 = Icges(6,5) * t152 - Icges(6,6) * t150;
t113 = Icges(6,5) * t150 + Icges(6,6) * t152;
t172 = Icges(6,3) * t143 - qJD(5) * t113;
t183 = t195 * t135;
t66 = Icges(6,6) * t134 + t183;
t250 = t150 * t66;
t247 = Icges(6,4) * t150;
t118 = Icges(6,1) * t152 - t247;
t184 = t118 * t135;
t68 = Icges(6,5) * t134 + t184;
t197 = -t152 * t68 + t250;
t242 = t134 * t143;
t276 = -t114 * t242 + t135 * t172 + t143 * t197;
t182 = t114 * t135;
t65 = Icges(6,4) * t240 - Icges(6,2) * t241 - Icges(6,6) * t135;
t251 = t150 * t65;
t109 = Icges(6,4) * t241;
t67 = Icges(6,1) * t240 - Icges(6,5) * t135 - t109;
t198 = -t152 * t67 + t251;
t275 = t172 * t134 + (t182 + t198) * t143;
t115 = Icges(6,2) * t152 + t247;
t194 = t115 * t150 - t117 * t152;
t273 = qJD(5) * t114 + t194 * t143;
t63 = Icges(6,5) * t240 - Icges(6,6) * t241 - Icges(6,3) * t135;
t24 = -t134 * t198 - t135 * t63;
t260 = -Icges(6,2) * t240 - t109 + t67;
t262 = t117 * t134 + t65;
t272 = -t150 * t260 - t152 * t262;
t269 = t143 / 0.2e1;
t267 = pkin(3) * t138;
t266 = pkin(3) * t148 ^ 2;
t132 = t139 * rSges(4,1);
t95 = -rSges(4,2) * t138 + t132;
t62 = t148 * t95 + t189;
t94 = rSges(4,1) * t138 + rSges(4,2) * t139;
t265 = t62 * t94;
t264 = -t134 * t63 - t237 * t67;
t64 = Icges(6,3) * t134 + t182;
t263 = t134 * t64 + t237 * t68;
t261 = -t117 * t135 - t66;
t259 = -t115 * t135 + t68;
t258 = rSges(6,1) * t152;
t229 = qJD(5) * t135;
t222 = t119 * t229;
t161 = t169 - t222;
t28 = (-t69 - t90) * t143 + t161;
t255 = t135 * t28;
t252 = t148 * t94;
t244 = t113 * t135;
t49 = -t134 * t194 - t244;
t249 = t49 * t143;
t245 = t113 * t134;
t243 = t114 * t143;
t239 = t135 * t143;
t233 = -t115 + t118;
t232 = t117 + t195;
t228 = qJD(5) * t150;
t223 = qJD(5) * t256;
t224 = t143 * t225 + (rSges(6,1) * t228 + t223) * t134;
t221 = t135 * t228;
t218 = -pkin(4) - t258;
t217 = -t230 / 0.2e1;
t214 = t229 / 0.2e1;
t53 = t68 * t240;
t212 = t135 * t64 - t53;
t211 = -t63 + t250;
t72 = rSges(5,1) * t239 - rSges(5,2) * t242;
t133 = pkin(3) * t139;
t207 = t133 + t88;
t201 = -rSges(6,2) * t150 + t258;
t29 = t168 - t204;
t200 = -t134 * t29 - t255;
t199 = t134 * t69 + t135 * t70;
t41 = t150 * t67 + t152 * t65;
t42 = t150 * t68 + t152 * t66;
t193 = -t222 + t281;
t192 = t206 * t154;
t191 = t278 * t154;
t188 = rSges(6,3) * t239 + t110 * t143 - t135 * t223;
t25 = -t241 * t66 - t212;
t186 = (t134 * t25 - t135 * t24) * qJD(5);
t26 = -t238 * t65 - t264;
t27 = -t238 * t66 + t263;
t185 = (t134 * t27 - t135 * t26) * qJD(5);
t180 = -t87 - t267;
t178 = -t150 * t259 + t152 * t261;
t177 = t134 * t218 + t130 + t234;
t176 = t133 + t179;
t61 = t190 - t252;
t175 = (-t150 * t232 + t152 * t233) * t143;
t174 = Icges(6,5) * t143 - qJD(5) * t117;
t173 = Icges(6,6) * t143 - qJD(5) * t115;
t171 = -t138 * t266 + t192;
t170 = -t139 * t266 - t191;
t167 = t177 - t267;
t43 = (-t143 * t240 - t221) * rSges(6,1) + t188;
t44 = t143 * t280 - t224;
t166 = (t43 + t60) * t135 + (-t143 * t70 + t44) * t134;
t50 = -t135 * t194 + t245;
t46 = t50 * t143;
t10 = t46 + t185;
t103 = t195 * qJD(5);
t104 = t118 * qJD(5);
t14 = -t198 * qJD(5) + t150 * (t134 * t174 + t143 * t184) + t152 * (t134 * t173 + t143 * t183);
t15 = -t197 * qJD(5) + t150 * (-t118 * t242 + t135 * t174) + t152 * (t135 * t173 - t195 * t242);
t158 = -t103 * t150 + t104 * t152 + t113 * t143 + (-t115 * t152 - t117 * t150) * qJD(5);
t18 = t134 * t273 + t135 * t158;
t19 = t134 * t158 - t135 * t273;
t9 = t186 + t249;
t165 = (t46 + ((t25 - t53 + (t64 + t251) * t135 + t264) * t135 + t263 * t134) * qJD(5)) * t214 + (-qJD(5) * t194 + t103 * t152 + t104 * t150) * t143 + (t9 - t249 + ((t135 * t211 - t263 + t27) * t135 + (t134 * t211 + t212 + t26) * t134) * qJD(5)) * t217 + (t18 + t15) * t230 / 0.2e1 - (t10 + t14 + t19) * t229 / 0.2e1 + ((t49 + t41) * t134 + (t42 + t50) * t135) * qJD(5) * t269;
t99 = pkin(8) * t239;
t164 = -rSges(6,1) * t221 + t188 + t99;
t157 = (t218 * t255 + (t28 * (-rSges(6,3) - pkin(8)) + t29 * t218) * t134) * t143;
t156 = t51 * (-t72 - t227) + t284;
t155 = t28 * (t224 - t227) + t29 * (t164 - t226) + t157;
t112 = rSges(4,2) * t236;
t105 = t201 * qJD(5);
t84 = rSges(4,1) * t235 - t112;
t81 = t119 * t135;
t80 = t119 * t134;
t57 = -t148 * t84 - t191;
t56 = -t148 * t252 + t192;
t48 = -t143 * t72 + t170;
t47 = -t143 * t79 + t171;
t32 = qJD(5) * t199 + qJD(2);
t21 = -t105 * t229 + (t143 * t279 - t44 + t89) * t143 + t170;
t20 = -t105 * t230 + (-pkin(4) * t242 - t222 + t43 + t99) * t143 + t171;
t11 = t166 * qJD(5);
t1 = [m(4) * (t57 * (t206 - t94) + t61 * t112 + t56 * (t95 + t278) + (-t132 * t61 - t265) * t148 + (t206 * t62 - t278 * t61) * qJD(1)) + t165 + m(5) * (t48 * (t180 + t206) + t47 * (t207 + t278) + (t206 * t52 - t278 * t51) * qJD(1) + t156) + (-(-t28 + t161 + t281) * t29 + t21 * (t167 + t206) + t20 * (t176 + t278) + (t206 * t29 - t278 * t28) * qJD(1) + t155) * m(6); m(6) * t11; t165 + (-t28 * (t204 - t227) - t29 * (t193 - t226) + t21 * t167 + t20 * t176 + t155) * m(6) + (-t51 * (-t227 - t253) - t284 + t48 * t180 + t47 * t207 + t156) * m(5) + (-(-t61 * t95 - t265) * t148 + t56 * t95 - t57 * t94 - t61 * t84 - t62 * t252) * m(4); t165 + (t21 * t177 + t20 * t179 + t157 + (t164 - t193) * t29 + (-t204 + t224) * t28) * m(6) + (-(-t51 * t88 - t52 * t87) * t143 + t47 * t88 - t48 * t87 - t51 * t72 - t52 * t79) * m(5); ((t143 * t42 - t14) * t135 + (t143 * t41 + t15) * t134) * t269 + ((-t230 * t244 + t243) * t134 + (t175 + (-t272 * t135 + (t245 + t178) * t134) * qJD(5)) * t135) * t217 + ((-t229 * t245 - t243) * t135 + (t175 + (t178 * t134 + (-t272 + t244) * t135) * qJD(5)) * t134) * t214 - t143 * ((t233 * t150 + t232 * t152) * t143 + ((t134 * t259 - t135 * t260) * t152 + (t134 * t261 + t135 * t262) * t150) * qJD(5)) / 0.2e1 + (t143 * t18 + ((-t134 * t275 + t143 * t27) * t135 + (t134 * t276 + t143 * t26) * t134) * t282) * t134 / 0.2e1 - (t143 * t19 + ((t135 * t275 + t143 * t25) * t135 + (-t135 * t276 + t143 * t24) * t134) * t282) * t135 / 0.2e1 + (t9 + t186) * t242 / 0.2e1 + (t10 + t185) * t239 / 0.2e1 + (t11 * t199 + t32 * t166 + t200 * t105 + ((-t143 * t29 - t21) * t135 + (t143 * t28 - t20) * t134) * t119 - (t28 * t80 - t29 * t81) * t143 - (t32 * (-t134 * t80 - t135 * t81) + t200 * t201) * qJD(5)) * m(6);];
tauc = t1(:);
