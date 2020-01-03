% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:08
% EndTime: 2019-12-31 17:22:12
% DurationCPUTime: 2.82s
% Computational Cost: add. (7351->265), mult. (5325->352), div. (0->0), fcn. (4048->8), ass. (0->177)
t142 = sin(qJ(1));
t241 = pkin(1) * qJD(1);
t208 = t142 * t241;
t140 = qJ(1) + qJ(2);
t133 = sin(t140);
t139 = qJD(1) + qJD(2);
t221 = t133 * t139;
t213 = pkin(2) * t221;
t132 = qJD(3) + t139;
t136 = qJ(3) + t140;
t130 = sin(t136);
t131 = cos(t136);
t85 = rSges(4,1) * t130 + rSges(4,2) * t131;
t78 = t132 * t85;
t53 = -t208 - t213 - t78;
t268 = 2 * qJD(4);
t224 = t131 * t132;
t100 = pkin(7) * t224;
t141 = sin(qJ(4));
t226 = t130 * t141;
t106 = rSges(5,2) * t226;
t143 = cos(qJ(4));
t242 = rSges(5,2) * t143;
t205 = qJD(4) * t242;
t179 = rSges(5,3) * t224 + t132 * t106 - t131 * t205;
t214 = qJD(4) * t141;
t203 = t131 * t214;
t219 = t131 * rSges(5,3) + t106;
t225 = t130 * t143;
t66 = rSges(5,1) * t225 - t219;
t57 = t132 * t66;
t126 = t131 * pkin(7);
t87 = pkin(3) * t130 - t126;
t267 = -rSges(5,1) * t203 + t132 * t87 + t100 + t179 + t57;
t134 = cos(t140);
t91 = rSges(3,1) * t133 + rSges(3,2) * t134;
t237 = t139 * t91;
t70 = -t208 - t237;
t115 = rSges(5,1) * t141 + t242;
t215 = qJD(4) * t131;
t204 = t115 * t215;
t168 = -t204 - t213;
t156 = t168 - t208;
t27 = (-t66 - t87) * t132 + t156;
t240 = t131 * t27;
t144 = cos(qJ(1));
t207 = t144 * t241;
t220 = t134 * t139;
t212 = pkin(2) * t220;
t169 = t207 + t212;
t262 = -t131 * pkin(3) - t130 * pkin(7);
t223 = t131 * t141;
t209 = rSges(5,2) * t223;
t222 = t131 * t143;
t263 = rSges(5,1) * t222 + t130 * rSges(5,3);
t67 = -t209 + t263;
t166 = t67 - t262;
t216 = qJD(4) * t130;
t84 = t115 * t216;
t187 = -t132 * t166 + t84;
t28 = t169 - t187;
t266 = -t130 * t28 - t240;
t135 = Icges(5,4) * t143;
t181 = -Icges(5,2) * t141 + t135;
t112 = Icges(5,1) * t141 + t135;
t261 = -t213 + t267;
t109 = Icges(5,5) * t143 - Icges(5,6) * t141;
t232 = Icges(5,4) * t141;
t110 = Icges(5,2) * t143 + t232;
t180 = t141 * t110 - t143 * t112;
t260 = t109 * qJD(4) + t132 * t180;
t108 = Icges(5,5) * t141 + Icges(5,6) * t143;
t158 = Icges(5,3) * t132 - qJD(4) * t108;
t172 = t181 * t131;
t63 = Icges(5,6) * t130 + t172;
t235 = t141 * t63;
t113 = Icges(5,1) * t143 - t232;
t173 = t113 * t131;
t65 = Icges(5,5) * t130 + t173;
t183 = -t143 * t65 + t235;
t227 = t130 * t132;
t259 = -t109 * t227 + t131 * t158 + t132 * t183;
t171 = t109 * t131;
t62 = Icges(5,4) * t225 - Icges(5,2) * t226 - Icges(5,6) * t131;
t236 = t141 * t62;
t105 = Icges(5,4) * t226;
t64 = Icges(5,1) * t225 - Icges(5,5) * t131 - t105;
t184 = -t143 * t64 + t236;
t258 = t130 * t158 + (t171 + t184) * t132;
t60 = Icges(5,5) * t225 - Icges(5,6) * t226 - Icges(5,3) * t131;
t23 = -t130 * t184 - t131 * t60;
t243 = rSges(5,1) * t143;
t200 = -pkin(3) - t243;
t146 = (t200 * t240 + (t27 * (-rSges(5,3) - pkin(7)) + t28 * t200) * t130) * t132;
t206 = t132 * t209 + (rSges(5,1) * t214 + t205) * t130;
t257 = t146 + (t206 - t187) * t27;
t245 = -Icges(5,2) * t225 - t105 + t64;
t247 = t112 * t130 + t62;
t256 = -t141 * t245 - t143 * t247;
t253 = t132 / 0.2e1;
t252 = pkin(1) * t142;
t251 = pkin(2) * t133;
t250 = pkin(2) * t139 ^ 2;
t137 = t144 * pkin(1);
t249 = -t130 * t60 - t64 * t222;
t61 = Icges(5,3) * t130 + t171;
t248 = t130 * t61 + t65 * t222;
t246 = -t112 * t131 - t63;
t244 = -t110 * t131 + t65;
t86 = t131 * rSges(4,1) - rSges(4,2) * t130;
t238 = t132 * t86;
t229 = t108 * t131;
t46 = -t130 * t180 - t229;
t234 = t46 * t132;
t230 = t108 * t130;
t228 = t109 * t132;
t218 = -t110 + t113;
t217 = t112 + t181;
t145 = qJD(1) ^ 2;
t211 = t145 * t252;
t210 = t145 * t137;
t199 = -t216 / 0.2e1;
t196 = t215 / 0.2e1;
t50 = t65 * t225;
t194 = t131 * t61 - t50;
t191 = -t60 + t235;
t92 = t134 * rSges(3,1) - rSges(3,2) * t133;
t69 = rSges(4,1) * t224 - rSges(4,2) * t227;
t83 = rSges(3,1) * t220 - rSges(3,2) * t221;
t129 = pkin(2) * t134;
t188 = t129 + t86;
t185 = -rSges(5,2) * t141 + t243;
t43 = t141 * t64 + t143 * t62;
t44 = t141 * t65 + t143 * t63;
t24 = -t226 * t63 - t194;
t177 = (t130 * t24 - t131 * t23) * qJD(4);
t25 = -t223 * t62 - t249;
t26 = -t223 * t63 + t248;
t176 = (t130 * t26 - t131 * t25) * qJD(4);
t31 = (t130 * t66 + t131 * t67) * qJD(4);
t175 = -t133 * t250 - t211;
t174 = -t134 * t250 - t210;
t167 = -t85 - t251;
t165 = -t69 - t212;
t164 = -t141 * t244 + t143 * t246;
t163 = t130 * t200 + t126 + t219;
t162 = t129 + t166;
t161 = (-t141 * t217 + t143 * t218) * t132;
t160 = Icges(5,5) * t132 - qJD(4) * t112;
t159 = Icges(5,6) * t132 - qJD(4) * t110;
t155 = t163 - t251;
t47 = -t131 * t180 + t230;
t45 = t47 * t132;
t10 = t45 + t176;
t13 = -qJD(4) * t184 + t141 * (t130 * t160 + t132 * t173) + t143 * (t130 * t159 + t132 * t172);
t14 = -qJD(4) * t183 + t141 * (-t113 * t227 + t131 * t160) + t143 * (t131 * t159 - t181 * t227);
t94 = t181 * qJD(4);
t95 = t113 * qJD(4);
t148 = t108 * t132 - t141 * t94 + t143 * t95 + (-t110 * t143 - t112 * t141) * qJD(4);
t17 = t130 * t260 + t148 * t131;
t18 = t148 * t130 - t131 * t260;
t9 = t177 + t234;
t154 = (t45 + ((t24 - t50 + (t61 + t236) * t131 + t249) * t131 + t248 * t130) * qJD(4)) * t196 + (-qJD(4) * t180 + t141 * t95 + t143 * t94) * t132 + (-t234 + ((t131 * t191 - t248 + t26) * t131 + (t130 * t191 + t194 + t25) * t130) * qJD(4) + t9) * t199 + (t14 + t17) * t216 / 0.2e1 - (t13 + t18 + t10) * t215 / 0.2e1 + ((t43 + t46) * t130 + (t44 + t47) * t131) * qJD(4) * t253;
t98 = t185 * qJD(4);
t80 = t115 * t131;
t79 = t115 * t130;
t71 = t139 * t92 + t207;
t59 = -t139 * t83 - t210;
t58 = -t139 * t237 - t211;
t54 = t169 + t238;
t49 = -t132 * t69 + t174;
t48 = -t132 * t78 + t175;
t42 = t132 * t263 - t206;
t41 = (-t132 * t225 - t203) * rSges(5,1) + t179;
t20 = -t98 * t215 + (t132 * t262 - t42 + t84) * t132 + t174;
t19 = -t98 * t216 + (-pkin(3) * t227 + t100 - t204 + t41) * t132 + t175;
t1 = [t154 + m(3) * (t59 * (-t91 - t252) + t58 * (t137 + t92) + (-t83 - t207 + t71) * t70) + (t20 * (t155 - t252) + t27 * (-t169 + t206) + t19 * (t137 + t162) + t146 + (-t208 - t156 + t27 + t261) * t28) * m(5) + m(4) * (t49 * (t167 - t252) + t48 * (t137 + t188) + (t165 - t207 + t54) * t53); t154 + (t20 * t155 + t19 * t162 + (-t168 + t261) * t28 + t257) * m(5) + (t49 * t167 + t48 * t188 + (t165 + t212 + t238) * t53) * m(4) + (t58 * t92 - t59 * t91 - t70 * t83 - t71 * t237 - (-t70 * t92 - t71 * t91) * t139) * m(3); t154 + (t20 * t163 + t19 * t166 + (t204 + t267) * t28 + t257) * m(5) + (t48 * t86 - t49 * t85 - t53 * t69 - t54 * t78 - (-t53 * t86 - t54 * t85) * t132) * m(4); ((t132 * t44 - t13) * t131 + (t132 * t43 + t14) * t130) * t253 + ((-t216 * t229 + t228) * t130 + (t161 + (-t256 * t131 + (t230 + t164) * t130) * qJD(4)) * t131) * t199 + ((-t215 * t230 - t228) * t131 + (t161 + (t164 * t130 + (-t256 + t229) * t131) * qJD(4)) * t130) * t196 - t132 * ((t218 * t141 + t217 * t143) * t132 + ((t130 * t244 - t131 * t245) * t143 + (t130 * t246 + t131 * t247) * t141) * qJD(4)) / 0.2e1 + (t132 * t17 + ((-t130 * t258 + t132 * t26) * t131 + (t130 * t259 + t132 * t25) * t130) * t268) * t130 / 0.2e1 - (t132 * t18 + ((t131 * t258 + t132 * t24) * t131 + (-t131 * t259 + t132 * t23) * t130) * t268) * t131 / 0.2e1 + (t9 + t177) * t227 / 0.2e1 + (t10 + t176) * t224 / 0.2e1 + (((-t132 * t28 - t20) * t131 + (t132 * t27 - t19) * t130) * t115 + t266 * t98 + 0.2e1 * ((-t132 * t67 + t42) * t130 + t131 * (t41 + t57)) * t31 - (t27 * t79 - t28 * t80) * t132 - (t31 * (-t130 * t79 - t131 * t80) + t266 * t185) * qJD(4)) * m(5);];
tauc = t1(:);
