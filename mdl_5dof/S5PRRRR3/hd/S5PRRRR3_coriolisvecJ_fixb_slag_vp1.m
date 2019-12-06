% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:04
% EndTime: 2019-12-05 17:06:10
% DurationCPUTime: 2.72s
% Computational Cost: add. (9775->268), mult. (5379->354), div. (0->0), fcn. (4080->8), ass. (0->181)
t143 = pkin(9) + qJ(2);
t137 = sin(t143);
t244 = pkin(2) * qJD(2);
t211 = t137 * t244;
t140 = qJ(3) + t143;
t134 = sin(t140);
t144 = qJD(2) + qJD(3);
t224 = t134 * t144;
t216 = pkin(3) * t224;
t139 = qJD(4) + t144;
t136 = qJ(4) + t140;
t131 = sin(t136);
t132 = cos(t136);
t85 = rSges(5,1) * t131 + rSges(5,2) * t132;
t79 = t139 * t85;
t54 = -t211 - t216 - t79;
t145 = sin(qJ(5));
t229 = t131 * t145;
t109 = rSges(6,2) * t229;
t146 = cos(qJ(5));
t245 = rSges(6,2) * t146;
t208 = qJD(5) * t245;
t227 = t132 * t139;
t182 = rSges(6,3) * t227 + t139 * t109 - t132 * t208;
t217 = qJD(5) * t145;
t206 = t132 * t217;
t222 = t132 * rSges(6,3) + t109;
t228 = t131 * t146;
t67 = rSges(6,1) * t228 - t222;
t58 = t139 * t67;
t127 = t132 * pkin(8);
t88 = pkin(4) * t131 - t127;
t99 = pkin(8) * t227;
t270 = -rSges(6,1) * t206 + t139 * t88 + t182 + t58 + t99;
t135 = cos(t140);
t93 = rSges(4,1) * t134 + rSges(4,2) * t135;
t240 = t144 * t93;
t69 = -t211 - t240;
t269 = 0.2e1 * qJD(5);
t225 = t132 * t146;
t267 = rSges(6,1) * t225 + t131 * rSges(6,3);
t266 = -t132 * pkin(4) - t131 * pkin(8);
t141 = Icges(6,4) * t146;
t184 = -Icges(6,2) * t145 + t141;
t116 = Icges(6,1) * t145 + t141;
t265 = -t216 + t270;
t226 = t132 * t145;
t212 = rSges(6,2) * t226;
t68 = -t212 + t267;
t169 = t68 - t266;
t118 = rSges(6,1) * t145 + t245;
t219 = qJD(5) * t131;
t87 = t118 * t219;
t192 = -t169 * t139 + t87;
t113 = Icges(6,5) * t146 - Icges(6,6) * t145;
t112 = Icges(6,5) * t145 + Icges(6,6) * t146;
t161 = Icges(6,3) * t139 - t112 * qJD(5);
t175 = t184 * t132;
t64 = Icges(6,6) * t131 + t175;
t238 = t145 * t64;
t235 = Icges(6,4) * t145;
t117 = Icges(6,1) * t146 - t235;
t176 = t117 * t132;
t66 = Icges(6,5) * t131 + t176;
t186 = -t146 * t66 + t238;
t230 = t131 * t139;
t264 = -t113 * t230 + t161 * t132 + t186 * t139;
t174 = t113 * t132;
t63 = Icges(6,4) * t228 - Icges(6,2) * t229 - Icges(6,6) * t132;
t239 = t145 * t63;
t108 = Icges(6,4) * t229;
t65 = Icges(6,1) * t228 - Icges(6,5) * t132 - t108;
t187 = -t146 * t65 + t239;
t263 = t161 * t131 + (t174 + t187) * t139;
t114 = Icges(6,2) * t146 + t235;
t183 = t145 * t114 - t146 * t116;
t262 = t113 * qJD(5) + t183 * t139;
t61 = Icges(6,5) * t228 - Icges(6,6) * t229 - Icges(6,3) * t132;
t24 = -t187 * t131 - t132 * t61;
t246 = rSges(6,1) * t146;
t203 = -pkin(4) - t246;
t218 = qJD(5) * t132;
t207 = t118 * t218;
t171 = -t207 - t216;
t159 = t171 - t211;
t28 = (-t67 - t88) * t139 + t159;
t243 = t132 * t28;
t138 = cos(t143);
t210 = t138 * t244;
t223 = t135 * t144;
t215 = pkin(3) * t223;
t172 = t210 + t215;
t29 = t172 - t192;
t148 = (t203 * t243 + (t28 * (-rSges(6,3) - pkin(8)) + t29 * t203) * t131) * t139;
t209 = t139 * t212 + (rSges(6,1) * t217 + t208) * t131;
t261 = t148 + (-t192 + t209) * t28;
t248 = -Icges(6,2) * t228 - t108 + t65;
t250 = t116 * t131 + t63;
t260 = -t248 * t145 - t250 * t146;
t257 = t139 / 0.2e1;
t256 = pkin(2) * t137;
t255 = pkin(2) * qJD(2) ^ 2;
t254 = pkin(3) * t134;
t253 = pkin(3) * t144 ^ 2;
t252 = -t131 * t61 - t65 * t225;
t62 = Icges(6,3) * t131 + t174;
t251 = t131 * t62 + t66 * t225;
t249 = -t116 * t132 - t64;
t247 = -t114 * t132 + t66;
t86 = t132 * rSges(5,1) - rSges(5,2) * t131;
t241 = t139 * t86;
t232 = t112 * t132;
t47 = -t183 * t131 - t232;
t237 = t47 * t139;
t233 = t112 * t131;
t231 = t113 * t139;
t221 = -t114 + t117;
t220 = t116 + t184;
t214 = t137 * t255;
t213 = t138 * t255;
t202 = -t219 / 0.2e1;
t199 = t218 / 0.2e1;
t51 = t66 * t228;
t197 = t132 * t62 - t51;
t196 = -t61 + t238;
t94 = t135 * rSges(4,1) - rSges(4,2) * t134;
t72 = rSges(5,1) * t227 - rSges(5,2) * t230;
t84 = rSges(4,1) * t223 - rSges(4,2) * t224;
t130 = pkin(3) * t135;
t193 = t130 + t86;
t190 = -rSges(6,2) * t145 + t246;
t189 = -t131 * t29 - t243;
t188 = t131 * t67 + t132 * t68;
t41 = t145 * t65 + t146 * t63;
t42 = t145 * t66 + t146 * t64;
t25 = -t64 * t229 - t197;
t180 = (t131 * t25 - t132 * t24) * qJD(5);
t26 = -t63 * t226 - t252;
t27 = -t64 * t226 + t251;
t179 = (t131 * t27 - t132 * t26) * qJD(5);
t178 = -t134 * t253 - t214;
t177 = -t135 * t253 - t213;
t170 = -t85 - t254;
t168 = -t72 - t215;
t167 = -t247 * t145 + t249 * t146;
t166 = t203 * t131 + t127 + t222;
t165 = t130 + t169;
t164 = (-t220 * t145 + t221 * t146) * t139;
t163 = Icges(6,5) * t139 - qJD(5) * t116;
t162 = Icges(6,6) * t139 - t114 * qJD(5);
t158 = t166 - t254;
t43 = (-t139 * t228 - t206) * rSges(6,1) + t182;
t44 = t267 * t139 - t209;
t157 = (t43 + t58) * t132 + (-t139 * t68 + t44) * t131;
t48 = -t183 * t132 + t233;
t46 = t48 * t139;
t10 = t46 + t179;
t103 = t184 * qJD(5);
t104 = t117 * qJD(5);
t14 = -t187 * qJD(5) + t145 * (t163 * t131 + t139 * t176) + t146 * (t162 * t131 + t139 * t175);
t15 = -t186 * qJD(5) + t145 * (-t117 * t230 + t163 * t132) + t146 * (t162 * t132 - t184 * t230);
t149 = -t103 * t145 + t104 * t146 + t112 * t139 + (-t114 * t146 - t116 * t145) * qJD(5);
t18 = t131 * t262 + t149 * t132;
t19 = t149 * t131 - t132 * t262;
t9 = t180 + t237;
t156 = (t46 + ((t25 - t51 + (t62 + t239) * t132 + t252) * t132 + t251 * t131) * qJD(5)) * t199 + (-t183 * qJD(5) + t103 * t146 + t104 * t145) * t139 + (t9 - t237 + ((t196 * t132 - t251 + t27) * t132 + (t196 * t131 + t197 + t26) * t131) * qJD(5)) * t202 + (t15 + t18) * t219 / 0.2e1 - (t10 + t14 + t19) * t218 / 0.2e1 + ((t41 + t47) * t131 + (t42 + t48) * t132) * qJD(5) * t257;
t133 = pkin(2) * t138;
t105 = t190 * qJD(5);
t81 = t118 * t132;
t80 = t118 * t131;
t70 = t144 * t94 + t210;
t60 = -t144 * t84 - t213;
t59 = -t144 * t240 - t214;
t55 = t172 + t241;
t50 = -t139 * t72 + t177;
t49 = -t139 * t79 + t178;
t32 = t188 * qJD(5) + qJD(1);
t21 = -t105 * t218 + (t266 * t139 - t44 + t87) * t139 + t177;
t20 = -t105 * t219 + (-pkin(4) * t230 - t207 + t43 + t99) * t139 + t178;
t11 = t157 * qJD(5);
t1 = [m(6) * t11; t156 + m(4) * (t60 * (-t93 - t256) + t59 * (t133 + t94) + (-t84 - t210 + t70) * t69) + (t21 * (t158 - t256) + t28 * (-t172 + t209) + t20 * (t133 + t165) + t148 + (-t159 + t28 - t211 + t265) * t29) * m(6) + m(5) * (t50 * (t170 - t256) + t49 * (t133 + t193) + (t55 + t168 - t210) * t54); t156 + (t21 * t158 + t20 * t165 + (-t171 + t265) * t29 + t261) * m(6) + (t50 * t170 + t49 * t193 + (t168 + t215 + t241) * t54) * m(5) + (-(-t69 * t94 - t70 * t93) * t144 + t59 * t94 - t60 * t93 - t69 * t84 - t70 * t240) * m(4); t156 + (t21 * t166 + t20 * t169 + (t207 + t270) * t29 + t261) * m(6) + (-(-t54 * t86 - t55 * t85) * t139 + t49 * t86 - t50 * t85 - t54 * t72 - t55 * t79) * m(5); ((t139 * t42 - t14) * t132 + (t139 * t41 + t15) * t131) * t257 + ((-t219 * t232 + t231) * t131 + (t164 + (-t260 * t132 + (t233 + t167) * t131) * qJD(5)) * t132) * t202 + ((-t218 * t233 - t231) * t132 + (t164 + (t167 * t131 + (-t260 + t232) * t132) * qJD(5)) * t131) * t199 - t139 * ((t221 * t145 + t220 * t146) * t139 + ((t247 * t131 - t248 * t132) * t146 + (t249 * t131 + t250 * t132) * t145) * qJD(5)) / 0.2e1 + (t139 * t18 + ((-t263 * t131 + t139 * t27) * t132 + (t264 * t131 + t139 * t26) * t131) * t269) * t131 / 0.2e1 - (t139 * t19 + ((t263 * t132 + t139 * t25) * t132 + (-t264 * t132 + t139 * t24) * t131) * t269) * t132 / 0.2e1 + (t9 + t180) * t230 / 0.2e1 + (t10 + t179) * t227 / 0.2e1 + (t11 * t188 + t32 * t157 + t189 * t105 + ((-t139 * t29 - t21) * t132 + (t139 * t28 - t20) * t131) * t118 - (t28 * t80 - t29 * t81) * t139 - (t32 * (-t131 * t80 - t132 * t81) + t189 * t190) * qJD(5)) * m(6);];
tauc = t1(:);
