% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:28
% EndTime: 2019-12-31 17:01:32
% DurationCPUTime: 2.35s
% Computational Cost: add. (5619->265), mult. (4334->357), div. (0->0), fcn. (3282->8), ass. (0->173)
t138 = qJ(1) + qJ(2);
t131 = pkin(7) + t138;
t129 = cos(t131);
t133 = cos(t138);
t130 = pkin(2) * t133;
t255 = -t129 * rSges(4,1) - t130;
t140 = sin(qJ(1));
t230 = pkin(1) * qJD(1);
t199 = t140 * t230;
t137 = qJD(1) + qJD(2);
t132 = sin(t138);
t90 = rSges(3,1) * t132 + rSges(3,2) * t133;
t228 = t137 * t90;
t69 = -t199 - t228;
t254 = 2 * qJD(4);
t128 = sin(t131);
t253 = -rSges(4,2) * t128 - t255;
t141 = cos(qJ(4));
t213 = t129 * t141;
t106 = rSges(5,1) * t213;
t252 = t128 * rSges(5,3) + t106;
t126 = t129 * pkin(3);
t87 = t128 * pkin(6) + t126;
t134 = Icges(5,4) * t141;
t139 = sin(qJ(4));
t172 = -Icges(5,2) * t139 + t134;
t111 = Icges(5,1) * t139 + t134;
t214 = t129 * t139;
t200 = rSges(5,2) * t214;
t68 = -t200 + t252;
t251 = t130 + t68 + t87;
t215 = t129 * t137;
t101 = pkin(6) * t215;
t217 = t128 * t139;
t105 = rSges(5,2) * t217;
t231 = rSges(5,2) * t141;
t197 = qJD(4) * t231;
t170 = rSges(5,3) * t215 + t137 * t105 - t129 * t197;
t206 = qJD(4) * t139;
t195 = t129 * t206;
t211 = t129 * rSges(5,3) + t105;
t216 = t128 * t141;
t67 = rSges(5,1) * t216 - t211;
t58 = t137 * t67;
t125 = t129 * pkin(6);
t86 = pkin(3) * t128 - t125;
t250 = -rSges(5,1) * t195 + t137 * t86 + t101 + t170 + t58;
t108 = Icges(5,5) * t141 - Icges(5,6) * t139;
t223 = Icges(5,4) * t139;
t109 = Icges(5,2) * t141 + t223;
t171 = t139 * t109 - t141 * t111;
t249 = t108 * qJD(4) + t137 * t171;
t107 = Icges(5,5) * t139 + Icges(5,6) * t141;
t154 = Icges(5,3) * t137 - qJD(4) * t107;
t163 = t172 * t129;
t64 = Icges(5,6) * t128 + t163;
t226 = t139 * t64;
t112 = Icges(5,1) * t141 - t223;
t164 = t112 * t129;
t66 = Icges(5,5) * t128 + t164;
t174 = -t141 * t66 + t226;
t218 = t128 * t137;
t248 = -t108 * t218 + t129 * t154 + t137 * t174;
t162 = t108 * t129;
t63 = Icges(5,4) * t216 - Icges(5,2) * t217 - Icges(5,6) * t129;
t227 = t139 * t63;
t104 = Icges(5,4) * t217;
t65 = Icges(5,1) * t216 - Icges(5,5) * t129 - t104;
t175 = -t141 * t65 + t227;
t247 = t128 * t154 + (t162 + t175) * t137;
t114 = rSges(5,1) * t139 + t231;
t208 = qJD(4) * t128;
t83 = t114 * t208;
t246 = t137 * t251 - t83;
t61 = Icges(5,5) * t216 - Icges(5,6) * t217 - Icges(5,3) * t129;
t24 = -t128 * t175 - t129 * t61;
t234 = -Icges(5,2) * t216 - t104 + t65;
t236 = t111 * t128 + t63;
t245 = -t139 * t234 - t141 * t236;
t136 = t137 ^ 2;
t242 = t137 / 0.2e1;
t241 = pkin(1) * t140;
t240 = pkin(2) * t132;
t239 = pkin(2) * t136;
t142 = cos(qJ(1));
t135 = t142 * pkin(1);
t238 = -t128 * t61 - t65 * t213;
t62 = Icges(5,3) * t128 + t162;
t237 = t128 * t62 + t66 * t213;
t235 = -t111 * t129 - t64;
t233 = -t109 * t129 + t66;
t232 = rSges(5,1) * t141;
t127 = t133 * rSges(3,1);
t220 = t107 * t129;
t47 = -t128 * t171 - t220;
t225 = t47 * t137;
t221 = t107 * t128;
t219 = t108 * t137;
t212 = t132 * t137;
t210 = -t109 + t112;
t209 = t111 + t172;
t207 = qJD(4) * t129;
t198 = t142 * t230;
t29 = t198 + t246;
t205 = t29 * t240;
t204 = t137 * t200 + (rSges(5,1) * t206 + t197) * t128;
t203 = pkin(2) * t212;
t143 = qJD(1) ^ 2;
t202 = t143 * t241;
t201 = t143 * t135;
t196 = t114 * t207;
t192 = -pkin(3) - t232;
t84 = rSges(4,1) * t128 + rSges(4,2) * t129;
t160 = -t84 - t240;
t190 = -t208 / 0.2e1;
t187 = t207 / 0.2e1;
t51 = t66 * t216;
t185 = t129 * t62 - t51;
t44 = (-t137 * t216 - t195) * rSges(5,1) + t170;
t184 = t44 + t58;
t45 = t137 * t252 - t204;
t183 = -t137 * t68 + t45;
t182 = -t61 + t226;
t91 = -rSges(3,2) * t132 + t127;
t81 = -rSges(3,2) * t212 + t127 * t137;
t176 = -rSges(5,2) * t139 + t232;
t35 = t139 * t65 + t141 * t63;
t36 = t139 * t66 + t141 * t64;
t25 = -t217 * t64 - t185;
t168 = (t128 * t25 - t129 * t24) * qJD(4);
t26 = -t214 * t63 - t238;
t27 = -t214 * t64 + t237;
t167 = (t128 * t27 - t129 * t26) * qJD(4);
t166 = -t132 * t239 - t202;
t165 = -t133 * t239 - t201;
t161 = -t196 - t199;
t159 = -t139 * t233 + t141 * t235;
t157 = (-t139 * t209 + t141 * t210) * t137;
t156 = Icges(5,5) * t137 - qJD(4) * t111;
t155 = Icges(5,6) * t137 - qJD(4) * t109;
t153 = t128 * t192 + t125 + t211 - t240;
t48 = -t129 * t171 + t221;
t46 = t48 * t137;
t10 = t46 + t167;
t14 = -qJD(4) * t175 + t139 * (t128 * t156 + t137 * t164) + t141 * (t128 * t155 + t137 * t163);
t15 = -qJD(4) * t174 + t139 * (-t112 * t218 + t129 * t156) + t141 * (t129 * t155 - t172 * t218);
t93 = t172 * qJD(4);
t94 = t112 * qJD(4);
t146 = t107 * t137 - t139 * t93 + t141 * t94 + (-t109 * t141 - t111 * t139) * qJD(4);
t18 = t128 * t249 + t146 * t129;
t19 = t146 * t128 - t129 * t249;
t9 = t168 + t225;
t152 = (-qJD(4) * t171 + t139 * t94 + t141 * t93) * t137 + (t46 + ((t25 - t51 + (t62 + t227) * t129 + t238) * t129 + t237 * t128) * qJD(4)) * t187 + (t9 - t225 + ((t129 * t182 - t237 + t27) * t129 + (t128 * t182 + t185 + t26) * t128) * qJD(4)) * t190 + (t15 + t18) * t208 / 0.2e1 - (t14 + t19 + t10) * t207 / 0.2e1 + ((t35 + t47) * t128 + (t36 + t48) * t129) * qJD(4) * t242;
t56 = t137 * t160 - t199;
t57 = t137 * t253 + t198;
t145 = (t57 * t160 + t255 * t56) * t137;
t28 = (-t67 - t86 - t240) * t137 + t161;
t144 = (t28 * (-t106 - t126 - t130) - t205 + (t28 * (-rSges(5,3) - pkin(6)) + t29 * t192) * t128) * t137;
t99 = rSges(4,2) * t218;
t95 = t176 * qJD(4);
t79 = t137 * t84;
t78 = t114 * t129;
t77 = t114 * t128;
t70 = t137 * t91 + t198;
t60 = -t137 * t81 - t201;
t59 = -t137 * t228 - t202;
t50 = -t137 * (rSges(4,1) * t215 - t99) + t165;
t49 = -t136 * t84 + t166;
t32 = qJD(3) + (t128 * t67 + t129 * t68) * qJD(4);
t23 = -t95 * t207 + (-t137 * t87 - t45 + t83) * t137 + t165;
t22 = -t95 * t208 + (-pkin(3) * t218 + t101 - t196 + t44) * t137 + t166;
t11 = (t183 * t128 + t184 * t129) * qJD(4);
t1 = [m(3) * (t60 * (-t90 - t241) + t59 * (t135 + t91) + (-t81 - t198 + t70) * t69) + t152 + (t23 * (t153 - t241) + t28 * (-t198 + t204) + t22 * (t135 + t251) + t144 + (-t161 + t28 + t203 - t199 + t250) * t29) * m(5) + (-(-t56 - t79 - t203) * t57 + t50 * (t160 - t241) + t56 * (t99 - t198) + t49 * (t135 + t253) + t145) * m(4); t152 + (t205 * t137 + t23 * t153 + t22 * t251 + t144 + (t196 + t250) * t29 + (t204 + t246) * t28) * m(5) + (t57 * t79 - (-t240 * t57 - t253 * t56) * t137 + t50 * t160 + t49 * t253 + t56 * t99 + t145) * m(4) + (-(-t69 * t91 - t70 * t90) * t137 + t59 * t91 - t60 * t90 - t69 * t81 - t70 * t228) * m(3); m(5) * t11; ((t137 * t36 - t14) * t129 + (t137 * t35 + t15) * t128) * t242 + ((-t208 * t220 + t219) * t128 + (t157 + (-t245 * t129 + (t221 + t159) * t128) * qJD(4)) * t129) * t190 + ((-t207 * t221 - t219) * t129 + (t157 + (t159 * t128 + (-t245 + t220) * t129) * qJD(4)) * t128) * t187 - t137 * ((t210 * t139 + t209 * t141) * t137 + ((t233 * t128 - t234 * t129) * t141 + (t128 * t235 + t236 * t129) * t139) * qJD(4)) / 0.2e1 + (t137 * t18 + ((-t128 * t247 + t137 * t27) * t129 + (t128 * t248 + t137 * t26) * t128) * t254) * t128 / 0.2e1 - (t137 * t19 + ((t129 * t247 + t137 * t25) * t129 + (-t129 * t248 + t137 * t24) * t128) * t254) * t129 / 0.2e1 + (t9 + t168) * t218 / 0.2e1 + (t10 + t167) * t215 / 0.2e1 + ((t11 * t68 + t32 * t184 - t28 * t95) * t129 + (t11 * t67 + t32 * t183 - t29 * t95) * t128 + ((-t137 * t29 - t23) * t129 + (t137 * t28 - t22) * t128) * t114 - (t28 * t77 - t29 * t78) * t137 - (t32 * (-t128 * t77 - t129 * t78) + (-t128 * t29 - t129 * t28) * t176) * qJD(4)) * m(5);];
tauc = t1(:);
