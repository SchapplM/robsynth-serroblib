% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:17
% EndTime: 2019-12-31 17:32:21
% DurationCPUTime: 3.40s
% Computational Cost: add. (4260->279), mult. (7335->378), div. (0->0), fcn. (7728->8), ass. (0->155)
t192 = sin(pkin(7));
t193 = cos(pkin(7));
t224 = sin(qJ(3));
t225 = cos(qJ(3));
t113 = -t192 * t224 - t193 * t225;
t114 = -t192 * t225 + t193 * t224;
t123 = pkin(8) + qJ(5);
t121 = sin(t123);
t122 = cos(t123);
t190 = Icges(6,4) * t122;
t148 = -Icges(6,2) * t121 + t190;
t48 = -Icges(6,6) * t114 + t113 * t148;
t191 = Icges(6,4) * t121;
t150 = Icges(6,1) * t122 - t191;
t51 = -Icges(6,5) * t114 + t113 * t150;
t153 = -t121 * t48 + t122 * t51;
t146 = Icges(6,5) * t122 - Icges(6,6) * t121;
t45 = -Icges(6,3) * t114 + t113 * t146;
t17 = t113 * t153 - t114 * t45;
t46 = -Icges(6,3) * t113 - t114 * t146;
t241 = t114 * t46;
t147 = Icges(6,2) * t122 + t191;
t149 = Icges(6,1) * t121 + t190;
t144 = -t121 * t147 + t122 * t149;
t145 = Icges(6,5) * t121 + Icges(6,6) * t122;
t60 = t145 * t114;
t36 = t113 * t144 - t60;
t240 = qJD(3) * t36;
t24 = -t121 * t51 - t122 * t48;
t61 = t145 * t113;
t159 = rSges(6,1) * t121 + rSges(6,2) * t122;
t179 = qJD(5) * t113;
t238 = t159 * t179;
t125 = cos(pkin(8));
t237 = rSges(5,2) * sin(pkin(8)) - rSges(5,1) * t125;
t96 = t114 * qJD(3);
t236 = 0.2e1 * qJD(5);
t101 = t114 * qJ(4);
t223 = t113 * pkin(3);
t76 = t101 - t223;
t180 = qJD(4) * t114;
t58 = -t113 * rSges(5,3) + t237 * t114;
t104 = t114 * pkin(3);
t166 = -qJ(4) * t113 - t104;
t73 = qJD(3) * t166;
t234 = qJD(3) * t58 + t180 + t73;
t100 = t113 * qJD(4);
t97 = qJD(3) * t113;
t169 = t97 * pkin(3) - qJ(4) * t96;
t126 = -pkin(6) - qJ(4);
t196 = t97 * t126;
t119 = pkin(4) * t125 + pkin(3);
t219 = -pkin(3) + t119;
t177 = qJD(5) * t121;
t140 = t113 * t177 - t122 * t96;
t176 = qJD(5) * t122;
t141 = t113 * t176 + t121 * t96;
t220 = t97 * rSges(6,3);
t33 = rSges(6,1) * t140 + rSges(6,2) * t141 - t220;
t87 = t97 * qJ(4);
t208 = t87 - t180;
t59 = -t96 * pkin(3) - t208;
t84 = qJD(4) * t96;
t203 = rSges(6,2) * t121;
t160 = rSges(6,1) * t122 - t203;
t94 = t160 * qJD(5);
t10 = -t84 + (t114 * t94 - t159 * t97) * qJD(5) + (t219 * t96 - t196 - t33 - t59 - t87) * qJD(3);
t167 = qJD(2) * t193;
t151 = -t100 - t167;
t178 = qJD(5) * t114;
t54 = -t114 * rSges(6,3) + t113 * t160;
t19 = t159 * t178 + (t113 * t219 + t114 * t126 + t101 + t54 - t76) * qJD(3) + t151;
t233 = t10 * t113 + t19 * t96;
t207 = t113 * t126 - t114 * t119;
t185 = t114 * t121;
t184 = t114 * t122;
t195 = -rSges(6,1) * t184 - t113 * rSges(6,3);
t55 = rSges(6,2) * t185 + t195;
t232 = (-t166 + t207 + t55) * qJD(3) + t238 + t73;
t231 = qJD(3) * t76;
t52 = -Icges(6,5) * t113 - t114 * t150;
t213 = t147 * t114 + t52;
t49 = -Icges(6,6) * t113 - t114 * t148;
t215 = -t149 * t114 + t49;
t230 = t121 * t213 + t122 * t215;
t221 = t97 * rSges(5,3);
t218 = t113 * t46 + t52 * t184;
t198 = t122 * t52;
t217 = -t113 * t198 + t241;
t214 = -t149 * t113 - t48;
t212 = t147 * t113 - t51;
t211 = qJD(3) * (t169 - t100) - qJD(4) * t97;
t197 = t122 * t97;
t210 = rSges(6,1) * t197 - t96 * rSges(6,3);
t209 = t97 * t119 + t96 * t126;
t206 = m(4) * qJD(3);
t202 = t114 * rSges(5,3);
t200 = t121 * t49;
t199 = t121 * t97;
t35 = t114 * t144 + t61;
t183 = t35 * qJD(3);
t182 = -t147 + t150;
t181 = -t148 - t149;
t175 = t146 * qJD(3);
t174 = t114 * t177;
t171 = t179 / 0.2e1;
t170 = -t178 / 0.2e1;
t162 = rSges(4,1) * t113 + rSges(4,2) * t114;
t120 = qJD(2) * t192;
t18 = t120 + t180 + t232;
t158 = -t113 * t18 - t114 * t19;
t157 = -t113 * t54 + t114 * t55;
t155 = -t198 + t200;
t152 = -t96 * rSges(5,3) - t237 * t97;
t15 = t45 * t113 + t184 * t51 - t185 * t48;
t143 = pkin(3) - t237;
t139 = t114 * t176 - t199;
t138 = t174 + t197;
t16 = t113 * t200 + t217;
t137 = qJD(5) * (-t16 * t113 + t17 * t114);
t14 = t185 * t49 - t218;
t136 = (-t14 * t113 + t15 * t114) * qJD(5);
t135 = t121 * t212 + t122 * t214;
t134 = -t113 * t237 - t202;
t34 = rSges(6,1) * t174 + rSges(6,2) * t139 + t210;
t133 = t113 * t33 + t114 * t34 - t54 * t96 - t55 * t97;
t132 = (t121 * t181 + t122 * t182) * qJD(3);
t92 = t148 * qJD(5);
t93 = t150 * qJD(5);
t127 = -t121 * t92 + t122 * t93 + (-t121 * t149 - t122 * t147) * qJD(5);
t91 = t146 * qJD(5);
t75 = -rSges(4,1) * t114 + rSges(4,2) * t113;
t71 = rSges(4,1) * t97 + rSges(4,2) * t96;
t70 = -rSges(4,1) * t96 + rSges(4,2) * t97;
t69 = qJD(3) * t162 - t167;
t68 = qJD(3) * t75 + t120;
t67 = t159 * t113;
t66 = t159 * t114;
t28 = Icges(6,5) * t138 + Icges(6,6) * t139 - Icges(6,3) * t96;
t27 = Icges(6,5) * t140 + Icges(6,6) * t141 - Icges(6,3) * t97;
t26 = (-t76 + t134) * qJD(3) + t151;
t25 = t120 + t234;
t23 = t121 * t52 + t122 * t49;
t22 = qJD(3) * t152 + t211;
t21 = -t84 + (-t237 * t96 + t221 - t59) * qJD(3);
t20 = qJD(5) * t157 + qJD(1);
t13 = t113 * t91 + t114 * t127 - t144 * t97 + t145 * t96;
t12 = t113 * t127 - t114 * t91 + t144 * t96 + t145 * t97;
t11 = (t113 * t94 + t159 * t96) * qJD(5) + (-t169 + t34 + t209) * qJD(3) + t211;
t9 = qJD(5) * t153 - t121 * (Icges(6,1) * t140 + Icges(6,4) * t141 - Icges(6,5) * t97) - t122 * (Icges(6,4) * t140 + Icges(6,2) * t141 - Icges(6,6) * t97);
t8 = qJD(5) * t155 - t121 * (Icges(6,1) * t138 + Icges(6,4) * t139 - Icges(6,5) * t96) - t122 * (Icges(6,4) * t138 + Icges(6,2) * t139 - Icges(6,6) * t96);
t7 = t133 * qJD(5);
t6 = t137 - t240;
t5 = t136 - t183;
t1 = [m(6) * t7; (t192 * t71 + t193 * t70) * t206 + m(5) * (t192 * t22 - t193 * t21) + m(6) * (-t10 * t193 + t11 * t192); -qJD(3) * (-qJD(5) * t144 - t121 * t93 - t122 * t92) + t5 * t178 / 0.2e1 + m(4) * (t68 * t71 - t69 * t70 + (-t162 * t70 + t71 * t75) * qJD(3)) - (t162 * t68 - t69 * t75) * t206 + (t11 * (t195 + t207) + t233 * (t119 + t160) + (t11 * t203 + t10 * (-rSges(6,3) + t126)) * t114 + (-rSges(6,2) * t199 + t209 + t210 - (t113 * t119 - t223 + t54) * qJD(3) - (qJ(4) + t126) * t96 + t231) * t18 + (-t196 + t220 + t232 - t238) * t19) * m(6) + (t22 * (-t104 + t58) + t21 * (-t101 - t202) + (-t22 * qJ(4) + t143 * t21) * t113 + (t143 * t96 + t208 + t221 + t234) * t26 + (-qJD(3) * t134 + t152 + t169 + t231) * t25) * m(5) + (-t183 + ((t15 - t16 + t217) * t114 + t218 * t113) * qJD(5) + t9 + t12) * t170 + (t8 + t13 + t6 + ((-t14 + (t45 + t200) * t114 - t218) * t114 + (-t15 + (t155 + t45) * t113 + t241) * t113) * qJD(5) + t240) * t171 + ((-t23 + t35) * t96 + (-t24 + t36) * t97) * qJD(5) / 0.2e1; (-t113 * t21 + t114 * t22 - t25 * t97 - t26 * t96 - (-t113 * t25 - t114 * t26) * qJD(3)) * m(5) + (-qJD(3) * t158 + t11 * t114 - t18 * t97 - t233) * m(6); -qJD(3) * (-t8 * t113 + t9 * t114 + t23 * t96 + t24 * t97) / 0.2e1 + ((t61 * t178 + t175) * t114 + (-t132 + (-t230 * t113 + (-t60 + t135) * t114) * qJD(5)) * t113) * t170 + ((t60 * t179 - t175) * t113 + (-t132 + (t135 * t114 + (-t230 - t61) * t113) * qJD(5)) * t114) * t171 + qJD(3) * (-(t182 * t121 - t181 * t122) * qJD(3) + ((t113 * t213 - t114 * t212) * t122 + (-t113 * t215 + t114 * t214) * t121) * qJD(5)) / 0.2e1 - (-t13 * qJD(3) + (-t113 * (-t113 * t28 - t155 * t97 - t46 * t96) + (-t113 * t27 - t153 * t97 + t45 * t96) * t114 - t14 * t96 - t15 * t97) * t236) * t113 / 0.2e1 + (-qJD(3) * t12 + (-(t114 * t28 + t155 * t96 - t46 * t97) * t113 + t114 * (t114 * t27 + t153 * t96 + t45 * t97) - t16 * t96 - t17 * t97) * t236) * t114 / 0.2e1 - (t5 + t136) * t96 / 0.2e1 - (t6 + t137) * t97 / 0.2e1 + (t7 * t157 + t20 * t133 - t158 * t94 - (-t10 * t114 - t11 * t113 - t18 * t96 + t19 * t97) * t159 - (t18 * t66 - t19 * t67) * qJD(3) - (t20 * (t113 * t67 + t114 * t66) - t158 * t160) * qJD(5)) * m(6);];
tauc = t1(:);
