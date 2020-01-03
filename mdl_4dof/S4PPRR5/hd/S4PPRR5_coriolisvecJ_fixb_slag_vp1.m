% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:41
% EndTime: 2019-12-31 16:19:48
% DurationCPUTime: 5.29s
% Computational Cost: add. (3393->383), mult. (9558->615), div. (0->0), fcn. (9282->6), ass. (0->194)
t133 = sin(pkin(6));
t131 = t133 ^ 2;
t134 = cos(pkin(6));
t132 = t134 ^ 2;
t245 = t131 + t132;
t138 = cos(qJ(3));
t194 = qJD(4) * t138;
t199 = qJD(3) * t133;
t120 = t134 * t194 + t199;
t235 = t120 / 0.2e1;
t198 = qJD(3) * t134;
t121 = -t133 * t194 + t198;
t233 = t121 / 0.2e1;
t136 = sin(qJ(3));
t135 = sin(qJ(4));
t137 = cos(qJ(4));
t179 = rSges(5,1) * t137 - rSges(5,2) * t135;
t91 = t136 * rSges(5,3) + t138 * t179;
t246 = t120 * t91;
t226 = t134 * t91;
t159 = Icges(4,5) * t138 - Icges(4,6) * t136;
t108 = t159 * t133;
t92 = qJD(3) * t108;
t109 = t159 * t134;
t244 = qJD(3) * t109;
t125 = rSges(4,1) * t138 - t136 * rSges(4,2);
t243 = qJD(3) * t125;
t197 = qJD(3) * t136;
t196 = qJD(3) * t138;
t158 = Icges(5,5) * t137 - Icges(5,6) * t135;
t85 = Icges(5,3) * t136 + t138 * t158;
t208 = Icges(5,4) * t137;
t160 = -Icges(5,2) * t135 + t208;
t87 = Icges(5,6) * t136 + t138 * t160;
t209 = Icges(5,4) * t135;
t163 = Icges(5,1) * t137 - t209;
t89 = Icges(5,5) * t136 + t138 * t163;
t212 = Icges(4,4) * t138;
t162 = -Icges(4,2) * t136 + t212;
t111 = t162 * t134;
t213 = Icges(4,4) * t136;
t164 = Icges(4,1) * t138 - t213;
t113 = t164 * t134;
t161 = Icges(4,2) * t138 + t213;
t77 = Icges(4,6) * t133 - t134 * t161;
t79 = Icges(4,5) * t133 + (-Icges(4,1) * t136 - t212) * t134;
t242 = (t136 * t77 - t138 * t79) * qJD(3) + t111 * t196 + t113 * t197;
t112 = t164 * t133;
t76 = Icges(4,6) * t134 + t133 * t161;
t203 = t133 * t138;
t128 = Icges(4,4) * t203;
t204 = t133 * t136;
t78 = Icges(4,1) * t204 + Icges(4,5) * t134 + t128;
t241 = -t138 * t162 * t199 + qJD(3) * (t136 * t76 - t138 * t78) - t112 * t197;
t201 = t135 * t136;
t106 = t133 * t137 + t134 * t201;
t202 = t134 * t138;
t200 = t136 * t137;
t107 = -t133 * t135 + t134 * t200;
t210 = Icges(5,4) * t107;
t47 = Icges(5,2) * t106 + Icges(5,6) * t202 - t210;
t103 = Icges(5,4) * t106;
t49 = -Icges(5,1) * t107 + Icges(5,5) * t202 + t103;
t170 = t135 * t47 - t137 * t49;
t104 = -t133 * t201 + t134 * t137;
t105 = t133 * t200 + t134 * t135;
t211 = Icges(5,4) * t105;
t46 = Icges(5,2) * t104 - Icges(5,6) * t203 + t211;
t102 = Icges(5,4) * t104;
t48 = Icges(5,1) * t105 - Icges(5,5) * t203 + t102;
t171 = t135 * t46 - t137 * t48;
t139 = t120 * (-t134 * t85 + t170) + t121 * (t85 * t133 + t171);
t240 = t136 * (-t113 - t77) + t138 * (-t111 + t79);
t239 = -t136 * (t112 - t76) - t138 * (-Icges(4,2) * t204 + t128 + t78);
t118 = (-Icges(5,1) * t135 - t208) * t138;
t195 = qJD(4) * t136;
t238 = t120 * (-Icges(5,1) * t106 - t210 + t47) + t121 * (-Icges(5,1) * t104 + t211 + t46) - t195 * (t118 - t87);
t117 = (-Icges(5,2) * t137 - t209) * t138;
t143 = t120 * (Icges(5,2) * t107 + t103 + t49) + t121 * (-Icges(5,2) * t105 + t102 + t48) + t195 * (t117 + t89);
t192 = t134 * t197;
t191 = t135 * t196;
t66 = qJD(4) * t107 + t134 * t191;
t190 = t137 * t196;
t67 = qJD(4) * t106 - t134 * t190;
t35 = Icges(5,5) * t67 + Icges(5,6) * t66 - Icges(5,3) * t192;
t45 = -Icges(5,5) * t107 + Icges(5,6) * t106 + Icges(5,3) * t202;
t151 = -t138 * t35 + t197 * t45;
t37 = Icges(5,4) * t67 + Icges(5,2) * t66 - Icges(5,6) * t192;
t39 = Icges(5,1) * t67 + Icges(5,4) * t66 - Icges(5,5) * t192;
t10 = t106 * t37 - t107 * t39 - t134 * t151 + t66 * t47 + t67 * t49;
t116 = (-Icges(5,5) * t135 - Icges(5,6) * t137) * t138;
t84 = Icges(5,3) * t138 - t136 * t158;
t52 = qJD(3) * t84 + qJD(4) * t116;
t150 = -t138 * t52 + t197 * t85;
t21 = t106 * t47 - t107 * t49 + t202 * t45;
t44 = Icges(5,5) * t105 + Icges(5,6) * t104 - Icges(5,3) * t203;
t20 = t106 * t46 - t107 * t48 + t202 * t44;
t228 = t133 * t20;
t176 = -t134 * t21 + t228;
t32 = t106 * t87 - t107 * t89 + t202 * t85;
t222 = t32 * t138;
t86 = Icges(5,6) * t138 - t136 * t160;
t53 = qJD(3) * t86 + qJD(4) * t117;
t88 = Icges(5,5) * t138 - t136 * t163;
t54 = qJD(3) * t88 + qJD(4) * t118;
t141 = (t106 * t53 - t107 * t54 - t134 * t150 + t66 * t87 + t67 * t89) * t136 + (t136 * t176 + t222) * qJD(3);
t193 = t133 * t197;
t64 = -qJD(4) * t105 - t133 * t191;
t65 = qJD(4) * t104 + t133 * t190;
t34 = Icges(5,5) * t65 + Icges(5,6) * t64 + Icges(5,3) * t193;
t152 = -t138 * t34 + t197 * t44;
t36 = Icges(5,4) * t65 + Icges(5,2) * t64 + Icges(5,6) * t193;
t38 = Icges(5,1) * t65 + Icges(5,4) * t64 + Icges(5,5) * t193;
t9 = t106 * t36 - t107 * t38 - t134 * t152 + t66 * t46 + t67 * t48;
t237 = qJD(4) * t141 / 0.2e1 + t10 * t235 + t9 * t233;
t236 = -t120 / 0.2e1;
t234 = -t121 / 0.2e1;
t19 = t104 * t47 + t105 * t49 - t203 * t45;
t227 = t134 * t19;
t224 = t136 * t85;
t31 = t104 * t87 + t105 * t89 - t203 * t85;
t223 = t31 * t138;
t181 = t136 * pkin(3) - pkin(5) * t138;
t123 = t181 * qJD(3);
t119 = (-rSges(5,1) * t135 - rSges(5,2) * t137) * t138;
t90 = rSges(5,3) * t138 - t136 * t179;
t55 = qJD(3) * t90 + qJD(4) * t119;
t215 = -t123 + t55;
t127 = pkin(3) * t138 + t136 * pkin(5);
t214 = t127 + t91;
t148 = qJD(3) * t127;
t98 = t125 * t199;
t189 = -t198 / 0.2e1;
t188 = t196 / 0.2e1;
t187 = -t195 / 0.2e1;
t186 = t195 / 0.2e1;
t185 = qJD(3) * t245;
t184 = t193 / 0.2e1;
t183 = t136 * t189;
t182 = qJD(4) * t188;
t180 = t136 * rSges(4,1) + rSges(4,2) * t138;
t178 = -t45 * t120 - t44 * t121;
t18 = t104 * t46 + t105 * t48 - t203 * t44;
t177 = t133 * t18 - t227;
t23 = t136 * t44 - t138 * t171;
t24 = t136 * t45 - t138 * t170;
t175 = t23 * t133 - t24 * t134;
t50 = t105 * rSges(5,1) + t104 * rSges(5,2) - rSges(5,3) * t203;
t51 = -t107 * rSges(5,1) + t106 * rSges(5,2) + rSges(5,3) * t202;
t174 = t133 * t51 + t134 * t50;
t173 = t245 * t180;
t172 = -t125 * t134 * t198 - t133 * t98;
t169 = t135 * t87 - t137 * t89;
t157 = qJD(4) * t184;
t156 = qJD(4) * t183;
t153 = t169 + t84;
t114 = t181 * t133;
t115 = t181 * t134;
t22 = -t120 * t50 + t121 * t51 + qJD(1) + (-t114 * t133 - t115 * t134) * qJD(3);
t147 = t22 * t174;
t146 = t153 * t136;
t145 = t116 * t195 + t120 * (Icges(5,5) * t106 + Icges(5,6) * t107) + t121 * (Icges(5,5) * t104 - Icges(5,6) * t105);
t142 = (t104 * t53 + t105 * t54 + t133 * t150 + t64 * t87 + t65 * t89) * t136 + (t136 * t177 + t223) * qJD(3);
t33 = -t138 * t169 + t224;
t140 = ((qJD(3) * t169 + t52) * t136 + (qJD(3) * t85 - t135 * t53 + t137 * t54 + (-t135 * t89 - t137 * t87) * qJD(4)) * t138) * t136 + (t136 * t175 + t33 * t138) * qJD(3);
t130 = qJD(2) * t133;
t122 = t180 * qJD(3);
t101 = t134 * t148;
t100 = t133 * t148;
t81 = (-qJD(2) - t243) * t134;
t80 = t130 + t98;
t73 = t89 * t134;
t72 = t89 * t133;
t71 = t87 * t134;
t70 = t87 * t133;
t63 = rSges(5,1) * t106 + rSges(5,2) * t107;
t62 = rSges(5,1) * t104 - rSges(5,2) * t105;
t43 = t172 * qJD(3);
t41 = rSges(5,1) * t67 + rSges(5,2) * t66 - rSges(5,3) * t192;
t40 = rSges(5,1) * t65 + rSges(5,2) * t64 + rSges(5,3) * t193;
t30 = t50 * t195 - t121 * t91 + (-qJD(2) - t148) * t134;
t29 = t127 * t199 - t195 * t51 + t130 + t246;
t17 = t40 * t195 - t121 * t55 + (t123 * t134 + (t138 * t50 - t204 * t91) * qJD(4)) * qJD(3);
t16 = -t41 * t195 + t120 * t55 + (-t123 * t133 + (-t136 * t226 - t138 * t51) * qJD(4)) * qJD(3);
t12 = -t120 * t40 + t121 * t41 + (-t100 * t133 - t101 * t134 + t174 * t195) * qJD(3);
t11 = t120 * t24 + t121 * t23 + t195 * t33;
t8 = t104 * t37 + t105 * t39 + t133 * t151 + t64 * t47 + t65 * t49;
t7 = t104 * t36 + t105 * t38 + t133 * t152 + t64 * t46 + t65 * t48;
t6 = t120 * t21 + t121 * t20 + t195 * t32;
t5 = t120 * t19 + t121 * t18 + t195 * t31;
t4 = (qJD(3) * t170 + t35) * t136 + (qJD(3) * t45 - t135 * t37 + t137 * t39 + (-t135 * t49 - t137 * t47) * qJD(4)) * t138;
t3 = (qJD(3) * t171 + t34) * t136 + (qJD(3) * t44 - t135 * t36 + t137 * t38 + (-t135 * t48 - t137 * t46) * qJD(4)) * t138;
t1 = qJD(4) * t142 + t8 * t120 + t7 * t121;
t2 = [m(4) * t43 + m(5) * t12; m(5) * (t133 * t16 - t134 * t17) - m(4) * t122 * t185; (t133 * t19 + t134 * t18) * t157 + (t133 * t21 + t134 * t20) * t156 + (t133 * t24 + t134 * t23) * t182 + (t133 * t8 + t134 * t7) * t233 + (t133 * (-t133 * t242 - t134 * t244) + t134 * (-t133 * t241 + t134 * t92)) * t198 + (t133 * (-t133 * t244 + t134 * t242) + t134 * (t133 * t92 + t134 * t241)) * t199 - (-t131 * t244 + (t239 * t134 + (t108 - t240) * t133) * t198) * t199 / 0.2e1 - t11 * t194 / 0.2e1 + (t132 * t92 + (t240 * t133 + (-t109 - t239) * t134) * t199) * t189 + (t10 * t133 + t134 * t9) * t235 + ((t104 * t70 + t105 * t72) * t121 + (-t104 * t71 - t105 * t73) * t120 + (t223 + (t104 * t86 + t105 * t88 - t227) * t136) * qJD(4) + (((t18 + t224) * qJD(4) - t178) * t136 + (-t153 * t195 - t139) * t138) * t133) * t234 + ((t106 * t70 - t107 * t72) * t121 + (-t106 * t71 + t107 * t73) * t120 + (t222 + (t106 * t86 - t107 * t88 + t228) * t136) * qJD(4) + (((-t21 - t224) * qJD(4) + t178) * t136 + (qJD(4) * t146 + t139) * t138) * t134) * t236 + t134 * t1 / 0.2e1 + t133 * t237 + (t133 * t5 + ((-t135 * t70 + t137 * t72 + t44) * t121 + (t135 * t71 - t137 * t73 + t45) * t120 + t33 * qJD(4)) * t138 + ((t146 + (-t135 * t86 + t137 * t88 + t85) * t138 + t175) * qJD(4) + t139) * t136) * t187 + (t133 * t4 + (t3 + t6) * t134) * t186 + (-t29 * (t120 * t90 - t181 * t199) - t30 * (-t121 * t90 + t181 * t198) - t22 * (-t121 * t226 - t133 * t246 - t245 * t148) - ((-t29 * t51 + t30 * t50) * t138 + t147 * t136) * qJD(4) + (-t17 * t214 - t30 * t215 + t12 * (-t115 + t51) + t22 * (-t101 + t41)) * t134 + (t16 * t214 + t29 * t215 + t12 * (-t114 - t50) + t22 * (-t100 - t40)) * t133) * m(5) + ((-t81 * t198 + t80 * t199) * t180 - t43 * t173 - (t125 * t185 + t80 * t133 - t81 * t134) * t122 + (t245 * t243 + t172) * (-qJD(3) * t173 + qJD(1))) * m(4); t5 * t184 - t1 * t203 / 0.2e1 + (t31 * t136 - t138 * t177) * t157 + ((-t133 * t7 + t134 * t8) * t138 + t142) * t233 + t6 * t183 + t202 * t237 + (t32 * t136 - t138 * t176) * t156 + ((t10 * t134 - t133 * t9) * t138 + t141) * t235 + t11 * t188 + t136 * (t140 * qJD(4) + t4 * t120 + t3 * t121) / 0.2e1 + (t33 * t136 - t138 * t175) * t182 + ((-t133 * t3 + t134 * t4) * t138 + t140) * t186 + (t104 * t143 - t105 * t238 - t145 * t203) * t234 + (t143 * t106 + t107 * t238 + t145 * t202) * t236 + (t145 * t136 + (-t135 * t143 - t137 * t238) * t138) * t187 + ((-t16 * t51 + t17 * t50 - t29 * t41 + t30 * t40 + (t147 + (-t133 * t30 - t134 * t29) * t91) * qJD(3)) * t136 + (t29 * (-qJD(3) * t51 + t134 * t55) + t30 * (qJD(3) * t50 + t133 * t55) - t12 * t174 + t22 * (-t133 * t41 - t134 * t40) + (t133 * t17 + t134 * t16) * t91) * t138 - t29 * (t119 * t120 - t63 * t195) - t30 * (-t119 * t121 + t62 * t195) - t22 * (-t120 * t62 + t121 * t63)) * m(5);];
tauc = t2(:);
