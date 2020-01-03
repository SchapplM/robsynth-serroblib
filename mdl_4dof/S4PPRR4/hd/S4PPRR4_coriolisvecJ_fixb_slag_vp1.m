% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:26
% EndTime: 2019-12-31 16:18:37
% DurationCPUTime: 6.32s
% Computational Cost: add. (6445->367), mult. (9540->593), div. (0->0), fcn. (9276->6), ass. (0->187)
t133 = sin(pkin(6));
t132 = pkin(7) + qJ(3);
t130 = sin(t132);
t131 = cos(t132);
t135 = sin(qJ(4));
t136 = cos(qJ(4));
t185 = rSges(5,1) * t136 - rSges(5,2) * t135;
t82 = -t131 * rSges(5,3) + t130 * t185;
t249 = t133 * t82;
t134 = cos(pkin(6));
t248 = t134 * t82;
t237 = t134 ^ 2;
t238 = t133 ^ 2;
t240 = t237 + t238;
t124 = rSges(4,1) * t130 + rSges(4,2) * t131;
t247 = qJD(3) * t124 * t240;
t200 = qJD(3) * t131;
t201 = qJD(3) * t130;
t217 = Icges(4,4) * t131;
t218 = Icges(4,4) * t130;
t246 = -(-Icges(4,2) * t131 - t218) * t201 + (-Icges(4,1) * t130 - t217) * t200;
t170 = -Icges(4,2) * t130 + t217;
t173 = Icges(4,1) * t131 - t218;
t245 = t130 * (-Icges(4,5) * t134 + t133 * t173) + t131 * (-Icges(4,6) * t134 + t133 * t170);
t244 = -t130 * (Icges(4,5) * t133 + t134 * t173) - t131 * (Icges(4,6) * t133 + t134 * t170);
t197 = qJD(4) * t130;
t199 = qJD(3) * t133;
t122 = t134 * t197 + t199;
t234 = t122 / 0.2e1;
t198 = qJD(3) * t134;
t123 = t133 * t197 - t198;
t232 = t123 / 0.2e1;
t243 = qJD(4) / 0.2e1;
t166 = Icges(5,5) * t136 - Icges(5,6) * t135;
t76 = -Icges(5,3) * t131 + t130 * t166;
t213 = Icges(5,4) * t136;
t168 = -Icges(5,2) * t135 + t213;
t78 = -Icges(5,6) * t131 + t130 * t168;
t214 = Icges(5,4) * t135;
t171 = Icges(5,1) * t136 - t214;
t80 = -Icges(5,5) * t131 + t130 * t171;
t205 = t134 * t135;
t206 = t133 * t136;
t118 = -t131 * t205 + t206;
t208 = t130 * t134;
t204 = t134 * t136;
t207 = t133 * t135;
t119 = t131 * t204 + t207;
t215 = Icges(5,4) * t119;
t51 = Icges(5,2) * t118 + Icges(5,6) * t208 + t215;
t109 = Icges(5,4) * t118;
t53 = Icges(5,1) * t119 + Icges(5,5) * t208 + t109;
t175 = -t135 * t51 + t136 * t53;
t116 = -t131 * t207 - t204;
t209 = t130 * t133;
t117 = t131 * t206 - t205;
t216 = Icges(5,4) * t117;
t50 = Icges(5,2) * t116 + Icges(5,6) * t209 + t216;
t108 = Icges(5,4) * t116;
t52 = Icges(5,1) * t117 + Icges(5,5) * t209 + t108;
t176 = -t135 * t50 + t136 * t52;
t239 = -(-t134 * t76 - t175) * t122 - (-t133 * t76 - t176) * t123;
t113 = (-Icges(5,2) * t136 - t214) * t130;
t196 = qJD(4) * t131;
t142 = t122 * (-Icges(5,2) * t119 + t109 + t53) + t123 * (-Icges(5,2) * t117 + t108 + t52) - t196 * (t113 + t80);
t192 = t131 * t198;
t195 = t135 * t201;
t74 = -qJD(4) * t119 + t134 * t195;
t194 = t136 * t201;
t75 = qJD(4) * t118 - t134 * t194;
t35 = Icges(5,5) * t75 + Icges(5,6) * t74 + Icges(5,3) * t192;
t49 = Icges(5,5) * t119 + Icges(5,6) * t118 + Icges(5,3) * t208;
t159 = t130 * t35 + t200 * t49;
t37 = Icges(5,4) * t75 + Icges(5,2) * t74 + Icges(5,6) * t192;
t39 = Icges(5,1) * t75 + Icges(5,4) * t74 + Icges(5,5) * t192;
t10 = t118 * t37 + t119 * t39 + t134 * t159 + t51 * t74 + t53 * t75;
t112 = (-Icges(5,5) * t135 - Icges(5,6) * t136) * t130;
t77 = Icges(5,3) * t130 + t131 * t166;
t44 = qJD(3) * t77 + qJD(4) * t112;
t158 = t130 * t44 + t200 * t76;
t21 = t118 * t51 + t119 * t53 + t208 * t49;
t48 = Icges(5,5) * t117 + Icges(5,6) * t116 + Icges(5,3) * t209;
t20 = t118 * t50 + t119 * t52 + t208 * t48;
t227 = t133 * t20;
t181 = t134 * t21 + t227;
t30 = t118 * t78 + t119 * t80 + t208 * t76;
t222 = t30 * t130;
t79 = Icges(5,6) * t130 + t131 * t168;
t45 = qJD(3) * t79 + qJD(4) * t113;
t114 = (-Icges(5,1) * t135 - t213) * t130;
t81 = Icges(5,5) * t130 + t131 * t171;
t46 = qJD(3) * t81 + qJD(4) * t114;
t140 = -t131 * (t118 * t45 + t119 * t46 + t134 * t158 + t74 * t78 + t75 * t80) + (t131 * t181 + t222) * qJD(3);
t193 = t131 * t199;
t72 = -qJD(4) * t117 + t133 * t195;
t73 = qJD(4) * t116 - t133 * t194;
t34 = Icges(5,5) * t73 + Icges(5,6) * t72 + Icges(5,3) * t193;
t160 = t130 * t34 + t200 * t48;
t36 = Icges(5,4) * t73 + Icges(5,2) * t72 + Icges(5,6) * t193;
t38 = Icges(5,1) * t73 + Icges(5,4) * t72 + Icges(5,5) * t193;
t9 = t118 * t36 + t119 * t38 + t134 * t160 + t50 * t74 + t52 * t75;
t236 = t10 * t234 + t140 * t243 + t232 * t9;
t235 = -t122 / 0.2e1;
t233 = -t123 / 0.2e1;
t228 = t131 * t76;
t19 = t116 * t51 + t117 * t53 + t209 * t49;
t225 = t134 * t19;
t29 = t116 * t78 + t117 * t80 + t209 * t76;
t223 = t29 * t130;
t127 = pkin(3) * t131 + pkin(5) * t130;
t121 = qJD(3) * t127;
t115 = (-rSges(5,1) * t135 - rSges(5,2) * t136) * t130;
t83 = t130 * rSges(5,3) + t131 * t185;
t47 = qJD(3) * t83 + qJD(4) * t115;
t220 = -t121 - t47;
t126 = pkin(3) * t130 - pkin(5) * t131;
t219 = -t126 - t82;
t203 = qJD(2) * t134;
t191 = t198 / 0.2e1;
t190 = -t196 / 0.2e1;
t189 = t196 / 0.2e1;
t188 = qJD(3) * t243;
t187 = t130 * t188;
t186 = t131 * t188;
t184 = t49 * t122 + t48 * t123;
t41 = rSges(5,1) * t75 + rSges(5,2) * t74 + rSges(5,3) * t192;
t55 = rSges(5,1) * t119 + rSges(5,2) * t118 + rSges(5,3) * t208;
t16 = -t41 * t196 - t122 * t47 + (-t121 * t133 + (t130 * t55 - t131 * t248) * qJD(4)) * qJD(3);
t40 = rSges(5,1) * t73 + rSges(5,2) * t72 + rSges(5,3) * t193;
t54 = rSges(5,1) * t117 + rSges(5,2) * t116 + rSges(5,3) * t209;
t17 = t40 * t196 + t123 * t47 + (-t121 * t134 + (-t130 * t54 + t131 * t249) * qJD(4)) * qJD(3);
t183 = t133 * t17 - t134 * t16;
t18 = t116 * t50 + t117 * t52 + t209 * t48;
t182 = t133 * t18 + t225;
t23 = t130 * t176 - t131 * t48;
t24 = t130 * t175 - t131 * t49;
t180 = t23 * t133 + t24 * t134;
t179 = -t133 * t55 + t134 * t54;
t174 = -t135 * t78 + t136 * t80;
t167 = -Icges(4,5) * t130 - Icges(4,6) * t131;
t165 = t133 * t186;
t164 = t134 * t186;
t161 = -t174 + t77;
t157 = qJD(3) * t126;
t110 = t127 * t133;
t111 = t127 * t134;
t22 = t122 * t54 - t123 * t55 + qJD(1) + (t110 * t133 + t111 * t134) * qJD(3);
t156 = t22 * t179;
t152 = qJD(3) * t167;
t151 = t112 * t196 - t122 * (Icges(5,5) * t118 - Icges(5,6) * t119) - t123 * (Icges(5,5) * t116 - Icges(5,6) * t117);
t146 = t130 * t151;
t145 = -qJD(3) * t245 + t133 * t246;
t144 = qJD(3) * t244 + t134 * t246;
t143 = (Icges(5,1) * t118 - t215 - t51) * t122 + (Icges(5,1) * t116 - t216 - t50) * t123 - (t114 - t78) * t196;
t141 = -(t116 * t45 + t117 * t46 + t133 * t158 + t72 * t78 + t73 * t80) * t131 + (t131 * t182 + t223) * qJD(3);
t33 = t130 * t174 - t228;
t139 = -t131 * ((qJD(3) * t174 - t44) * t131 + (qJD(3) * t76 - t135 * t45 + t136 * t46 + (-t135 * t80 - t136 * t78) * qJD(4)) * t130) + (t33 * t130 + t131 * t180) * qJD(3);
t138 = t133 * t244 + t134 * t245;
t137 = (-t161 * t196 - t239) * t130;
t129 = qJD(2) * t133;
t103 = t167 * t134;
t102 = t167 * t133;
t101 = t134 * t157;
t100 = t133 * t157;
t93 = t134 * t152;
t92 = t133 * t152;
t89 = -t124 * t199 - t203;
t88 = -t124 * t198 + t129;
t69 = t80 * t134;
t68 = t80 * t133;
t67 = t78 * t134;
t66 = t78 * t133;
t63 = rSges(5,1) * t118 - rSges(5,2) * t119;
t62 = rSges(5,1) * t116 - rSges(5,2) * t117;
t43 = t247 * qJD(3);
t32 = -t122 * t82 - t126 * t199 - t196 * t55 - t203;
t31 = t123 * t82 - t126 * t198 + t196 * t54 + t129;
t12 = t122 * t40 - t123 * t41 + (-t100 * t133 - t101 * t134 + t179 * t196) * qJD(3);
t11 = t122 * t24 + t123 * t23 - t196 * t33;
t8 = t116 * t37 + t117 * t39 + t133 * t159 + t51 * t72 + t53 * t73;
t7 = t116 * t36 + t117 * t38 + t133 * t160 + t50 * t72 + t52 * t73;
t6 = t122 * t21 + t123 * t20 - t196 * t30;
t5 = t122 * t19 + t123 * t18 - t196 * t29;
t4 = (qJD(3) * t175 - t35) * t131 + (qJD(3) * t49 - t135 * t37 + t136 * t39 + (-t135 * t53 - t136 * t51) * qJD(4)) * t130;
t3 = (qJD(3) * t176 - t34) * t131 + (qJD(3) * t48 - t135 * t36 + t136 * t38 + (-t135 * t52 - t136 * t50) * qJD(4)) * t130;
t1 = qJD(4) * t141 + t122 * t8 + t123 * t7;
t2 = [-m(4) * t43 + m(5) * t12; m(5) * t183; (t133 * t21 - t134 * t20) * t164 + (t133 * t19 - t134 * t18) * t165 + (t133 * t24 - t134 * t23) * t187 + (t10 * t133 - t134 * t9) * t234 + (t133 * (t133 * t93 + t134 * t144) - t134 * (t133 * t92 + t134 * t145)) * t199 - (t133 * (t133 * t144 - t134 * t93) - t134 * (t133 * t145 - t134 * t92)) * t198 + (t102 * qJD(3) * t237 + (-t134 * t103 + t138) * t199) * t191 - t11 * t197 / 0.2e1 - (t103 * qJD(3) * t238 + (-t133 * t102 + t138) * t198) * t199 / 0.2e1 + (((t135 * t67 - t136 * t69 + t49) * t122 + (t135 * t66 - t136 * t68 + t48) * t123 + t33 * qJD(4)) * t130 + ((t161 * t131 + (t135 * t79 - t136 * t81 - t76) * t130 + t180) * qJD(4) + t239) * t131) * t189 + (t133 * t8 - t134 * t7) * t232 + ((-t118 * t67 - t119 * t69) * t122 + (-t118 * t66 - t119 * t68) * t123 + (t222 + (-t118 * t79 - t119 * t81 + t227) * t131) * qJD(4) + (((t21 - t228) * qJD(4) + t184) * t131 + t137) * t134) * t235 + ((-t116 * t67 - t117 * t69) * t122 + (-t116 * t66 - t117 * t68) * t123 + (t223 + (-t116 * t79 - t117 * t81 + t225) * t131) * qJD(4) + (((t18 - t228) * qJD(4) + t184) * t131 + t137) * t133) * t233 + t133 * t236 - t134 * t1 / 0.2e1 + (-t31 * (t123 * t83 - t127 * t198) - t32 * (-t122 * t83 - t127 * t199) - t22 * (-t122 * t249 + t123 * t248 - t157 * t240) - ((-t31 * t54 + t32 * t55) * t130 + t156 * t131) * qJD(4) + (t17 * t219 + t31 * t220 + t12 * (t111 + t55) + t22 * (-t101 + t41)) * t134 + (t16 * t219 + t32 * t220 + t12 * (t110 + t54) + t22 * (-t100 + t40)) * t133) * m(5) + (-t43 * t240 + t88 * t198 + t89 * t199 + (-t89 * t133 - t88 * t134 + t247) * qJD(3)) * (rSges(4,1) * t131 - rSges(4,2) * t130) * m(4) + ((-t3 + t6) * t134 + (t4 + t5) * t133) * t190; t131 * t6 * t191 + t208 * t236 + (t130 * t181 - t131 * t30) * t164 + ((t10 * t134 + t133 * t9) * t130 + t140) * t234 + t5 * t193 / 0.2e1 + t1 * t209 / 0.2e1 + (t130 * t182 - t131 * t29) * t165 + ((t133 * t7 + t134 * t8) * t130 + t141) * t232 + t11 * t201 / 0.2e1 - t131 * (t139 * qJD(4) + t122 * t4 + t123 * t3) / 0.2e1 + (t130 * t180 - t131 * t33) * t187 + ((t133 * t3 + t134 * t4) * t130 + t139) * t190 + (t118 * t142 + t119 * t143 - t134 * t146) * t235 + (t116 * t142 + t117 * t143 - t133 * t146) * t233 + (t151 * t131 + (-t135 * t142 + t143 * t136) * t130) * t189 + ((-t16 * t55 + t17 * t54 + t31 * t40 - t32 * t41 + (t156 + (t133 * t31 - t134 * t32) * t82) * qJD(3)) * t131 + (t31 * (-qJD(3) * t54 + t133 * t47) + t32 * (qJD(3) * t55 - t134 * t47) + t12 * t179 + t22 * (-t133 * t41 + t134 * t40) + t183 * t82) * t130 - t31 * (t115 * t123 + t62 * t196) - t32 * (-t115 * t122 - t63 * t196) - t22 * (t122 * t62 - t123 * t63)) * m(5);];
tauc = t2(:);
