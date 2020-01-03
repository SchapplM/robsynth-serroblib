% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRPR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR5_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR5_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:59
% EndTime: 2019-12-31 16:23:05
% DurationCPUTime: 3.81s
% Computational Cost: add. (7461->292), mult. (10269->478), div. (0->0), fcn. (11193->8), ass. (0->183)
t253 = Icges(4,1) - Icges(4,2);
t252 = -2 * Icges(4,4);
t141 = sin(pkin(6));
t140 = qJ(2) + pkin(7);
t136 = sin(t140);
t137 = cos(t140);
t144 = sin(qJ(4));
t146 = cos(qJ(4));
t179 = rSges(5,1) * t146 - rSges(5,2) * t144;
t94 = -rSges(5,3) * t137 + t179 * t136;
t84 = t141 * t94;
t142 = cos(pkin(6));
t85 = t142 * t94;
t138 = t141 ^ 2;
t139 = t142 ^ 2;
t242 = t138 + t139;
t145 = sin(qJ(2));
t147 = cos(qJ(2));
t251 = 0.2e1 * (Icges(3,1) - Icges(3,2)) * t147 * t145 + (-0.2e1 * t145 ^ 2 + 0.2e1 * t147 ^ 2) * Icges(3,4);
t197 = t136 * t141;
t163 = -Icges(3,5) * t145 - Icges(3,6) * t147;
t250 = t163 * t141;
t249 = t163 * t142;
t248 = t136 * t252 + t253 * t137;
t247 = -t253 * t136 + t137 * t252;
t204 = Icges(5,4) * t144;
t169 = Icges(5,1) * t146 - t204;
t92 = -Icges(5,5) * t137 + t169 * t136;
t212 = (-Icges(5,2) * t146 - t204) * t136 + t92;
t203 = Icges(5,4) * t146;
t164 = -Icges(5,2) * t144 + t203;
t90 = -Icges(5,6) * t137 + t164 * t136;
t211 = (-Icges(5,1) * t144 - t203) * t136 - t90;
t199 = (-Icges(5,5) * t144 - Icges(5,6) * t146) * t136 * t137;
t162 = -Icges(4,5) * t136 - Icges(4,6) * t137;
t244 = t162 * t141 + t250;
t243 = t162 * t142 + t249;
t161 = Icges(5,5) * t146 - Icges(5,6) * t144;
t88 = -Icges(5,3) * t137 + t161 * t136;
t240 = (-Icges(4,5) * t141 - t248 * t142) * t197 + (-Icges(4,6) * t141 + t247 * t142) * t141 * t137;
t239 = (-Icges(4,5) * t142 + t248 * t141) * t136 - (Icges(4,6) * t142 + t247 * t141) * t137 + t249 + t251 * t141;
t236 = 2 * qJD(2);
t235 = m(4) / 0.2e1;
t234 = m(5) / 0.2e1;
t233 = -t137 / 0.2e1;
t232 = t141 / 0.2e1;
t231 = -t142 / 0.2e1;
t230 = t142 / 0.2e1;
t229 = pkin(2) * t145;
t228 = pkin(2) * t147;
t192 = t142 * t146;
t195 = t141 * t144;
t118 = -t137 * t195 - t192;
t193 = t142 * t144;
t194 = t141 * t146;
t119 = t137 * t194 - t193;
t206 = Icges(5,4) * t119;
t63 = Icges(5,2) * t118 + Icges(5,6) * t197 + t206;
t225 = Icges(5,1) * t118 - t206 - t63;
t120 = -t137 * t193 + t194;
t121 = t137 * t192 + t195;
t205 = Icges(5,4) * t121;
t196 = t136 * t142;
t64 = Icges(5,2) * t120 + Icges(5,6) * t196 + t205;
t224 = Icges(5,1) * t120 - t205 - t64;
t112 = Icges(5,4) * t118;
t65 = Icges(5,1) * t119 + Icges(5,5) * t197 + t112;
t223 = -Icges(5,2) * t119 + t112 + t65;
t113 = Icges(5,4) * t120;
t66 = Icges(5,1) * t121 + Icges(5,5) * t196 + t113;
t222 = -Icges(5,2) * t121 + t113 + t66;
t221 = t242 * t228;
t220 = m(5) * qJD(4);
t61 = Icges(5,5) * t119 + Icges(5,6) * t118 + Icges(5,3) * t197;
t218 = t137 * t61;
t62 = Icges(5,5) * t121 + Icges(5,6) * t120 + Icges(5,3) * t196;
t217 = t137 * t62;
t216 = t137 * t88;
t37 = t118 * t64 + t119 * t66 + t197 * t62;
t215 = t37 * t142;
t38 = t120 * t63 + t121 * t65 + t196 * t61;
t214 = t38 * t141;
t67 = rSges(5,1) * t119 + rSges(5,2) * t118 + rSges(5,3) * t197;
t95 = rSges(5,3) * t136 + t137 * t179;
t42 = (t141 * t95 - t67) * t136;
t68 = rSges(5,1) * t121 + rSges(5,2) * t120 + rSges(5,3) * t196;
t43 = (-t142 * t95 + t68) * t136;
t18 = t141 * t42 - t142 * t43;
t191 = t18 * qJD(3);
t177 = -t141 * t68 + t142 * t67;
t35 = t177 * t137 + (t141 * t85 - t142 * t84) * t136;
t75 = rSges(5,1) * t118 - rSges(5,2) * t119;
t76 = rSges(5,1) * t120 - rSges(5,2) * t121;
t53 = t141 * t75 + t142 * t76;
t27 = 0.2e1 * (t35 / 0.4e1 - t53 / 0.4e1) * m(5);
t190 = t27 * qJD(1);
t188 = qJD(4) * t136;
t187 = t18 * t234;
t128 = rSges(4,1) * t136 + rSges(4,2) * t137;
t186 = -t128 - t229;
t129 = rSges(4,1) * t137 - rSges(4,2) * t136;
t185 = -t129 - t228;
t130 = pkin(3) * t136 - pkin(5) * t137;
t182 = -t130 - t94 - t229;
t131 = pkin(3) * t137 + pkin(5) * t136;
t181 = -t131 - t95 - t228;
t180 = t242 * t229;
t132 = t145 * rSges(3,1) + rSges(3,2) * t147;
t176 = -t144 * t63 + t146 * t65;
t40 = t136 * t176 - t218;
t175 = -t144 * t64 + t146 * t66;
t41 = t136 * t175 - t217;
t178 = t40 * t141 + t41 * t142;
t174 = -t144 * t90 + t146 * t92;
t160 = -t141 * t88 - t176;
t159 = -t142 * t88 - t175;
t69 = Icges(5,5) * t118 - Icges(5,6) * t119;
t23 = t118 * t223 + t119 * t225 + t197 * t69;
t70 = Icges(5,5) * t120 - Icges(5,6) * t121;
t24 = t118 * t222 + t119 * t224 + t197 * t70;
t11 = t141 * t24 - t142 * t23;
t25 = t120 * t223 + t121 * t225 + t196 * t69;
t26 = t120 * t222 + t121 * t224 + t196 * t70;
t12 = t141 * t26 - t142 * t25;
t158 = t11 * t231 + t12 * t232;
t56 = -t242 * t128 - t180;
t157 = t56 * t129;
t156 = (Icges(5,3) * t136 + t137 * t161 - t174) * t137;
t154 = -t251 * t142 + t250;
t149 = t160 * t136 + t218;
t148 = t159 * t136 + t217;
t117 = (-rSges(5,1) * t144 - rSges(5,2) * t146) * t136;
t97 = t185 * t142;
t96 = t185 * t141;
t93 = Icges(5,5) * t136 + t137 * t169;
t91 = Icges(5,6) * t136 + t137 * t164;
t83 = t92 * t142;
t82 = t92 * t141;
t81 = t90 * t142;
t80 = t90 * t141;
t77 = t242 * t132;
t60 = t181 * t142;
t59 = t182 * t142;
t58 = t181 * t141;
t57 = t182 * t141;
t55 = -t117 * t196 - t137 * t76;
t54 = t117 * t197 + t137 * t75;
t52 = -t137 * t68 - t196 * t94;
t51 = t137 * t67 + t197 * t94;
t49 = (-t141 * t76 + t142 * t75) * t136;
t48 = t136 * t174 - t216;
t47 = t177 * t136;
t46 = t120 * t90 + t121 * t92 + t196 * t88;
t45 = t118 * t90 + t119 * t92 + t197 * t88;
t44 = -t180 + (-t130 * t142 - t85) * t142 + (-t130 * t141 - t84) * t141;
t39 = t120 * t64 + t121 * t66 + t196 * t62;
t36 = t118 * t63 + t119 * t65 + t197 * t61;
t33 = (t131 * t142 + t68) * t142 + (t131 * t141 + t67) * t141 + t221;
t32 = -t137 * t70 + (-t144 * t222 + t146 * t224) * t136;
t31 = -t137 * t69 + (-t144 * t223 + t146 * t225) * t136;
t30 = -t159 * t137 + (t144 * t81 - t146 * t83 + t62) * t136;
t29 = -t160 * t137 + (t144 * t80 - t146 * t82 + t61) * t136;
t28 = (t35 + t53) * t234;
t22 = -t120 * t81 - t121 * t83 + t142 * t148;
t21 = -t120 * t80 - t121 * t82 + t142 * t149;
t20 = -t118 * t81 - t119 * t83 + t141 * t148;
t19 = -t118 * t80 - t119 * t82 + t141 * t149;
t17 = qJD(2) * t187;
t16 = t33 * t53 + (-t141 * t57 - t142 * t59) * t117;
t15 = t136 * t178 - t137 * t48;
t14 = -t137 * t46 + (t142 * t39 + t214) * t136;
t13 = -t137 * t45 + (t141 * t36 + t215) * t136;
t10 = t141 * t22 - t142 * t21;
t9 = t141 * t20 - t142 * t19;
t8 = -(t120 * t212 + t121 * t211) * t137 + (t25 * t141 + (t26 - t199) * t142) * t136;
t7 = -(t118 * t212 + t119 * t211) * t137 + (t24 * t142 + (t23 - t199) * t141) * t136;
t6 = t35 * t47 + t42 * t51 + t43 * t52;
t5 = (t156 + t178) * t137 + (t30 * t142 + t29 * t141 - (-t144 * t91 + t146 * t93 + t88) * t137 + t48) * t136;
t4 = (-t120 * t91 - t121 * t93 + t214 + (t39 - t216) * t142) * t137 + (t21 * t141 + t46 + (t22 - t156) * t142) * t136;
t3 = (-t118 * t91 - t119 * t93 + t215 + (t36 - t216) * t141) * t137 + (t20 * t142 + t45 + (t19 - t156) * t141) * t136;
t2 = m(5) * t16 + t158;
t1 = m(5) * t6 + (t14 * t230 + t13 * t232 - t5 / 0.2e1) * t137 + (t4 * t230 + t3 * t232 + t15 / 0.2e1) * t136;
t34 = [0, t28 * qJD(4) + (-m(3) * t77 / 0.2e1 + t56 * t235 + t44 * t234) * t236, 0, t28 * qJD(2) + t220 * t49; -t27 * qJD(4), t2 * qJD(4) + (m(4) * (t221 * t56 + (t142 * t157 + t186 * t97) * t142 + (t141 * t157 + t186 * t96) * t141) + m(5) * (t33 * t44 + t57 * t58 + t59 * t60) + m(3) * (t132 - t77) * t242 * (rSges(3,1) * t147 - t145 * rSges(3,2)) + (t10 + (t239 * t142 + (t154 - t244) * t141 + t240) * t142 + t243 * t138) * t232 + (t9 + (t154 * t141 + (t239 - t243) * t142 + t240) * t141 + t244 * t139) * t231) * qJD(2), -t18 * t220 / 0.2e1, -t190 + t2 * qJD(2) - t191 * m(5) / 0.2e1 + (-t15 / 0.2e1 + (t12 / 0.2e1 - t4 / 0.2e1) * t142 + (t11 / 0.2e1 - t3 / 0.2e1) * t141) * t188 + (t231 * t7 + t232 * t8 + (t5 / 0.2e1 + (t31 / 0.2e1 - t14 / 0.2e1) * t142 + (-t32 / 0.2e1 - t13 / 0.2e1) * t141) * t137 + (t49 * t33 + t47 * t53 + t54 * t59 + t55 * t57 + (-t141 * t52 - t142 * t51) * t117 - t6) * m(5)) * qJD(4); 0, ((t141 * t97 - t142 * t96) * t235 + (t141 * t60 - t142 * t58) * t234) * t236 + qJD(4) * t187, 0, t17 + (t141 * t54 - t142 * t55) * t220; t27 * qJD(2), t190 + (t136 * (t141 * t41 - t142 * t40) / 0.2e1 + (t141 * t30 - t142 * t29) * t233 + t10 * t196 / 0.2e1 + t9 * t197 / 0.2e1 + t4 * t232 + t3 * t231 - t158 + ((t141 * t39 - t142 * t38) * t230 + (t141 * t37 - t142 * t36) * t232) * t137) * qJD(2) + t1 * qJD(4) + ((t33 * t35 + t42 * t59 + t43 * t57 + t44 * t47 + t51 * t60 + t52 * t58 - t16) * qJD(2) + t191 / 0.2e1) * m(5), t17, t1 * qJD(2) + (m(5) * (t47 * t49 + t51 * t54 + t52 * t55) - t137 ^ 2 * t199 / 0.2e1) * qJD(4) + (t8 * t230 + t7 * t232 + (t32 * t142 + t31 * t141 - (-t212 * t144 + t211 * t146) * t137) * t233) * t188;];
Cq = t34;
