% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPPRR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:51
% EndTime: 2019-12-05 14:58:02
% DurationCPUTime: 3.91s
% Computational Cost: add. (17314->369), mult. (24269->618), div. (0->0), fcn. (29467->8), ass. (0->188)
t236 = -2 * Icges(5,4);
t235 = 2 * Icges(5,4);
t237 = -Icges(5,1) + Icges(5,2);
t166 = cos(pkin(8));
t218 = -t166 / 0.2e1;
t163 = pkin(9) + qJ(4);
t161 = sin(t163);
t162 = cos(t163);
t164 = sin(pkin(8));
t234 = Icges(5,6) * t166 + (t237 * t161 + t162 * t236) * t164;
t233 = Icges(5,5) * t166 + (t161 * t235 + t237 * t162) * t164;
t167 = cos(pkin(7));
t165 = sin(pkin(7));
t194 = t165 * t166;
t145 = t161 * t194 + t162 * t167;
t191 = t167 * t161;
t146 = t162 * t194 - t191;
t147 = -t165 * t162 + t166 * t191;
t192 = t166 * t167;
t148 = t161 * t165 + t162 * t192;
t197 = t164 * t167;
t198 = t164 * t165;
t232 = (-Icges(5,6) * t197 + t237 * t147 + t148 * t236) * t167 + (-Icges(5,6) * t198 + t237 * t145 + t146 * t236) * t165;
t231 = (-Icges(5,5) * t197 + t147 * t235 + t237 * t148) * t167 + (-Icges(5,5) * t198 + t145 * t235 + t237 * t146) * t165;
t228 = 2 * qJD(4);
t227 = m(5) / 0.2e1;
t226 = m(6) / 0.2e1;
t168 = sin(qJ(5));
t169 = cos(qJ(5));
t195 = t164 * t169;
t138 = -t148 * t168 + t167 * t195;
t196 = t164 * t168;
t139 = t148 * t169 + t167 * t196;
t85 = rSges(6,1) * t139 + rSges(6,2) * t138 + rSges(6,3) * t147;
t207 = -pkin(4) * t148 - pkin(6) * t147 - t85;
t136 = -t146 * t168 + t165 * t195;
t137 = t146 * t169 + t165 * t196;
t84 = rSges(6,1) * t137 + rSges(6,2) * t136 + rSges(6,3) * t145;
t209 = pkin(4) * t146 + pkin(6) * t145 + t84;
t57 = (t165 * t207 + t167 * t209) * t164;
t156 = -t162 * t196 - t166 * t169;
t157 = t162 * t195 - t166 * t168;
t199 = t161 * t164;
t113 = rSges(6,1) * t157 + rSges(6,2) * t156 + rSges(6,3) * t199;
t188 = t113 + (pkin(4) * t162 + pkin(6) * t161) * t164;
t63 = t166 * t209 + t188 * t198;
t64 = t166 * t207 - t188 * t197;
t102 = rSges(6,1) * t136 - rSges(6,2) * t137;
t103 = rSges(6,1) * t138 - rSges(6,2) * t139;
t72 = (t102 * t167 - t103 * t165) * t164;
t129 = rSges(6,1) * t156 - rSges(6,2) * t157;
t75 = t102 * t166 + t129 * t198;
t76 = -t103 * t166 - t129 * t197;
t225 = m(6) * (t57 * t72 + t63 * t75 + t64 * t76);
t224 = t145 / 0.2e1;
t223 = t146 / 0.2e1;
t222 = t147 / 0.2e1;
t221 = t148 / 0.2e1;
t220 = t161 / 0.2e1;
t219 = t165 / 0.2e1;
t217 = t167 / 0.2e1;
t130 = Icges(6,4) * t136;
t82 = Icges(6,1) * t137 + Icges(6,5) * t145 + t130;
t216 = -Icges(6,2) * t137 + t130 + t82;
t131 = Icges(6,4) * t138;
t83 = Icges(6,1) * t139 + Icges(6,5) * t147 + t131;
t215 = -Icges(6,2) * t139 + t131 + t83;
t214 = m(6) * qJD(5);
t78 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t145;
t202 = Icges(6,4) * t137;
t80 = Icges(6,2) * t136 + Icges(6,6) * t145 + t202;
t52 = t156 * t80 + t157 * t82 + t199 * t78;
t79 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t147;
t201 = Icges(6,4) * t139;
t81 = Icges(6,2) * t138 + Icges(6,6) * t147 + t201;
t53 = t156 * t81 + t157 * t83 + t199 * t79;
t110 = Icges(6,5) * t157 + Icges(6,6) * t156 + Icges(6,3) * t199;
t200 = Icges(6,4) * t157;
t111 = Icges(6,2) * t156 + Icges(6,6) * t199 + t200;
t152 = Icges(6,4) * t156;
t112 = Icges(6,1) * t157 + Icges(6,5) * t199 + t152;
t65 = t110 * t199 + t111 * t156 + t112 * t157;
t213 = t162 * (t145 * t52 + t147 * t53 + t199 * t65);
t212 = Icges(6,1) * t136 - t202 - t80;
t211 = Icges(6,1) * t138 - t201 - t81;
t184 = -rSges(6,1) * t169 + rSges(6,2) * t168;
t94 = t146 * rSges(6,3) + t145 * t184;
t210 = -pkin(4) * t145 + pkin(6) * t146 + t94;
t95 = t148 * rSges(6,3) + t147 * t184;
t208 = pkin(4) * t147 - pkin(6) * t148 - t95;
t193 = t165 * t167;
t190 = Icges(6,1) * t156 - t111 - t200;
t189 = -Icges(6,2) * t157 + t112 + t152;
t135 = (rSges(6,3) * t162 + t161 * t184) * t164;
t187 = t135 + (-pkin(4) * t161 + pkin(6) * t162) * t164;
t186 = qJD(4) * t164;
t96 = Icges(6,5) * t136 - Icges(6,6) * t137;
t38 = t156 * t216 + t157 * t212 + t199 * t96;
t97 = Icges(6,5) * t138 - Icges(6,6) * t139;
t39 = t156 * t215 + t157 * t211 + t199 * t97;
t126 = Icges(6,5) * t156 - Icges(6,6) * t157;
t56 = t126 * t199 + t156 * t189 + t157 * t190;
t17 = -t166 * t56 + (t165 * t38 + t167 * t39) * t164;
t180 = -Icges(6,5) * t169 + Icges(6,6) * t168;
t179 = Icges(6,3) * t146 + t145 * t180 + t168 * t80 - t169 * t82;
t181 = -Icges(6,4) * t169 + Icges(6,2) * t168;
t90 = Icges(6,6) * t146 + t145 * t181;
t182 = -Icges(6,1) * t169 + Icges(6,4) * t168;
t92 = Icges(6,5) * t146 + t145 * t182;
t32 = t156 * t90 + t157 * t92 + (t161 * t179 + t162 * t78) * t164;
t178 = Icges(6,3) * t148 + t147 * t180 + t168 * t81 - t169 * t83;
t91 = Icges(6,6) * t148 + t147 * t181;
t93 = Icges(6,5) * t148 + t147 * t182;
t33 = t156 * t91 + t157 * t93 + (t161 * t178 + t162 * t79) * t164;
t133 = (Icges(6,6) * t162 + t161 * t181) * t164;
t134 = (Icges(6,5) * t162 + t161 * t182) * t164;
t175 = t111 * t168 - t112 * t169 + (Icges(6,3) * t162 + t161 * t180) * t164;
t51 = t156 * t133 + t157 * t134 + (t110 * t162 + t161 * t175) * t164;
t5 = t145 * t32 + t146 * t52 + t147 * t33 + t148 * t53 + (t161 * t51 + t162 * t65) * t164;
t185 = t17 / 0.2e1 - t5 / 0.2e1;
t120 = -rSges(5,1) * t145 - rSges(5,2) * t146;
t153 = (-rSges(5,1) * t161 - rSges(5,2) * t162) * t164;
t86 = t120 * t166 + t153 * t198;
t121 = -rSges(5,1) * t147 - rSges(5,2) * t148;
t87 = -t121 * t166 - t153 * t197;
t183 = t165 * t86 - t167 * t87;
t34 = t136 * t216 + t137 * t212 + t145 * t96;
t35 = t136 * t215 + t137 * t211 + t145 * t97;
t45 = t126 * t145 + t136 * t189 + t137 * t190;
t11 = -t166 * t45 + (t165 * t34 + t167 * t35) * t164;
t36 = t138 * t216 + t139 * t212 + t147 * t96;
t37 = t138 * t215 + t139 * t211 + t147 * t97;
t46 = t126 * t147 + t138 * t189 + t139 * t190;
t12 = -t166 * t46 + (t165 * t36 + t167 * t37) * t164;
t177 = t11 * t219 + t12 * t217;
t176 = m(6) * (t165 * t75 - t167 * t76);
t54 = t113 * t146 + t135 * t145 + (-t161 * t94 - t162 * t84) * t164;
t55 = -t113 * t148 - t135 * t147 + (t161 * t95 + t162 * t85) * t164;
t174 = (t165 * t54 - t167 * t55) * t226;
t42 = -t145 * t95 - t146 * t85 + t147 * t94 + t148 * t84;
t171 = (-t166 * t42 + (t165 * t55 + t167 * t54) * t164) * t226;
t172 = m(6) * (-t166 * t72 + (t165 * t76 + t167 * t75) * t164);
t14 = t171 - t172 / 0.2e1;
t23 = t174 - t176 / 0.2e1;
t30 = 0.2e1 * (t42 / 0.4e1 - t72 / 0.4e1) * m(6);
t173 = -t30 * qJD(1) - t23 * qJD(2) - t14 * qJD(3);
t25 = t136 * t90 + t137 * t92 + t145 * t179 + t146 * t78;
t26 = t136 * t91 + t137 * t93 + t145 * t178 + t146 * t79;
t43 = t146 * t110 + t136 * t133 + t137 * t134 + t145 * t175;
t47 = t136 * t80 + t137 * t82 + t145 * t78;
t48 = t136 * t81 + t137 * t83 + t145 * t79;
t60 = t110 * t145 + t111 * t136 + t112 * t137;
t3 = t145 * t25 + t146 * t47 + t147 * t26 + t148 * t48 + (t161 * t43 + t162 * t60) * t164;
t27 = t138 * t90 + t139 * t92 + t147 * t179 + t148 * t78;
t28 = t138 * t91 + t139 * t93 + t147 * t178 + t148 * t79;
t44 = t148 * t110 + t138 * t133 + t139 * t134 + t147 * t175;
t49 = t138 * t80 + t139 * t82 + t147 * t78;
t50 = t138 * t81 + t139 * t83 + t147 * t79;
t61 = t110 * t147 + t111 * t138 + t112 * t139;
t4 = t145 * t27 + t146 * t49 + t147 * t28 + t148 * t50 + (t161 * t44 + t162 * t61) * t164;
t62 = -t145 * t85 + t147 * t84;
t69 = t113 * t145 - t199 * t84;
t70 = -t113 * t147 + t199 * t85;
t170 = m(6) * (t42 * t62 + t54 * t69 + t55 * t70) + (t145 * t47 + t147 * t48 + t199 * t60) * t223 + (t145 * t49 + t147 * t50 + t199 * t61) * t221 + t222 * t4 + t224 * t3;
t149 = (-Icges(5,5) * t161 - Icges(5,6) * t162) * t164;
t115 = -Icges(5,5) * t147 - Icges(5,6) * t148;
t114 = -Icges(5,5) * t145 - Icges(5,6) * t146;
t109 = rSges(5,1) * t148 - rSges(5,2) * t147 + rSges(5,3) * t197;
t108 = rSges(5,1) * t146 - rSges(5,2) * t145 + rSges(5,3) * t198;
t77 = (t120 * t167 - t121 * t165) * t164;
t74 = t103 * t199 - t129 * t147;
t73 = -t102 * t199 + t129 * t145;
t68 = t102 * t147 - t103 * t145;
t67 = t166 * t208 - t187 * t197;
t66 = t166 * t210 + t187 * t198;
t58 = (t165 * t208 + t167 * t210) * t164;
t31 = (t42 + t72) * t226;
t24 = t176 / 0.2e1 + t174;
t16 = t145 * t38 + t147 * t39 + t199 * t56;
t15 = t172 / 0.2e1 + t171;
t13 = -t166 * t51 + (t165 * t32 + t167 * t33) * t164;
t10 = t145 * t36 + t147 * t37 + t199 * t46;
t9 = t145 * t34 + t147 * t35 + t199 * t45;
t7 = -t166 * t44 + (t165 * t27 + t167 * t28) * t164;
t6 = -t166 * t43 + (t165 * t25 + t167 * t26) * t164;
t2 = t164 * t177 + t17 * t218 + t225;
t1 = (t213 / 0.2e1 + t5 * t220) * t164 + t170;
t8 = [0, 0, 0, t31 * qJD(5) + (t226 * t58 + t227 * t77) * t228, t31 * qJD(4) + t214 * t68; 0, 0, 0, t24 * qJD(5) + (t183 * t227 + (t165 * t66 - t167 * t67) * t226) * t228, t24 * qJD(4) + (t165 * t73 - t167 * t74) * t214; 0, 0, 0, t15 * qJD(5) + ((-t166 * t77 + (t165 * t87 + t167 * t86) * t164) * t227 + (-t166 * t58 + (t165 * t67 + t167 * t66) * t164) * t226) * t228, t15 * qJD(4) + (-t166 * t68 + (t165 * t74 + t167 * t73) * t164) * t214; -t30 * qJD(5), -t23 * qJD(5), -t14 * qJD(5), t2 * qJD(5) + (-(t233 * t147 + t234 * t148) * t192 / 0.2e1 - (t233 * t145 + t234 * t146) * t194 / 0.2e1 + (-t115 * t192 - t114 * t194 - (t233 * t161 + t234 * t162) * t166) * t218 + t7 * t217 + t6 * t219 + ((t231 * t147 + t232 * t148 - t149 * t192) * t217 + (t231 * t145 + t232 * t146 - t149 * t194) * t219 + (t231 * t161 + t232 * t162) * t218 + ((t115 * t167 ^ 2 + t114 * t193) * t217 + (t114 * t165 ^ 2 + t115 * t193) * t219) * t164) * t164) * t186 + (m(5) * ((t108 * t86 - t109 * t87) * t166 + ((t108 * t167 - t109 * t165) * t77 + t183 * (-rSges(5,3) * t166 + (rSges(5,1) * t162 - rSges(5,2) * t161) * t164)) * t164) + m(6) * (t57 * t58 + t63 * t66 + t64 * t67) + (t149 * t166 ^ 2 + t13) * t218) * qJD(4), t2 * qJD(4) + t173 + (m(6) * (t57 * t68 + t62 * t72 + t63 * t73 + t64 * t74 + t69 * t75 + t70 * t76) + t12 * t222 + t11 * t224 + t16 * t218 - t170 + (t10 * t217 + t9 * t219 - t213 / 0.2e1 + t185 * t161) * t164) * qJD(5); t30 * qJD(4), t23 * qJD(4), t14 * qJD(4), t1 * qJD(5) + (t4 * t217 + (t165 * t49 + t167 * t50) * t221 + t3 * t219 + (t165 * t47 + t167 * t48) * t223 + t13 * t220 + (t65 * t218 + (t165 * t52 + t167 * t53) * t164 / 0.2e1) * t162 - t177) * t186 - t173 + (m(6) * (t42 * t57 + t54 * t63 + t55 * t64 + t58 * t62 + t66 * t69 + t67 * t70) + t7 * t222 + t6 * t224 - t225 + (-t148 * t61 / 0.2e1 - t146 * t60 / 0.2e1 + t185) * t166) * qJD(4), t1 * qJD(4) + (m(6) * (t62 * t68 + t69 * t73 + t70 * t74) + t10 * t222 + t9 * t224 + t16 * t199 / 0.2e1) * qJD(5);];
Cq = t8;
