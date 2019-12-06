% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPPRR2
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPPRR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:28
% EndTime: 2019-12-05 14:59:33
% DurationCPUTime: 3.40s
% Computational Cost: add. (15868->381), mult. (41779->580), div. (0->0), fcn. (53656->10), ass. (0->195)
t239 = -Icges(5,1) + Icges(5,2);
t238 = -2 * Icges(5,4);
t237 = 2 * Icges(5,4);
t170 = cos(pkin(8));
t169 = cos(pkin(9));
t213 = sin(pkin(7));
t190 = t213 * t169;
t168 = sin(pkin(9));
t214 = cos(pkin(7));
t193 = t214 * t168;
t160 = t170 * t190 - t193;
t212 = sin(pkin(8));
t187 = t213 * t212;
t220 = sin(qJ(4));
t221 = cos(qJ(4));
t152 = t160 * t220 - t221 * t187;
t153 = t160 * t221 + t220 * t187;
t191 = t213 * t168;
t192 = t214 * t169;
t159 = t170 * t191 + t192;
t171 = sin(qJ(5));
t172 = cos(qJ(5));
t136 = -t153 * t171 + t159 * t172;
t137 = t153 * t172 + t159 * t171;
t84 = rSges(6,1) * t137 + rSges(6,2) * t136 + rSges(6,3) * t152;
t236 = pkin(4) * t153 + pkin(6) * t152 + t84;
t162 = t170 * t192 + t191;
t188 = t214 * t212;
t154 = t162 * t220 - t221 * t188;
t155 = t162 * t221 + t220 * t188;
t161 = t170 * t193 - t190;
t138 = -t155 * t171 + t161 * t172;
t139 = t155 * t172 + t161 * t171;
t85 = rSges(6,1) * t139 + rSges(6,2) * t138 + rSges(6,3) * t154;
t235 = pkin(4) * t155 + pkin(6) * t154 + t85;
t195 = t168 * t212;
t186 = -rSges(6,1) * t172 + rSges(6,2) * t171;
t102 = t153 * rSges(6,3) + t152 * t186;
t234 = -pkin(4) * t152 + pkin(6) * t153 + t102;
t103 = t155 * rSges(6,3) + t154 * t186;
t233 = -pkin(4) * t154 + pkin(6) * t155 + t103;
t232 = 2 * qJD(4);
t231 = m(5) / 0.2e1;
t230 = m(6) / 0.2e1;
t229 = t152 / 0.2e1;
t228 = t153 / 0.2e1;
t227 = t154 / 0.2e1;
t226 = t155 / 0.2e1;
t225 = t159 / 0.2e1;
t224 = t161 / 0.2e1;
t194 = t169 * t212;
t163 = t170 * t221 + t220 * t194;
t223 = t163 / 0.2e1;
t164 = -t170 * t220 + t221 * t194;
t222 = t164 / 0.2e1;
t208 = Icges(6,4) * t137;
t80 = Icges(6,2) * t136 + Icges(6,6) * t152 + t208;
t219 = Icges(6,1) * t136 - t208 - t80;
t207 = Icges(6,4) * t139;
t81 = Icges(6,2) * t138 + Icges(6,6) * t154 + t207;
t218 = Icges(6,1) * t138 - t207 - t81;
t134 = Icges(6,4) * t136;
t82 = Icges(6,1) * t137 + Icges(6,5) * t152 + t134;
t217 = -Icges(6,2) * t137 + t134 + t82;
t135 = Icges(6,4) * t138;
t83 = Icges(6,1) * t139 + Icges(6,5) * t154 + t135;
t216 = -Icges(6,2) * t139 + t135 + t83;
t215 = m(6) * qJD(5);
t157 = t164 * t172 + t171 * t195;
t206 = Icges(6,4) * t157;
t205 = -Icges(5,6) * t159 + t239 * t152 + t153 * t238;
t204 = -Icges(5,6) * t161 + t239 * t154 + t155 * t238;
t203 = -Icges(5,5) * t159 + t152 * t237 + t239 * t153;
t202 = -Icges(5,5) * t161 + t154 * t237 + t239 * t155;
t156 = -t164 * t171 + t172 * t195;
t111 = Icges(6,2) * t156 + Icges(6,6) * t163 + t206;
t201 = Icges(6,1) * t156 - t111 - t206;
t151 = Icges(6,4) * t156;
t112 = Icges(6,1) * t157 + Icges(6,5) * t163 + t151;
t200 = -Icges(6,2) * t157 + t112 + t151;
t113 = rSges(6,1) * t157 + rSges(6,2) * t156 + rSges(6,3) * t163;
t199 = pkin(4) * t164 + pkin(6) * t163 + t113;
t133 = t164 * rSges(6,3) + t186 * t163;
t198 = -pkin(4) * t163 + pkin(6) * t164 + t133;
t197 = -Icges(5,6) * t195 + t239 * t163 + t164 * t238;
t196 = -Icges(5,5) * t195 + t163 * t237 + t239 * t164;
t189 = t195 / 0.2e1;
t185 = -Icges(6,1) * t172 + Icges(6,4) * t171;
t184 = -Icges(6,4) * t172 + Icges(6,2) * t171;
t183 = -Icges(6,5) * t172 + Icges(6,6) * t171;
t182 = Icges(6,3) * t153 + t152 * t183 + t171 * t80 - t172 * t82;
t181 = Icges(6,3) * t155 + t154 * t183 + t171 * t81 - t172 * t83;
t180 = Icges(6,3) * t164 + t111 * t171 - t112 * t172 + t183 * t163;
t44 = t102 * t154 - t103 * t152 - t153 * t85 + t155 * t84;
t56 = -t102 * t163 + t113 * t153 + t133 * t152 - t164 * t84;
t57 = t103 * t163 - t113 * t155 - t133 * t154 + t164 * t85;
t173 = (-t44 * t170 + t57 * t187 + t56 * t188) * t230;
t94 = rSges(6,1) * t136 - rSges(6,2) * t137;
t95 = rSges(6,1) * t138 - rSges(6,2) * t139;
t72 = -t159 * t95 + t161 * t94;
t129 = rSges(6,1) * t156 - rSges(6,2) * t157;
t75 = t159 * t129 - t195 * t94;
t76 = -t161 * t129 + t195 * t95;
t174 = m(6) * (-t72 * t170 + t76 * t187 + t75 * t188);
t14 = t173 - t174 / 0.2e1;
t176 = (t56 * t213 - t57 * t214) * t230;
t178 = m(6) * (t75 * t213 - t76 * t214);
t23 = t176 - t178 / 0.2e1;
t30 = 0.2e1 * (t44 / 0.4e1 - t72 / 0.4e1) * m(6);
t179 = -t30 * qJD(1) - t23 * qJD(2) - t14 * qJD(3);
t78 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t152;
t47 = t136 * t80 + t137 * t82 + t152 * t78;
t79 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t154;
t48 = t136 * t81 + t137 * t83 + t152 * t79;
t49 = t138 * t80 + t139 * t82 + t154 * t78;
t50 = t138 * t81 + t139 * t83 + t154 * t79;
t52 = t156 * t80 + t157 * t82 + t163 * t78;
t53 = t156 * t81 + t157 * t83 + t163 * t79;
t110 = Icges(6,5) * t157 + Icges(6,6) * t156 + Icges(6,3) * t163;
t60 = t110 * t152 + t111 * t136 + t112 * t137;
t61 = t110 * t154 + t111 * t138 + t112 * t139;
t65 = t110 * t163 + t111 * t156 + t112 * t157;
t177 = (t152 * t47 + t154 * t48 + t163 * t60) * t228 + (t152 * t49 + t154 * t50 + t163 * t61) * t226 + (t152 * t52 + t154 * t53 + t163 * t65) * t222;
t88 = Icges(6,5) * t136 - Icges(6,6) * t137;
t32 = t217 * t136 + t137 * t219 + t152 * t88;
t89 = Icges(6,5) * t138 - Icges(6,6) * t139;
t33 = t136 * t216 + t137 * t218 + t152 * t89;
t126 = Icges(6,5) * t156 - Icges(6,6) * t157;
t45 = t126 * t152 + t136 * t200 + t137 * t201;
t11 = t32 * t159 + t33 * t161 + t195 * t45;
t34 = t217 * t138 + t139 * t219 + t154 * t88;
t35 = t138 * t216 + t139 * t218 + t154 * t89;
t46 = t126 * t154 + t138 * t200 + t139 * t201;
t12 = t34 * t159 + t35 * t161 + t195 * t46;
t38 = t217 * t156 + t157 * t219 + t163 * t88;
t39 = t156 * t216 + t157 * t218 + t163 * t89;
t54 = t126 * t163 + t156 * t200 + t157 * t201;
t17 = t38 * t159 + t39 * t161 + t195 * t54;
t175 = t11 * t225 + t12 * t224 + t17 * t189;
t146 = -rSges(5,1) * t163 - rSges(5,2) * t164;
t143 = -Icges(5,5) * t163 - Icges(5,6) * t164;
t142 = t164 * rSges(5,1) - t163 * rSges(5,2) + rSges(5,3) * t195;
t132 = Icges(6,5) * t164 + t185 * t163;
t131 = Icges(6,6) * t164 + t184 * t163;
t121 = -rSges(5,1) * t154 - rSges(5,2) * t155;
t120 = -rSges(5,1) * t152 - rSges(5,2) * t153;
t115 = -Icges(5,5) * t154 - Icges(5,6) * t155;
t114 = -Icges(5,5) * t152 - Icges(5,6) * t153;
t109 = rSges(5,1) * t155 - rSges(5,2) * t154 + rSges(5,3) * t161;
t108 = rSges(5,1) * t153 - rSges(5,2) * t152 + rSges(5,3) * t159;
t101 = Icges(6,5) * t155 + t154 * t185;
t100 = Icges(6,5) * t153 + t152 * t185;
t99 = Icges(6,6) * t155 + t154 * t184;
t98 = Icges(6,6) * t153 + t152 * t184;
t87 = t121 * t195 - t161 * t146;
t86 = -t120 * t195 + t159 * t146;
t77 = t120 * t161 - t121 * t159;
t74 = -t129 * t154 + t163 * t95;
t73 = t129 * t152 - t163 * t94;
t70 = -t113 * t154 + t163 * t85;
t69 = t113 * t152 - t163 * t84;
t68 = -t152 * t95 + t154 * t94;
t67 = -t198 * t161 + t233 * t195;
t66 = t198 * t159 - t234 * t195;
t64 = -t152 * t85 + t154 * t84;
t63 = -t199 * t161 + t235 * t195;
t62 = t199 * t159 - t236 * t195;
t59 = -t233 * t159 + t234 * t161;
t55 = -t235 * t159 + t236 * t161;
t51 = t164 * t110 + t156 * t131 + t157 * t132 + t163 * t180;
t42 = t155 * t110 + t138 * t131 + t139 * t132 + t154 * t180;
t41 = t153 * t110 + t136 * t131 + t137 * t132 + t152 * t180;
t37 = t157 * t101 + t156 * t99 + t163 * t181 + t164 * t79;
t36 = t157 * t100 + t156 * t98 + t163 * t182 + t164 * t78;
t31 = (t44 + t72) * t230;
t28 = t139 * t101 + t138 * t99 + t154 * t181 + t155 * t79;
t27 = t139 * t100 + t138 * t98 + t154 * t182 + t155 * t78;
t26 = t137 * t101 + t136 * t99 + t152 * t181 + t153 * t79;
t25 = t137 * t100 + t136 * t98 + t152 * t182 + t153 * t78;
t24 = t178 / 0.2e1 + t176;
t18 = t55 * t72 + t62 * t75 + t63 * t76;
t16 = t152 * t38 + t154 * t39 + t163 * t54;
t15 = t174 / 0.2e1 + t173;
t13 = t36 * t159 + t37 * t161 + t195 * t51;
t10 = t152 * t34 + t154 * t35 + t163 * t46;
t9 = t152 * t32 + t154 * t33 + t163 * t45;
t8 = t44 * t64 + t56 * t69 + t57 * t70;
t7 = t27 * t159 + t28 * t161 + t195 * t42;
t6 = t25 * t159 + t26 * t161 + t195 * t41;
t5 = t152 * t36 + t153 * t52 + t154 * t37 + t155 * t53 + t163 * t51 + t164 * t65;
t4 = t152 * t27 + t153 * t49 + t154 * t28 + t155 * t50 + t163 * t42 + t164 * t61;
t3 = t152 * t25 + t153 * t47 + t154 * t26 + t155 * t48 + t163 * t41 + t164 * t60;
t2 = m(6) * t18 + t175;
t1 = m(6) * t8 + t5 * t223 + t4 * t227 + t3 * t229 + t177;
t19 = [0, 0, 0, t31 * qJD(5) + (t59 * t230 + t77 * t231) * t232, t31 * qJD(4) + t215 * t68; 0, 0, 0, t24 * qJD(5) + ((t213 * t86 - t214 * t87) * t231 + (t213 * t66 - t214 * t67) * t230) * t232, t24 * qJD(4) + (t213 * t73 - t214 * t74) * t215; 0, 0, 0, t15 * qJD(5) + ((-t77 * t170 + t187 * t87 + t188 * t86) * t231 + (-t59 * t170 + t187 * t67 + t188 * t66) * t230) * t232, t15 * qJD(4) + (-t68 * t170 + t187 * t74 + t188 * t73) * t215; -t30 * qJD(5), -t23 * qJD(5), -t14 * qJD(5), t2 * qJD(5) + (m(5) * ((-t108 * t195 + t159 * t142) * t86 + (t109 * t195 - t161 * t142) * t87 + (t108 * t161 - t109 * t159) * t77) + m(6) * (t55 * t59 + t62 * t66 + t63 * t67) + ((t115 * t159 + t152 * t202 + t153 * t204) * t161 + (t114 * t159 + t152 * t203 + t153 * t205) * t159 + (t143 * t159 + t152 * t196 + t153 * t197) * t195 + t6) * t225 + ((t115 * t161 + t154 * t202 + t155 * t204) * t161 + (t114 * t161 + t154 * t203 + t155 * t205) * t159 + (t143 * t161 + t154 * t196 + t155 * t197) * t195 + t7) * t224 + ((t115 * t195 + t163 * t202 + t164 * t204) * t161 + (t114 * t195 + t163 * t203 + t164 * t205) * t159 + (t143 * t195 + t163 * t196 + t164 * t197) * t195 + t13) * t189) * qJD(4), t2 * qJD(4) + (t10 * t224 + t9 * t225 + t16 * t189 + (t17 / 0.2e1 - t5 / 0.2e1) * t163 + (t12 / 0.2e1 - t4 / 0.2e1) * t154 + (t11 / 0.2e1 - t3 / 0.2e1) * t152 + (t55 * t68 + t62 * t73 + t63 * t74 + t64 * t72 + t69 * t75 + t70 * t76 - t8) * m(6) - t177) * qJD(5) + t179; t30 * qJD(4), t23 * qJD(4), t14 * qJD(4), (t4 * t224 + (t49 * t159 + t50 * t161 + t195 * t61) * t226 + t7 * t227 + t3 * t225 + (t47 * t159 + t48 * t161 + t195 * t60) * t228 + t6 * t229 + t5 * t189 + (t52 * t159 + t53 * t161 + t195 * t65) * t222 + t13 * t223 + (t44 * t55 + t56 * t62 + t57 * t63 + t59 * t64 + t66 * t69 + t67 * t70 - t18) * m(6) - t175) * qJD(4) + t1 * qJD(5) - t179, t1 * qJD(4) + (m(6) * (t64 * t68 + t69 * t73 + t70 * t74) + t10 * t227 + t9 * t229 + t16 * t223) * qJD(5);];
Cq = t19;
