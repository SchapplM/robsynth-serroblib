% Calculate joint inertia matrix for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:54
% EndTime: 2019-12-31 20:16:59
% DurationCPUTime: 1.86s
% Computational Cost: add. (4687->292), mult. (4413->436), div. (0->0), fcn. (4595->10), ass. (0->150)
t221 = Icges(3,3) + Icges(4,3);
t141 = qJ(2) + pkin(9);
t132 = sin(t141);
t133 = cos(t141);
t146 = sin(qJ(2));
t149 = cos(qJ(2));
t220 = Icges(3,5) * t149 + Icges(4,5) * t133 - Icges(3,6) * t146 - Icges(4,6) * t132;
t147 = sin(qJ(1));
t142 = t147 ^ 2;
t219 = t147 * pkin(6);
t134 = qJ(4) + t141;
t130 = cos(t134);
t150 = cos(qJ(1));
t191 = t130 * t150;
t129 = sin(t134);
t192 = t129 * t150;
t148 = cos(qJ(5));
t188 = t147 * t148;
t145 = sin(qJ(5));
t190 = t145 * t150;
t100 = -t130 * t190 + t188;
t187 = t148 * t150;
t189 = t147 * t145;
t101 = t130 * t187 + t189;
t53 = rSges(6,1) * t101 + rSges(6,2) * t100 + rSges(6,3) * t192;
t218 = pkin(4) * t191 + pkin(8) * t192 + t53;
t174 = rSges(4,1) * t133 - rSges(4,2) * t132;
t143 = t150 ^ 2;
t194 = Icges(5,4) * t130;
t161 = -Icges(5,2) * t129 + t194;
t77 = Icges(5,6) * t147 + t150 * t161;
t195 = Icges(5,4) * t129;
t164 = Icges(5,1) * t130 - t195;
t79 = Icges(5,5) * t147 + t150 * t164;
t171 = -t129 * t77 + t130 * t79;
t76 = -Icges(5,6) * t150 + t147 * t161;
t78 = -Icges(5,5) * t150 + t147 * t164;
t172 = t129 * t76 - t130 * t78;
t158 = Icges(5,5) * t130 - Icges(5,6) * t129;
t74 = -Icges(5,3) * t150 + t147 * t158;
t75 = Icges(5,3) * t147 + t150 * t158;
t193 = t129 * t147;
t98 = -t130 * t189 - t187;
t99 = t130 * t188 - t190;
t46 = Icges(6,5) * t99 + Icges(6,6) * t98 + Icges(6,3) * t193;
t48 = Icges(6,4) * t99 + Icges(6,2) * t98 + Icges(6,6) * t193;
t50 = Icges(6,1) * t99 + Icges(6,4) * t98 + Icges(6,5) * t193;
t14 = t193 * t46 + t48 * t98 + t50 * t99;
t47 = Icges(6,5) * t101 + Icges(6,6) * t100 + Icges(6,3) * t192;
t49 = Icges(6,4) * t101 + Icges(6,2) * t100 + Icges(6,6) * t192;
t51 = Icges(6,1) * t101 + Icges(6,4) * t100 + Icges(6,5) * t192;
t15 = t193 * t47 + t49 * t98 + t51 * t99;
t8 = -t14 * t150 + t147 * t15;
t217 = -t143 * t74 - (t171 * t147 + (t172 - t75) * t150) * t147 - t8;
t216 = -t147 * t220 + t150 * t221;
t215 = t147 * t221 + t150 * t220;
t214 = t147 / 0.2e1;
t213 = -t150 / 0.2e1;
t16 = t100 * t48 + t101 * t50 + t192 * t46;
t17 = t100 * t49 + t101 * t51 + t192 * t47;
t9 = t147 * t17 - t150 * t16;
t212 = (t142 * t75 + t9 + (t172 * t150 + (t171 - t74) * t147) * t150) * t147;
t211 = pkin(2) * t146;
t210 = pkin(4) * t130;
t144 = -qJ(3) - pkin(6);
t131 = t149 * pkin(2) + pkin(1);
t154 = rSges(5,1) * t191 - rSges(5,2) * t192 + t147 * rSges(5,3);
t173 = rSges(5,1) * t130 - rSges(5,2) * t129;
t41 = t147 * (-rSges(5,3) * t150 + t147 * t173) + t150 * t154;
t126 = t150 * t131;
t139 = t150 * pkin(6);
t209 = t147 * (t139 + (-pkin(1) + t131) * t147) + t150 * (-pkin(1) * t150 + t126 - t219);
t208 = rSges(3,1) * t149;
t206 = rSges(3,2) * t146;
t65 = -Icges(6,6) * t130 + (Icges(6,4) * t148 - Icges(6,2) * t145) * t129;
t204 = t145 * t65;
t203 = t150 * rSges(3,3);
t21 = -t130 * t46 + (-t145 * t48 + t148 * t50) * t129;
t202 = t21 * t150;
t22 = -t130 * t47 + (-t145 * t49 + t148 * t51) * t129;
t201 = t22 * t147;
t67 = -rSges(6,3) * t130 + (rSges(6,1) * t148 - rSges(6,2) * t145) * t129;
t200 = -pkin(4) * t129 + pkin(8) * t130 - t67;
t199 = Icges(3,4) * t146;
t198 = Icges(3,4) * t149;
t197 = Icges(4,4) * t132;
t196 = Icges(4,4) * t133;
t185 = t147 * rSges(3,3) + t150 * t208;
t183 = t142 + t143;
t182 = -rSges(4,1) * t132 - rSges(4,2) * t133 - t211;
t113 = pkin(3) * t133 + t131;
t108 = t150 * t113;
t181 = t150 * (t108 - t126) + t209 + (t113 - t131) * t142;
t177 = -t99 * rSges(6,1) - t98 * rSges(6,2);
t52 = rSges(6,3) * t193 - t177;
t23 = t147 * t52 + t142 * (pkin(8) * t129 + t210) + t218 * t150;
t140 = -pkin(7) + t144;
t180 = -t140 * t147 + t108;
t64 = -Icges(6,3) * t130 + (Icges(6,5) * t148 - Icges(6,6) * t145) * t129;
t66 = -Icges(6,5) * t130 + (Icges(6,1) * t148 - Icges(6,4) * t145) * t129;
t26 = t193 * t64 + t65 * t98 + t66 * t99;
t3 = -t26 * t130 + (t14 * t147 + t15 * t150) * t129;
t27 = t100 * t65 + t101 * t66 + t192 * t64;
t4 = -t27 * t130 + (t147 * t16 + t150 * t17) * t129;
t179 = t8 * t193 / 0.2e1 + t3 * t213 + t4 * t214 - t130 * (t201 - t202) / 0.2e1 + t9 * t192 / 0.2e1;
t176 = -pkin(3) * t132 - t211;
t175 = -t206 + t208;
t166 = Icges(3,1) * t149 - t199;
t165 = Icges(4,1) * t133 - t197;
t163 = -Icges(3,2) * t146 + t198;
t162 = -Icges(4,2) * t132 + t196;
t104 = Icges(5,2) * t130 + t195;
t105 = Icges(5,1) * t129 + t194;
t157 = -t104 * t129 + t105 * t130;
t156 = t150 * t217 + t212;
t155 = t147 * rSges(4,3) + t150 * t174;
t106 = rSges(5,1) * t129 + rSges(5,2) * t130;
t153 = -t106 + t176;
t152 = t176 + t200;
t103 = Icges(5,5) * t129 + Icges(5,6) * t130;
t151 = -t202 / 0.2e1 + t201 / 0.2e1 + (t103 * t147 + t129 * t79 + t130 * t77 + t150 * t157 + t27) * t214 + (-t103 * t150 + t129 * t78 + t130 * t76 + t147 * t157 + t26) * t213;
t124 = rSges(2,1) * t150 - rSges(2,2) * t147;
t123 = -rSges(2,1) * t147 - rSges(2,2) * t150;
t122 = rSges(3,1) * t146 + rSges(3,2) * t149;
t81 = t182 * t150;
t80 = t182 * t147;
t73 = t219 + (pkin(1) - t206) * t150 + t185;
t72 = t203 + t139 + (-pkin(1) - t175) * t147;
t63 = t153 * t150;
t62 = t153 * t147;
t61 = -t147 * t144 + t126 + t155;
t60 = (rSges(4,3) - t144) * t150 + (-t131 - t174) * t147;
t59 = t129 * t148 * t66;
t56 = t150 * (-t150 * t206 + t185) + (t147 * t175 - t203) * t147;
t55 = t154 + t180;
t54 = (rSges(5,3) - t140) * t150 + (-t113 - t173) * t147;
t45 = t200 * t150;
t44 = t200 * t147;
t40 = t152 * t150;
t39 = t152 * t147;
t34 = t180 + t218;
t33 = -t140 * t150 + (-t210 - t113 + (-rSges(6,3) - pkin(8)) * t129) * t147 + t177;
t32 = -t130 * t53 - t192 * t67;
t31 = t130 * t52 + t193 * t67;
t30 = t150 * t155 + (-t150 * rSges(4,3) + t147 * t174) * t147 + t209;
t29 = -t129 * t204 - t130 * t64 + t59;
t28 = (-t147 * t53 + t150 * t52) * t129;
t18 = t181 + t41;
t11 = t23 + t181;
t1 = [t133 * (Icges(4,2) * t133 + t197) + t132 * (Icges(4,1) * t132 + t196) + t149 * (Icges(3,2) * t149 + t199) + t146 * (Icges(3,1) * t146 + t198) + Icges(2,3) + t59 + (-t64 + t104) * t130 + (t105 - t204) * t129 + m(6) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t54 ^ 2 + t55 ^ 2) + m(4) * (t60 ^ 2 + t61 ^ 2) + m(3) * (t72 ^ 2 + t73 ^ 2) + m(2) * (t123 ^ 2 + t124 ^ 2); t151 + m(3) * (-t147 * t73 - t150 * t72) * t122 + m(6) * (t33 * t40 + t34 * t39) + m(5) * (t54 * t63 + t55 * t62) + m(4) * (t60 * t81 + t61 * t80) + (t132 * (Icges(4,5) * t147 + t150 * t165) + t133 * (Icges(4,6) * t147 + t150 * t162) + t146 * (Icges(3,5) * t147 + t150 * t166) + t149 * (Icges(3,6) * t147 + t150 * t163)) * t214 + (t132 * (-Icges(4,5) * t150 + t147 * t165) + t133 * (-Icges(4,6) * t150 + t147 * t162) + t146 * (-Icges(3,5) * t150 + t147 * t166) + t149 * (-Icges(3,6) * t150 + t147 * t163)) * t213 + (Icges(3,5) * t146 + Icges(4,5) * t132 + Icges(3,6) * t149 + Icges(4,6) * t133) * (t142 / 0.2e1 + t143 / 0.2e1); m(6) * (t11 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t18 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(4) * (t30 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(3) * (t122 ^ 2 * t183 + t56 ^ 2) + t212 + t215 * t147 * t142 + (t216 * t143 + (t147 * t216 + t150 * t215) * t147 + t217) * t150; m(6) * (t147 * t33 - t150 * t34) + m(5) * (t147 * t54 - t150 * t55) + m(4) * (t147 * t60 - t150 * t61); m(6) * (t147 * t40 - t150 * t39) + m(5) * (t147 * t63 - t150 * t62) + m(4) * (t147 * t81 - t150 * t80); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t183; m(6) * (t33 * t45 + t34 * t44) + m(5) * (-t147 * t55 - t150 * t54) * t106 + t151; m(6) * (t11 * t23 + t39 * t44 + t40 * t45) + m(5) * (t41 * t18 + (-t147 * t62 - t150 * t63) * t106) + t156; m(6) * (t147 * t45 - t150 * t44); m(5) * (t106 ^ 2 * t183 + t41 ^ 2) + m(6) * (t23 ^ 2 + t44 ^ 2 + t45 ^ 2) + t156; m(6) * (t31 * t33 + t32 * t34) - t29 * t130 + ((t27 / 0.2e1 + t22 / 0.2e1) * t150 + (t21 / 0.2e1 + t26 / 0.2e1) * t147) * t129; m(6) * (t11 * t28 + t31 * t40 + t32 * t39) + t179; m(6) * (t147 * t31 - t150 * t32); m(6) * (t23 * t28 + t31 * t45 + t32 * t44) + t179; m(6) * (t28 ^ 2 + t31 ^ 2 + t32 ^ 2) + t130 ^ 2 * t29 + (t150 * t4 + t147 * t3 - t130 * (t147 * t21 + t150 * t22)) * t129;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
