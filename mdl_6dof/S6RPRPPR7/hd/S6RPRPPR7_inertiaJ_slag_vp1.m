% Calculate joint inertia matrix for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:56:14
% EndTime: 2019-03-09 02:56:18
% DurationCPUTime: 2.05s
% Computational Cost: add. (2546->275), mult. (3427->404), div. (0->0), fcn. (3403->8), ass. (0->141)
t223 = -Icges(5,4) - Icges(6,6);
t222 = Icges(5,1) + Icges(6,2);
t221 = Icges(5,2) + Icges(6,3);
t125 = qJ(3) + pkin(9);
t117 = cos(t125);
t220 = t223 * t117;
t116 = sin(t125);
t219 = t223 * t116;
t218 = t117 * t221 - t219;
t217 = t116 * t222 - t220;
t216 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t130 = sin(qJ(3));
t133 = cos(qJ(3));
t215 = Icges(4,5) * t130 + Icges(4,6) * t133 + (-Icges(6,5) + Icges(5,6)) * t117 + (-Icges(6,4) + Icges(5,5)) * t116;
t209 = t117 / 0.2e1;
t213 = Icges(5,5) * t209 - Icges(6,4) * t117 / 0.2e1 + (-Icges(5,6) / 0.2e1 + Icges(6,5) / 0.2e1) * t116 + Icges(4,5) * t133 - Icges(4,6) * t130;
t131 = sin(qJ(1));
t212 = -t131 / 0.2e1;
t201 = t131 / 0.2e1;
t134 = cos(qJ(1));
t211 = -t134 / 0.2e1;
t210 = t134 / 0.2e1;
t195 = rSges(6,2) * t116;
t136 = -t195 + (-rSges(6,3) - qJ(5)) * t117;
t185 = t116 * t134;
t107 = pkin(4) * t185;
t119 = t134 * qJ(2);
t128 = -qJ(4) - pkin(7);
t176 = pkin(3) * t130 * t134 + t128 * t131;
t166 = t119 + t176;
t163 = t107 + t166;
t22 = (-rSges(6,1) - pkin(1)) * t131 + t136 * t134 + t163;
t186 = t116 * t131;
t106 = pkin(4) * t186;
t182 = t130 * t131;
t111 = pkin(3) * t182;
t175 = pkin(1) * t134 + qJ(2) * t131;
t137 = -t128 * t134 + t111 + t175;
t135 = t106 + t137;
t23 = t134 * rSges(6,1) + t131 * t136 + t135;
t208 = m(6) * (t131 * t22 - t134 * t23);
t132 = cos(qJ(6));
t178 = t132 * t134;
t129 = sin(qJ(6));
t181 = t131 * t129;
t77 = t117 * t178 - t181;
t180 = t131 * t132;
t183 = t129 * t134;
t78 = t117 * t183 + t180;
t161 = -t78 * rSges(7,1) - t77 * rSges(7,2);
t187 = qJ(5) * t117;
t17 = (-pkin(1) - pkin(5)) * t131 + (-t187 + (rSges(7,3) + pkin(8)) * t116) * t134 + t161 + t163;
t184 = t117 * t131;
t169 = qJ(5) * t184;
t123 = t134 * pkin(5);
t79 = -t117 * t180 - t183;
t80 = -t117 * t181 + t178;
t33 = rSges(7,1) * t80 + rSges(7,2) * t79 + rSges(7,3) * t186;
t206 = pkin(8) * t186 + t123 + t33;
t18 = t135 - t169 + t206;
t207 = m(7) * (t131 * t17 - t134 * t18);
t205 = (rSges(4,1) * t130 + rSges(4,2) * t133) * t134;
t204 = t131 * t215 + t134 * t216;
t203 = t131 * t216 - t134 * t215;
t126 = t131 ^ 2;
t127 = t134 ^ 2;
t100 = rSges(4,1) * t133 - rSges(4,2) * t130;
t199 = m(4) * t100;
t198 = pkin(3) * t133;
t64 = t134 * (-t131 * pkin(7) - t176);
t197 = t134 * (t134 * t187 - t107) + t64;
t76 = t111 + (-pkin(7) - t128) * t134;
t196 = -t106 + t169 - t76;
t120 = t134 * rSges(5,3);
t179 = t131 * t133;
t112 = pkin(3) * t179;
t89 = pkin(4) * t117 + qJ(5) * t116;
t194 = t131 * t89 + t112;
t108 = t126 + t127;
t174 = m(6) / 0.2e1 + m(7) / 0.2e1;
t44 = Icges(7,3) * t117 + (Icges(7,5) * t129 + Icges(7,6) * t132) * t116;
t45 = Icges(7,6) * t117 + (Icges(7,4) * t129 + Icges(7,2) * t132) * t116;
t46 = Icges(7,5) * t117 + (Icges(7,1) * t129 + Icges(7,4) * t132) * t116;
t173 = t117 * t44 + (t129 * t46 + t132 * t45) * t116;
t172 = (m(5) + m(6) + m(7)) * t108;
t27 = Icges(7,5) * t80 + Icges(7,6) * t79 + Icges(7,3) * t186;
t29 = Icges(7,4) * t80 + Icges(7,2) * t79 + Icges(7,6) * t186;
t31 = Icges(7,1) * t80 + Icges(7,4) * t79 + Icges(7,5) * t186;
t11 = t117 * t27 + (t129 * t31 + t132 * t29) * t116;
t14 = t186 * t44 + t45 * t79 + t46 * t80;
t171 = t11 / 0.2e1 + t14 / 0.2e1;
t26 = Icges(7,5) * t78 + Icges(7,6) * t77 - Icges(7,3) * t185;
t28 = Icges(7,4) * t78 + Icges(7,2) * t77 - Icges(7,6) * t185;
t30 = Icges(7,1) * t78 + Icges(7,4) * t77 - Icges(7,5) * t185;
t10 = t117 * t26 + (t129 * t30 + t132 * t28) * t116;
t13 = -t185 * t44 + t45 * t77 + t46 * t78;
t170 = -t13 / 0.2e1 - t10 / 0.2e1;
t168 = -rSges(5,1) * t186 - rSges(5,2) * t184 - t120;
t167 = rSges(4,1) * t182 + rSges(4,2) * t179 + rSges(4,3) * t134;
t165 = -t89 - t198;
t47 = rSges(7,3) * t117 + (rSges(7,1) * t129 + rSges(7,2) * t132) * t116;
t164 = pkin(8) * t117 + t47;
t159 = rSges(5,1) * t116 + rSges(5,2) * t117;
t158 = rSges(6,3) * t117 + t195;
t20 = t117 * t33 - t186 * t47;
t32 = -rSges(7,3) * t185 - t161;
t21 = -t117 * t32 - t185 * t47;
t150 = t131 * t21 - t134 * t20;
t24 = t131 * t164 + t194;
t25 = (-t164 + t165) * t134;
t148 = t131 * t24 - t134 * t25;
t90 = -rSges(6,2) * t117 + rSges(6,3) * t116;
t37 = t131 * t90 + t194;
t38 = (t165 - t90) * t134;
t147 = t131 * t37 - t134 * t38;
t101 = rSges(2,1) * t134 - rSges(2,2) * t131;
t99 = -rSges(2,1) * t131 - rSges(2,2) * t134;
t91 = rSges(5,1) * t117 - rSges(5,2) * t116;
t75 = -rSges(3,2) * t134 + rSges(3,3) * t131 + t175;
t74 = rSges(3,3) * t134 + t119 + (rSges(3,2) - pkin(1)) * t131;
t49 = (-t91 - t198) * t134;
t48 = t131 * t91 + t112;
t42 = pkin(7) * t134 + t167 + t175;
t41 = t119 + t205 + (-rSges(4,3) - pkin(1) - pkin(7)) * t131;
t36 = t137 - t168;
t35 = t159 * t134 + (-rSges(5,3) - pkin(1)) * t131 + t166;
t34 = -t131 * t167 + (t131 * rSges(4,3) - t205) * t134;
t19 = t64 - t159 * t127 + (t168 - t76 + t120) * t131;
t16 = t173 * t117;
t15 = (t131 * t32 + t134 * t33) * t116;
t12 = t158 * t127 + (t131 * t158 + t196) * t131 + t197;
t9 = t186 * t27 + t29 * t79 + t31 * t80;
t8 = t186 * t26 + t28 * t79 + t30 * t80;
t7 = -t185 * t27 + t29 * t77 + t31 * t78;
t6 = -t185 * t26 + t28 * t77 + t30 * t78;
t5 = (-pkin(8) * t185 + t32) * t134 + (t196 + t123 - t206) * t131 + t197;
t4 = t131 * t8 + t134 * t9;
t3 = t131 * t6 + t134 * t7;
t2 = t14 * t117 + (t131 * t9 - t134 * t8) * t116;
t1 = t13 * t117 + (t131 * t7 - t134 * t6) * t116;
t39 = [Icges(4,1) * t133 ^ 2 + Icges(3,1) + Icges(2,3) + m(7) * (t17 ^ 2 + t18 ^ 2) + m(5) * (t35 ^ 2 + t36 ^ 2) + m(6) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t41 ^ 2 + t42 ^ 2) + m(3) * (t74 ^ 2 + t75 ^ 2) + m(2) * (t101 ^ 2 + t99 ^ 2) + t173 + (t117 * t222 + t219) * t117 + (t116 * t221 + t220) * t116 + (-0.2e1 * Icges(4,4) * t133 + Icges(4,2) * t130) * t130; t207 + m(5) * (t131 * t35 - t134 * t36) + t208 + m(4) * (t131 * t41 - t134 * t42) + m(3) * (t131 * t74 - t134 * t75); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t108 + t172; m(7) * (t17 * t24 + t18 * t25) + m(5) * (t35 * t48 + t36 * t49) + m(6) * (t22 * t37 + t23 * t38) + (-t42 * t199 + (Icges(6,4) * t211 + Icges(5,5) * t210 + t201 * t217) * t117 + t171 + (Icges(6,5) * t210 + Icges(5,6) * t211 + t212 * t218) * t116 + t213 * t134) * t134 + (t41 * t199 + (Icges(6,4) * t212 + Icges(5,5) * t201 + t211 * t217) * t117 - t170 + (Icges(6,5) * t201 + Icges(5,6) * t212 + t210 * t218) * t116 + t213 * t131) * t131; m(5) * (t131 * t48 - t134 * t49) + m(6) * t147 + m(7) * t148 + t108 * t199; m(7) * (t24 ^ 2 + t25 ^ 2 + t5 ^ 2) + m(6) * (t12 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(5) * (t19 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(4) * (t100 ^ 2 * t108 + t34 ^ 2) + (t127 * t204 + t4) * t134 + (t3 + t203 * t126 + (t131 * t204 + t134 * t203) * t134) * t131; m(7) * (t131 * t18 + t134 * t17) + m(5) * (t131 * t36 + t134 * t35) + m(6) * (t131 * t23 + t134 * t22); 0; m(7) * (t131 * t25 + t134 * t24) + m(6) * (t131 * t38 + t134 * t37) + m(5) * (t131 * t49 + t134 * t48); t172; 0.2e1 * (-t207 / 0.2e1 - t208 / 0.2e1) * t117; -0.2e1 * t174 * t108 * t117; m(7) * (t116 * t5 - t117 * t148) + m(6) * (t116 * t12 - t117 * t147); 0; 0.2e1 * t174 * (t108 * t117 ^ 2 + t116 ^ 2); m(7) * (t17 * t21 + t18 * t20) + t16 + (t131 * t171 + t134 * t170) * t116; m(7) * t150; m(7) * (t15 * t5 + t20 * t25 + t21 * t24) + t1 * t201 + t2 * t210 + (t10 * t131 + t11 * t134) * t209 + (t4 * t201 + t211 * t3) * t116; m(7) * (t131 * t20 + t134 * t21); m(7) * (t116 * t15 - t117 * t150); t117 * t16 + m(7) * (t15 ^ 2 + t20 ^ 2 + t21 ^ 2) + (t131 * t2 - t134 * t1 + t117 * (-t10 * t134 + t11 * t131)) * t116;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t39(1) t39(2) t39(4) t39(7) t39(11) t39(16); t39(2) t39(3) t39(5) t39(8) t39(12) t39(17); t39(4) t39(5) t39(6) t39(9) t39(13) t39(18); t39(7) t39(8) t39(9) t39(10) t39(14) t39(19); t39(11) t39(12) t39(13) t39(14) t39(15) t39(20); t39(16) t39(17) t39(18) t39(19) t39(20) t39(21);];
Mq  = res;
