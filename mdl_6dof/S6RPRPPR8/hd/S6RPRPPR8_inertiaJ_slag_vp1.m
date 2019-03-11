% Calculate joint inertia matrix for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:47
% EndTime: 2019-03-09 02:58:52
% DurationCPUTime: 2.40s
% Computational Cost: add. (1522->275), mult. (3535->398), div. (0->0), fcn. (3503->6), ass. (0->141)
t217 = -Icges(4,4) - Icges(6,4);
t216 = Icges(6,1) + Icges(4,2);
t211 = Icges(4,6) + Icges(6,5);
t123 = sin(qJ(3));
t215 = t217 * t123;
t214 = Icges(4,1) + Icges(5,1) + Icges(6,2);
t210 = Icges(4,5) + Icges(5,4) + Icges(6,6);
t126 = cos(qJ(3));
t213 = t216 * t126 - t215;
t212 = (Icges(5,5) + t217) * t126;
t209 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t208 = (-Icges(5,6) + t211) * t126 + t210 * t123;
t207 = t210 * t126;
t206 = (t214 * t123 - t212) * t126;
t124 = sin(qJ(1));
t190 = t124 / 0.2e1;
t127 = cos(qJ(1));
t204 = -t127 / 0.2e1;
t203 = t127 / 0.2e1;
t201 = t126 / 0.2e1;
t159 = (-rSges(5,3) - qJ(4)) * t126;
t172 = t123 * t127;
t111 = pkin(3) * t172;
t114 = t127 * qJ(2);
t167 = t111 + t114;
t185 = rSges(5,1) * t123;
t192 = -pkin(1) - pkin(7);
t26 = (t159 + t185) * t127 + (-rSges(5,2) + t192) * t124 + t167;
t173 = t123 * t124;
t108 = pkin(3) * t173;
t166 = t127 * pkin(1) + t124 * qJ(2);
t162 = t127 * pkin(7) + t166;
t157 = t108 + t162;
t116 = t127 * rSges(5,2);
t169 = rSges(5,1) * t173 + t116;
t27 = t124 * t159 + t157 + t169;
t200 = m(5) * (t124 * t26 - t127 * t27);
t184 = rSges(6,2) * t123;
t128 = -t184 + (-rSges(6,1) - qJ(4)) * t126;
t168 = -pkin(4) * t172 - t124 * qJ(5);
t156 = t167 - t168;
t22 = t128 * t127 + (rSges(6,3) + t192) * t124 + t156;
t107 = pkin(4) * t173;
t138 = t107 + t157;
t23 = (-rSges(6,3) - qJ(5)) * t127 + t128 * t124 + t138;
t199 = m(6) * (t124 * t22 - t127 * t23);
t122 = sin(qJ(6));
t125 = cos(qJ(6));
t170 = t126 * t127;
t76 = -t122 * t170 - t124 * t125;
t77 = -t124 * t122 + t125 * t170;
t154 = -t77 * rSges(7,1) - t76 * rSges(7,2);
t160 = (-pkin(5) - qJ(4)) * t126;
t16 = t192 * t124 + (t160 + (rSges(7,3) + pkin(8)) * t123) * t127 + t154 + t156;
t174 = qJ(5) * t127;
t171 = t124 * t126;
t74 = t122 * t171 - t125 * t127;
t75 = -t122 * t127 - t125 * t171;
t35 = t75 * rSges(7,1) + t74 * rSges(7,2) + rSges(7,3) * t173;
t197 = -pkin(8) * t173 - t35;
t17 = t124 * t160 + t138 - t174 - t197;
t198 = m(7) * (t124 * t16 - t127 * t17);
t196 = (rSges(4,1) * t123 + rSges(4,2) * t126) * t127;
t195 = t124 * t208 + t209 * t127;
t194 = t124 * t209 - t208 * t127;
t119 = t124 ^ 2;
t121 = t127 ^ 2;
t193 = m(5) / 0.2e1;
t98 = rSges(4,1) * t126 - rSges(4,2) * t123;
t191 = m(4) * t98;
t45 = Icges(7,3) * t126 + (Icges(7,5) * t125 - Icges(7,6) * t122) * t123;
t59 = Icges(7,5) * t126 + (Icges(7,1) * t125 - Icges(7,4) * t122) * t123;
t189 = t123 * t125 * t59 + t126 * t45;
t69 = t127 * (qJ(4) * t170 - t111);
t188 = t127 * t168 + t69;
t78 = -qJ(4) * t171 + t108;
t187 = -t78 - t107 + t174;
t101 = t119 + t121;
t186 = (m(6) + m(7)) * t101;
t52 = Icges(7,6) * t126 + (Icges(7,4) * t125 - Icges(7,2) * t122) * t123;
t183 = t122 * t52;
t66 = rSges(7,3) * t126 + (rSges(7,1) * t125 - rSges(7,2) * t122) * t123;
t182 = pkin(5) * t123 + pkin(8) * t126 + t66;
t96 = pkin(3) * t126 + qJ(4) * t123;
t79 = t124 * t96;
t181 = pkin(4) * t171 + t79;
t176 = Icges(5,5) * t123;
t29 = Icges(7,5) * t75 + Icges(7,6) * t74 + Icges(7,3) * t173;
t31 = Icges(7,4) * t75 + Icges(7,2) * t74 + Icges(7,6) * t173;
t33 = Icges(7,1) * t75 + Icges(7,4) * t74 + Icges(7,5) * t173;
t10 = t126 * t29 + (-t122 * t31 + t125 * t33) * t123;
t13 = t45 * t173 + t52 * t74 + t59 * t75;
t165 = t10 / 0.2e1 + t13 / 0.2e1;
t30 = Icges(7,5) * t77 + Icges(7,6) * t76 - Icges(7,3) * t172;
t32 = Icges(7,4) * t77 + Icges(7,2) * t76 - Icges(7,6) * t172;
t34 = Icges(7,1) * t77 + Icges(7,4) * t76 - Icges(7,5) * t172;
t11 = t126 * t30 + (-t122 * t32 + t125 * t34) * t123;
t14 = -t45 * t172 + t76 * t52 + t77 * t59;
t164 = -t11 / 0.2e1 - t14 / 0.2e1;
t163 = rSges(4,1) * t173 + rSges(4,2) * t171 + t127 * rSges(4,3);
t161 = -pkin(4) * t126 - t96;
t158 = t193 + m(6) / 0.2e1 + m(7) / 0.2e1;
t155 = t210 * t201 + (Icges(5,6) / 0.2e1 - t211 / 0.2e1) * t123;
t152 = rSges(6,1) * t126 + t184;
t20 = t126 * t35 - t66 * t173;
t36 = -rSges(7,3) * t172 - t154;
t21 = -t126 * t36 - t66 * t172;
t144 = t21 * t124 - t127 * t20;
t24 = t182 * t124 + t181;
t25 = (t161 - t182) * t127;
t142 = t24 * t124 - t127 * t25;
t94 = rSges(6,1) * t123 - rSges(6,2) * t126;
t37 = t124 * t94 + t181;
t38 = (t161 - t94) * t127;
t140 = t37 * t124 - t127 * t38;
t97 = rSges(5,1) * t126 + rSges(5,3) * t123;
t41 = t124 * t97 + t79;
t42 = (-t96 - t97) * t127;
t139 = t41 * t124 - t127 * t42;
t130 = -Icges(5,3) * t126 + t176;
t99 = rSges(2,1) * t127 - t124 * rSges(2,2);
t95 = -t124 * rSges(2,1) - rSges(2,2) * t127;
t68 = -rSges(3,2) * t127 + t124 * rSges(3,3) + t166;
t67 = rSges(3,3) * t127 + t114 + (rSges(3,2) - pkin(1)) * t124;
t40 = t162 + t163;
t39 = t114 + t196 + (-rSges(4,3) + t192) * t124;
t28 = -t124 * t163 + (t124 * rSges(4,3) - t196) * t127;
t19 = (-t123 * t183 + t189) * t126;
t18 = t69 + (rSges(5,3) * t126 - t185) * t121 + (rSges(5,3) * t171 + t116 - t169 - t78) * t124;
t15 = (t124 * t36 + t127 * t35) * t123;
t12 = t152 * t121 + (t152 * t124 + t187) * t124 + t188;
t9 = -t30 * t172 + t76 * t32 + t77 * t34;
t8 = -t29 * t172 + t76 * t31 + t77 * t33;
t7 = t30 * t173 + t32 * t74 + t34 * t75;
t6 = t29 * t173 + t31 * t74 + t33 * t75;
t5 = (t36 + (pkin(5) * t126 - pkin(8) * t123) * t127) * t127 + (pkin(5) * t171 + t187 + t197) * t124 + t188;
t4 = t9 * t124 + t127 * t8;
t3 = t7 * t124 + t127 * t6;
t2 = t14 * t126 + (t124 * t8 - t127 * t9) * t123;
t1 = t13 * t126 + (t124 * t6 - t127 * t7) * t123;
t43 = [Icges(3,1) + Icges(2,3) + m(7) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t39 ^ 2 + t40 ^ 2) + m(3) * (t67 ^ 2 + t68 ^ 2) + m(2) * (t95 ^ 2 + t99 ^ 2) + t189 + (t214 * t126 + t176 + t215) * t126 + (-t183 + (Icges(5,3) + t216) * t123 + t212) * t123; t198 + t200 + t199 + m(4) * (t124 * t39 - t127 * t40) + m(3) * (t124 * t67 - t127 * t68); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t193) * t101 + t186; m(7) * (t16 * t24 + t17 * t25) + m(5) * (t26 * t41 + t27 * t42) + m(6) * (t22 * t37 + t23 * t38) + (t155 * t127 + t190 * t206 - t40 * t191 + t203 * t207 + t165) * t127 + (t155 * t124 + t190 * t207 + t39 * t191 + t204 * t206 - t164) * t124 + ((Icges(5,6) * t203 + t130 * t190 + t211 * t204) * t127 + (Icges(5,6) * t190 + t130 * t204 + t213 * t203) * t124 - (t211 * t124 + t213 * t127) * t124 / 0.2e1) * t123; m(5) * t139 + m(6) * t140 + m(7) * t142 + t101 * t191; m(7) * (t24 ^ 2 + t25 ^ 2 + t5 ^ 2) + m(5) * (t18 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t12 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(4) * (t101 * t98 ^ 2 + t28 ^ 2) + (t195 * t121 + t3) * t127 + (t4 + t194 * t119 + (t195 * t124 + t194 * t127) * t127) * t124; 0.2e1 * (-t198 / 0.2e1 - t200 / 0.2e1 - t199 / 0.2e1) * t126; -0.2e1 * t158 * t101 * t126; m(7) * (t123 * t5 - t142 * t126) + m(5) * (t123 * t18 - t139 * t126) + m(6) * (t123 * t12 - t140 * t126); 0.2e1 * t158 * (t101 * t126 ^ 2 + t123 ^ 2); m(7) * (-t124 * t17 - t127 * t16) + m(6) * (-t124 * t23 - t127 * t22); 0; m(7) * (-t124 * t25 - t127 * t24) + m(6) * (-t124 * t38 - t127 * t37); 0; t186; t19 + m(7) * (t16 * t21 + t17 * t20) + (t165 * t124 + t164 * t127) * t123; m(7) * t144; m(7) * (t15 * t5 + t20 * t25 + t21 * t24) + t1 * t203 + (t10 * t127 + t11 * t124) * t201 + t2 * t190 + (t3 * t190 + t4 * t204) * t123; m(7) * (t15 * t123 - t144 * t126); m(7) * (-t20 * t124 - t127 * t21); t126 * t19 + m(7) * (t15 ^ 2 + t20 ^ 2 + t21 ^ 2) + (t124 * t1 - t127 * t2 + t126 * (t10 * t124 - t11 * t127)) * t123;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t43(1) t43(2) t43(4) t43(7) t43(11) t43(16); t43(2) t43(3) t43(5) t43(8) t43(12) t43(17); t43(4) t43(5) t43(6) t43(9) t43(13) t43(18); t43(7) t43(8) t43(9) t43(10) t43(14) t43(19); t43(11) t43(12) t43(13) t43(14) t43(15) t43(20); t43(16) t43(17) t43(18) t43(19) t43(20) t43(21);];
Mq  = res;
