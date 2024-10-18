% Calculate joint inertia matrix for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR11_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR11_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:11
% EndTime: 2024-09-27 21:45:12
% DurationCPUTime: 0.47s
% Computational Cost: add. (13825->326), mult. (8140->450), div. (0->0), fcn. (7608->20), ass. (0->172)
t184 = pkin(7) + pkin(8);
t178 = qJ(3) + qJ(4);
t171 = pkin(5) - t178;
t163 = -qJ(5) + t171;
t155 = cos(t163) / 0.2e1;
t170 = pkin(5) + t178;
t162 = qJ(5) + t170;
t159 = cos(t162);
t139 = t155 + t159 / 0.2e1;
t174 = qJ(5) + t178;
t164 = sin(t174);
t176 = pkin(10) + qJ(2);
t168 = sin(t176);
t169 = cos(t176);
t113 = -t168 * t139 - t169 * t164;
t154 = sin(t162) / 0.2e1;
t158 = sin(t163);
t138 = t154 - t158 / 0.2e1;
t165 = cos(t174);
t114 = -t168 * t138 + t169 * t165;
t137 = t154 + t158 / 0.2e1;
t140 = t155 - t159 / 0.2e1;
t180 = cos(pkin(5));
t111 = t169 * t139 - t168 * t164;
t112 = t169 * t138 + t168 * t165;
t179 = sin(pkin(5));
t206 = t169 * t179;
t61 = Icges(6,5) * t112 + Icges(6,6) * t111 - Icges(6,3) * t206;
t63 = Icges(6,4) * t112 + Icges(6,2) * t111 - Icges(6,6) * t206;
t65 = Icges(6,1) * t112 + Icges(6,4) * t111 - Icges(6,5) * t206;
t13 = t137 * t63 + t140 * t65 + t180 * t61;
t207 = t168 * t179;
t62 = Icges(6,5) * t114 + Icges(6,6) * t113 + Icges(6,3) * t207;
t64 = Icges(6,4) * t114 + Icges(6,2) * t113 + Icges(6,6) * t207;
t66 = Icges(6,1) * t114 + Icges(6,4) * t113 + Icges(6,5) * t207;
t14 = t137 * t64 + t140 * t66 + t180 * t62;
t100 = Icges(6,4) * t140 + Icges(6,2) * t137 + Icges(6,6) * t180;
t101 = Icges(6,1) * t140 + Icges(6,4) * t137 + Icges(6,5) * t180;
t99 = Icges(6,5) * t140 + Icges(6,6) * t137 + Icges(6,3) * t180;
t25 = t113 * t100 + t114 * t101 + t99 * t207;
t200 = t137 * t100 + t140 * t101 + t180 * t99;
t35 = t200 * t180;
t213 = ((t113 * t64 + t114 * t66 + t62 * t207) * t207 - (t113 * t63 + t114 * t65 + t61 * t207) * t206 + t25 * t180) * t207 + t180 * (t35 + (-t13 * t169 + t14 * t168) * t179);
t212 = t169 * pkin(2);
t183 = cos(qJ(3));
t167 = t183 * pkin(3) + pkin(2);
t187 = -t112 * rSges(6,1) - t111 * rSges(6,2);
t67 = -rSges(6,3) * t206 - t187;
t68 = t114 * rSges(6,1) + t113 * rSges(6,2) + rSges(6,3) * t207;
t32 = t68 * t206 + t67 * t207;
t56 = t180 * t68;
t173 = cos(t178);
t148 = pkin(4) * t173 + t167;
t136 = t169 * t148;
t149 = t169 * t167;
t166 = cos(qJ(4)) * pkin(4) + pkin(3);
t177 = pkin(9) + t184;
t182 = sin(qJ(3));
t199 = pkin(4) * sin(qJ(4)) * t183;
t121 = -t179 * t177 + (t166 * t182 + t199) * t180;
t204 = t180 * t182;
t145 = pkin(3) * t204 - t179 * t184;
t201 = t121 - t145;
t58 = -t201 * t168 + t136 - t149;
t211 = t180 * t58 + t56;
t57 = t201 * t169 + (t148 - t167) * t168;
t210 = -t57 - t67;
t157 = cos(t171) / 0.2e1;
t161 = cos(t170);
t143 = t157 + t161 / 0.2e1;
t172 = sin(t178);
t117 = t169 * t143 - t168 * t172;
t156 = sin(t170) / 0.2e1;
t160 = sin(t171);
t142 = t156 - t160 / 0.2e1;
t118 = t169 * t142 + t168 * t173;
t188 = -t118 * rSges(5,1) - t117 * rSges(5,2);
t77 = -rSges(5,3) * t206 - t188;
t119 = -t168 * t143 - t169 * t172;
t120 = -t168 * t142 + t169 * t173;
t78 = t120 * rSges(5,1) + t119 * rSges(5,2) + rSges(5,3) * t207;
t36 = t78 * t206 + t77 * t207;
t153 = pkin(7) * t206;
t96 = t169 * t145 + t153 + (-pkin(2) + t167) * t168;
t97 = -t212 + t149 + (-pkin(7) * t179 - t145) * t168;
t209 = t97 * t206 + t96 * t207;
t208 = t168 * t169;
t205 = t179 * t182;
t203 = t180 * t183;
t103 = t140 * rSges(6,1) + t137 * rSges(6,2) + t180 * rSges(6,3);
t202 = -(t177 - t184) * t180 - (t199 + (-pkin(3) + t166) * t182) * t179 - t103;
t141 = t156 + t160 / 0.2e1;
t144 = t157 - t161 / 0.2e1;
t105 = Icges(5,5) * t144 + Icges(5,6) * t141 + Icges(5,3) * t180;
t106 = Icges(5,4) * t144 + Icges(5,2) * t141 + Icges(5,6) * t180;
t107 = Icges(5,1) * t144 + Icges(5,4) * t141 + Icges(5,5) * t180;
t198 = t180 * t105 + t141 * t106 + t144 * t107;
t71 = Icges(5,5) * t118 + Icges(5,6) * t117 - Icges(5,3) * t206;
t73 = Icges(5,4) * t118 + Icges(5,2) * t117 - Icges(5,6) * t206;
t75 = Icges(5,1) * t118 + Icges(5,4) * t117 - Icges(5,5) * t206;
t20 = t141 * t73 + t144 * t75 + t180 * t71;
t72 = Icges(5,5) * t120 + Icges(5,6) * t119 + Icges(5,3) * t207;
t74 = Icges(5,4) * t120 + Icges(5,2) * t119 + Icges(5,6) * t207;
t76 = Icges(5,1) * t120 + Icges(5,4) * t119 + Icges(5,5) * t207;
t21 = t141 * t74 + t144 * t76 + t180 * t72;
t31 = t105 * t207 + t119 * t106 + t120 * t107;
t37 = t198 * t180;
t197 = ((t119 * t74 + t120 * t76 + t72 * t207) * t207 - (t119 * t73 + t120 * t75 + t71 * t207) * t206 + t31 * t180) * t207 + t180 * (t37 + (t168 * t21 - t169 * t20) * t179) + t213;
t127 = Icges(4,3) * t180 + (Icges(4,5) * t182 + Icges(4,6) * t183) * t179;
t128 = Icges(4,6) * t180 + (Icges(4,4) * t182 + Icges(4,2) * t183) * t179;
t129 = Icges(4,5) * t180 + (Icges(4,1) * t182 + Icges(4,4) * t183) * t179;
t196 = t179 * t183 * t128 + t180 * t127 + t129 * t205;
t133 = -t168 * t203 - t169 * t182;
t134 = -t168 * t204 + t169 * t183;
t92 = t134 * rSges(4,1) + t133 * rSges(4,2) + rSges(4,3) * t207;
t195 = t207 / 0.2e1;
t194 = -t206 / 0.2e1;
t8 = t58 * t206 + t57 * t207 + t32;
t193 = t202 * t179;
t108 = t144 * rSges(5,1) + t141 * rSges(5,2) + t180 * rSges(5,3);
t135 = pkin(3) * t205 + (-pkin(7) + t184) * t180;
t192 = (-t108 - t135) * t179;
t24 = t111 * t100 + t112 * t101 - t99 * t206;
t191 = t35 + (t14 + t25) * t195 + (t13 + t24) * t194;
t190 = (-t135 + t202) * t179;
t2 = (t111 * t64 + t112 * t66 - t62 * t206) * t207 - (t111 * t63 + t112 * t65 - t61 * t206) * t206 + t24 * t180;
t189 = -t2 * t206 + t213;
t30 = -t105 * t206 + t117 * t106 + t118 * t107;
t4 = (t117 * t74 + t118 * t76 - t72 * t206) * t207 - (t117 * t73 + t118 * t75 - t71 * t206) * t206 + t30 * t180;
t186 = (-t2 - t4) * t206 + t197;
t131 = -t168 * t182 + t169 * t203;
t132 = t168 * t183 + t169 * t204;
t91 = t132 * rSges(4,1) + t131 * rSges(4,2) - rSges(4,3) * t206;
t185 = t191 + t37 + (t21 + t31) * t195 + (t20 + t30) * t194;
t147 = t169 * rSges(3,1) - t168 * rSges(3,2);
t146 = -t168 * rSges(3,1) - t169 * rSges(3,2);
t130 = t180 * rSges(4,3) + (rSges(4,1) * t182 + rSges(4,2) * t183) * t179;
t95 = t180 * t97;
t90 = Icges(4,1) * t134 + Icges(4,4) * t133 + Icges(4,5) * t207;
t89 = Icges(4,1) * t132 + Icges(4,4) * t131 - Icges(4,5) * t206;
t88 = Icges(4,4) * t134 + Icges(4,2) * t133 + Icges(4,6) * t207;
t87 = Icges(4,4) * t132 + Icges(4,2) * t131 - Icges(4,6) * t206;
t86 = Icges(4,5) * t134 + Icges(4,6) * t133 + Icges(4,3) * t207;
t85 = Icges(4,5) * t132 + Icges(4,6) * t131 - Icges(4,3) * t206;
t80 = pkin(7) * t207 + t212 + t92;
t79 = -t168 * pkin(2) + t153 - t91;
t70 = t180 * t78;
t69 = t196 * t180;
t55 = -t130 * t206 - t180 * t91;
t54 = -t130 * t207 + t180 * t92;
t48 = -t168 * t145 + t149 + t78;
t47 = -t168 * t167 + (rSges(5,3) * t179 - t145) * t169 + t188;
t46 = t127 * t207 + t133 * t128 + t134 * t129;
t45 = -t127 * t206 + t131 * t128 + t132 * t129;
t44 = -t108 * t206 - t180 * t77;
t43 = (t168 * t91 + t169 * t92) * t179;
t42 = -t108 * t207 + t70;
t41 = -t168 * t121 + t136 + t68;
t40 = -t168 * t148 + (rSges(6,3) * t179 - t121) * t169 + t187;
t39 = -t103 * t206 - t180 * t67;
t38 = -t103 * t207 + t56;
t34 = t180 * t86 + (t182 * t90 + t183 * t88) * t179;
t33 = t180 * t85 + (t182 * t89 + t183 * t87) * t179;
t27 = (-t77 - t96) * t180 + t169 * t192;
t26 = t168 * t192 + t70 + t95;
t17 = t169 * t193 + t210 * t180;
t16 = t168 * t193 + t211;
t15 = t209 + t36;
t10 = (-t96 + t210) * t180 + t169 * t190;
t9 = t168 * t190 + t211 + t95;
t7 = t8 + t209;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; Icges(3,3) + m(6) * (t40 ^ 2 + t41 ^ 2) + m(5) * (t47 ^ 2 + t48 ^ 2) + m(4) * (t79 ^ 2 + t80 ^ 2) + m(3) * (t146 ^ 2 + t147 ^ 2) + t196 + t198 + t200; m(4) * t43 + m(5) * t15 + m(6) * t7; t69 + t185 + ((-t33 / 0.2e1 - t45 / 0.2e1) * t169 + (t46 / 0.2e1 + t34 / 0.2e1) * t168) * t179 + m(6) * (t10 * t40 + t9 * t41) + m(5) * (t26 * t48 + t27 * t47) + m(4) * (t54 * t80 + t55 * t79); t180 * t69 + (-t169 * t2 - t169 * t4 + (t168 * ((t133 * t88 + t134 * t90) * t168 - (t133 * t87 + t134 * t89) * t169) - t169 * ((t131 * t88 + t132 * t90) * t168 - (t131 * t87 + t132 * t89) * t169) + (t168 * (t168 ^ 2 * t86 - t85 * t208) - t169 * (t169 ^ 2 * t85 - t86 * t208)) * t179) * t179 + ((-t33 - t45) * t169 + (t34 + t46) * t168) * t180) * t179 + m(6) * (t10 ^ 2 + t7 ^ 2 + t9 ^ 2) + m(5) * (t15 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(4) * (t43 ^ 2 + t54 ^ 2 + t55 ^ 2) + t197; m(5) * t36 + m(6) * t8; m(6) * (t16 * t41 + t17 * t40) + m(5) * (t42 * t48 + t44 * t47) + t185; m(6) * (t17 * t10 + t16 * t9 + t8 * t7) + m(5) * (t36 * t15 + t42 * t26 + t44 * t27) + t186; m(5) * (t36 ^ 2 + t42 ^ 2 + t44 ^ 2) + m(6) * (t16 ^ 2 + t17 ^ 2 + t8 ^ 2) + t186; m(6) * t32; m(6) * (t38 * t41 + t39 * t40) + t191; m(6) * (t39 * t10 + t32 * t7 + t38 * t9) + t189; m(6) * (t38 * t16 + t39 * t17 + t32 * t8) + t189; m(6) * (t32 ^ 2 + t38 ^ 2 + t39 ^ 2) + t189;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
