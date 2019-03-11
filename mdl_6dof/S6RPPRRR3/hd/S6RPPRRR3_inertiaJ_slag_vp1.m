% Calculate joint inertia matrix for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:56
% EndTime: 2019-03-09 02:22:59
% DurationCPUTime: 2.03s
% Computational Cost: add. (6470->324), mult. (6271->469), div. (0->0), fcn. (6728->10), ass. (0->164)
t154 = cos(qJ(4));
t212 = Icges(5,5) * t154;
t151 = sin(qJ(4));
t211 = Icges(5,6) * t151;
t210 = t212 / 0.2e1 - t211 / 0.2e1;
t148 = qJ(1) + pkin(10);
t144 = cos(t148);
t209 = (rSges(5,1) * t151 + rSges(5,2) * t154) * t144;
t143 = sin(t148);
t141 = t143 ^ 2;
t142 = t144 ^ 2;
t208 = -pkin(2) - pkin(7);
t186 = t144 * t154;
t149 = qJ(5) + qJ(6);
t145 = sin(t149);
t146 = cos(t149);
t188 = t143 * t154;
t184 = t145 * t151;
t96 = -t143 * t184 + t144 * t146;
t183 = t146 * t151;
t97 = t143 * t183 + t144 * t145;
t57 = Icges(7,5) * t97 + Icges(7,6) * t96 - Icges(7,3) * t188;
t59 = Icges(7,4) * t97 + Icges(7,2) * t96 - Icges(7,6) * t188;
t61 = Icges(7,1) * t97 + Icges(7,4) * t96 - Icges(7,5) * t188;
t28 = t151 * t57 + (-t145 * t59 + t146 * t61) * t154;
t98 = t143 * t146 + t144 * t184;
t99 = t143 * t145 - t144 * t183;
t58 = Icges(7,5) * t99 + Icges(7,6) * t98 + Icges(7,3) * t186;
t60 = Icges(7,4) * t99 + Icges(7,2) * t98 + Icges(7,6) * t186;
t62 = Icges(7,1) * t99 + Icges(7,4) * t98 + Icges(7,5) * t186;
t29 = t151 * t58 + (-t145 * t60 + t146 * t62) * t154;
t101 = Icges(7,6) * t151 + (Icges(7,4) * t146 - Icges(7,2) * t145) * t154;
t185 = t145 * t101;
t100 = Icges(7,3) * t151 + (Icges(7,5) * t146 - Icges(7,6) * t145) * t154;
t102 = Icges(7,5) * t151 + (Icges(7,1) * t146 - Icges(7,4) * t145) * t154;
t196 = t154 * t146 * t102 + t151 * t100;
t49 = (-t154 * t185 + t196) * t151;
t20 = t186 * t57 + t98 * t59 + t99 * t61;
t21 = t186 * t58 + t98 * t60 + t99 * t62;
t39 = t100 * t186 + t98 * t101 + t99 * t102;
t5 = t39 * t151 + (-t143 * t20 + t144 * t21) * t154;
t207 = t5 * t186 + t151 * (t49 + (-t143 * t28 + t144 * t29) * t154);
t206 = t143 / 0.2e1;
t205 = t144 / 0.2e1;
t203 = t151 / 0.2e1;
t127 = t154 * rSges(5,1) - t151 * rSges(5,2);
t201 = m(5) * t127;
t152 = sin(qJ(1));
t200 = t152 * pkin(1);
t156 = -pkin(9) - pkin(8);
t199 = -pkin(8) - t156;
t63 = t97 * rSges(7,1) + t96 * rSges(7,2) - rSges(7,3) * t188;
t189 = t143 * t151;
t132 = pkin(4) * t189;
t115 = -pkin(8) * t188 + t132;
t153 = cos(qJ(5));
t140 = t153 * pkin(5) + pkin(4);
t150 = sin(qJ(5));
t187 = t144 * t150;
t175 = pkin(5) * t187 + t140 * t189 + t156 * t188;
t73 = -t115 + t175;
t198 = -t63 - t73;
t134 = t144 * t151 * pkin(4);
t190 = t143 * t150;
t191 = t140 * t151;
t165 = -t99 * rSges(7,1) - t98 * rSges(7,2);
t64 = rSges(7,3) * t186 - t165;
t197 = t64 + pkin(5) * t190 + t134 + (t154 * t199 - t191) * t144;
t103 = t151 * rSges(7,3) + (rSges(7,1) * t146 - rSges(7,2) * t145) * t154;
t45 = t103 * t188 + t151 * t63;
t95 = (-pkin(4) + t140) * t154 + t199 * t151;
t195 = t103 + t95;
t111 = Icges(6,3) * t151 + (Icges(6,5) * t153 - Icges(6,6) * t150) * t154;
t113 = Icges(6,5) * t151 + (Icges(6,1) * t153 - Icges(6,4) * t150) * t154;
t194 = t154 * t153 * t113 + t151 * t111;
t112 = Icges(6,6) * t151 + (Icges(6,4) * t153 - Icges(6,2) * t150) * t154;
t182 = t150 * t112;
t181 = t150 * t151;
t180 = t151 * t153;
t107 = -t143 * t181 + t144 * t153;
t108 = t143 * t180 + t187;
t179 = t108 * rSges(6,1) + t107 * rSges(6,2);
t178 = t141 + t142;
t65 = Icges(6,5) * t108 + Icges(6,6) * t107 - Icges(6,3) * t188;
t67 = Icges(6,4) * t108 + Icges(6,2) * t107 - Icges(6,6) * t188;
t69 = Icges(6,1) * t108 + Icges(6,4) * t107 - Icges(6,5) * t188;
t32 = t151 * t65 + (-t150 * t67 + t153 * t69) * t154;
t43 = t107 * t112 + t108 * t113 - t111 * t188;
t177 = -t43 / 0.2e1 - t32 / 0.2e1;
t109 = t143 * t153 + t144 * t181;
t110 = -t144 * t180 + t190;
t66 = Icges(6,5) * t110 + Icges(6,6) * t109 + Icges(6,3) * t186;
t68 = Icges(6,4) * t110 + Icges(6,2) * t109 + Icges(6,6) * t186;
t70 = Icges(6,1) * t110 + Icges(6,4) * t109 + Icges(6,5) * t186;
t33 = t151 * t66 + (-t150 * t68 + t153 * t70) * t154;
t44 = t109 * t112 + t110 * t113 + t111 * t186;
t176 = t44 / 0.2e1 + t33 / 0.2e1;
t174 = rSges(5,1) * t189 + rSges(5,2) * t188 + t144 * rSges(5,3);
t155 = cos(qJ(1));
t147 = t155 * pkin(1);
t173 = t144 * pkin(2) + t143 * qJ(3) + t147;
t172 = (-rSges(6,3) - pkin(8)) * t154;
t171 = -t188 / 0.2e1;
t170 = t186 / 0.2e1;
t169 = t144 * qJ(3) - t200;
t18 = -t188 * t57 + t96 * t59 + t97 * t61;
t19 = -t188 * t58 + t96 * t60 + t97 * t62;
t11 = t19 * t143 + t18 * t144;
t12 = t21 * t143 + t20 * t144;
t38 = -t100 * t188 + t96 * t101 + t97 * t102;
t4 = t38 * t151 + (-t143 * t18 + t144 * t19) * t154;
t168 = t11 * t171 + t12 * t170 + t4 * t205 + t5 * t206 + (t29 * t143 + t28 * t144) * t203;
t167 = t144 * pkin(7) + t173;
t166 = t49 + (t28 + t38) * t171 + (t29 + t39) * t170;
t164 = -t188 * t4 + t207;
t162 = -t110 * rSges(6,1) - t109 * rSges(6,2);
t157 = Icges(5,5) * t151 + Icges(5,6) * t154;
t131 = t154 * pkin(4) + t151 * pkin(8);
t128 = t155 * rSges(2,1) - t152 * rSges(2,2);
t126 = -t152 * rSges(2,1) - t155 * rSges(2,2);
t123 = -t211 + t212;
t118 = t143 * t131;
t117 = t144 * rSges(3,1) - t143 * rSges(3,2) + t147;
t116 = -t143 * rSges(3,1) - t144 * rSges(3,2) - t200;
t114 = t151 * rSges(6,3) + (rSges(6,1) * t153 - rSges(6,2) * t150) * t154;
t94 = t144 * (pkin(8) * t186 - t134);
t86 = Icges(5,3) * t143 - t144 * t157;
t85 = Icges(5,3) * t144 + t143 * t157;
t83 = t103 * t186;
t80 = -t144 * rSges(4,2) + t143 * rSges(4,3) + t173;
t79 = t144 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t143 + t169;
t78 = (-t114 - t131) * t144;
t77 = t143 * t114 + t118;
t76 = t167 + t174;
t75 = t209 + (-rSges(5,3) + t208) * t143 + t169;
t72 = rSges(6,3) * t186 - t162;
t71 = -rSges(6,3) * t188 + t179;
t56 = -t143 * t174 + (t143 * rSges(5,3) - t209) * t144;
t54 = (-t131 - t195) * t144;
t53 = t143 * t195 + t118;
t52 = (-t154 * t182 + t194) * t151;
t51 = t114 * t186 - t151 * t72;
t50 = t114 * t188 + t151 * t71;
t48 = t143 * t172 + t132 + t167 + t179;
t47 = t143 * t208 + t144 * t172 + t134 + t162 + t169;
t46 = -t151 * t64 + t83;
t42 = t63 + t167 + t175;
t41 = (t191 + (-rSges(7,3) + t156) * t154) * t144 + (-pkin(5) * t150 + t208) * t143 + t165 + t169;
t40 = (-t143 * t72 - t144 * t71) * t154;
t35 = (-t143 * t64 - t144 * t63) * t154;
t34 = t144 * t72 + t94 + (-t115 - t71) * t143;
t31 = -t151 * t197 + t186 * t95 + t83;
t30 = t151 * t73 + t188 * t95 + t45;
t27 = t109 * t68 + t110 * t70 + t186 * t66;
t26 = t109 * t67 + t110 * t69 + t186 * t65;
t25 = t107 * t68 + t108 * t70 - t188 * t66;
t24 = t107 * t67 + t108 * t69 - t188 * t65;
t17 = (-t143 * t197 + t144 * t198) * t154;
t16 = t94 + t197 * t144 + (-t115 + t198) * t143;
t14 = t27 * t143 + t26 * t144;
t13 = t25 * t143 + t24 * t144;
t7 = t44 * t151 + (-t143 * t26 + t144 * t27) * t154;
t6 = t43 * t151 + (-t143 * t24 + t144 * t25) * t154;
t1 = [Icges(4,1) + Icges(2,3) + Icges(3,3) + (Icges(5,1) * t154 - t182 - t185) * t154 + m(7) * (t41 ^ 2 + t42 ^ 2) + m(6) * (t47 ^ 2 + t48 ^ 2) + m(5) * (t75 ^ 2 + t76 ^ 2) + m(4) * (t79 ^ 2 + t80 ^ 2) + m(3) * (t116 ^ 2 + t117 ^ 2) + m(2) * (t126 ^ 2 + t128 ^ 2) + t194 + t196 + (-0.2e1 * Icges(5,4) * t154 + Icges(5,2) * t151) * t151; 0; m(3) + m(4) + m(5) + m(6) + m(7); m(7) * (t143 * t41 - t144 * t42) + m(6) * (t143 * t47 - t144 * t48) + m(5) * (t143 * t75 - t144 * t76) + m(4) * (t143 * t79 - t144 * t80); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t178; m(7) * (t53 * t41 + t54 * t42) + m(6) * (t77 * t47 + t78 * t48) + (t28 / 0.2e1 + t38 / 0.2e1 - t76 * t201 + t123 * t205 - t177 + t210 * t144) * t144 + (t29 / 0.2e1 + t39 / 0.2e1 + t75 * t201 + t123 * t206 + t176 + t210 * t143) * t143; m(5) * t56 + m(6) * t34 + m(7) * t16; m(6) * (t77 * t143 - t78 * t144) + m(7) * (t53 * t143 - t54 * t144) + t178 * t201; m(7) * (t16 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(6) * (t34 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(5) * (t127 ^ 2 * t178 + t56 ^ 2) + (t141 * t86 + t12 + t14) * t143 + (t142 * t85 + t11 + t13 + (t143 * t85 + t144 * t86) * t143) * t144; t52 + m(7) * (t30 * t42 + t31 * t41) + m(6) * (t51 * t47 + t50 * t48) + (t143 * t177 + t144 * t176) * t154 + t166; m(6) * t40 + m(7) * t17; m(6) * (t51 * t143 - t50 * t144) + m(7) * (t31 * t143 - t30 * t144); t6 * t205 + t7 * t206 + (t33 * t143 + t32 * t144) * t203 + (t14 * t205 - t143 * t13 / 0.2e1) * t154 + m(7) * (t17 * t16 + t30 * t54 + t31 * t53) + m(6) * (t40 * t34 + t50 * t78 + t51 * t77) + t168; t151 * t52 + m(7) * (t17 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t40 ^ 2 + t50 ^ 2 + t51 ^ 2) + ((t151 * t33 + t7) * t144 + (-t151 * t32 - t4 - t6) * t143) * t154 + t207; m(7) * (t46 * t41 + t45 * t42) + t166; m(7) * t35; m(7) * (t46 * t143 - t45 * t144); m(7) * (t35 * t16 + t45 * t54 + t46 * t53) + t168; m(7) * (t35 * t17 + t45 * t30 + t46 * t31) + t164; m(7) * (t35 ^ 2 + t45 ^ 2 + t46 ^ 2) + t164;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
