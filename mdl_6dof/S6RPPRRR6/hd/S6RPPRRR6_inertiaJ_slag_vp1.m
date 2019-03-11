% Calculate joint inertia matrix for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:30:33
% EndTime: 2019-03-09 02:30:37
% DurationCPUTime: 2.15s
% Computational Cost: add. (4040->344), mult. (6441->495), div. (0->0), fcn. (6870->8), ass. (0->162)
t152 = cos(qJ(4));
t202 = Icges(5,5) * t152;
t201 = t202 / 0.2e1;
t150 = sin(qJ(1));
t145 = t150 ^ 2;
t153 = cos(qJ(1));
t146 = t153 ^ 2;
t149 = sin(qJ(4));
t200 = t149 / 0.2e1;
t199 = -t150 / 0.2e1;
t198 = t153 / 0.2e1;
t197 = -rSges(5,3) - pkin(7);
t129 = rSges(5,1) * t152 - rSges(5,2) * t149;
t196 = m(5) * t129;
t133 = t145 + t146;
t121 = m(5) * t133;
t195 = -pkin(1) - qJ(3);
t151 = cos(qJ(5));
t138 = pkin(5) * t151 + pkin(4);
t194 = -pkin(4) + t138;
t180 = t150 * t152;
t136 = pkin(8) * t180;
t154 = -pkin(9) - pkin(8);
t177 = t152 * t154;
t148 = sin(qJ(5));
t184 = t148 * t153;
t147 = qJ(5) + qJ(6);
t139 = sin(t147);
t140 = cos(t147);
t183 = t149 * t150;
t103 = -t139 * t183 + t140 * t153;
t104 = t139 * t153 + t140 * t183;
t161 = -t104 * rSges(7,1) - t103 * rSges(7,2);
t62 = -rSges(7,3) * t180 - t161;
t193 = -t62 - pkin(5) * t184 - t136 - (t149 * t194 + t177) * t150;
t84 = Icges(7,3) * t149 + (Icges(7,5) * t140 - Icges(7,6) * t139) * t152;
t86 = Icges(7,5) * t149 + (Icges(7,1) * t140 - Icges(7,4) * t139) * t152;
t192 = t140 * t152 * t86 + t149 * t84;
t178 = t152 * t153;
t182 = t149 * t153;
t105 = -t139 * t182 - t140 * t150;
t106 = -t139 * t150 + t140 * t182;
t63 = rSges(7,1) * t106 + rSges(7,2) * t105 - rSges(7,3) * t178;
t87 = rSges(7,3) * t149 + (rSges(7,1) * t140 - rSges(7,2) * t139) * t152;
t46 = t149 * t63 + t178 * t87;
t93 = Icges(6,3) * t149 + (Icges(6,5) * t151 - Icges(6,6) * t148) * t152;
t99 = Icges(6,5) * t149 + (Icges(6,1) * t151 - Icges(6,4) * t148) * t152;
t191 = t151 * t152 * t99 + t149 * t93;
t83 = t194 * t152 + (-pkin(8) - t154) * t149;
t190 = t83 + t87;
t85 = Icges(7,6) * t149 + (Icges(7,4) * t140 - Icges(7,2) * t139) * t152;
t189 = t139 * t85;
t96 = Icges(6,6) * t149 + (Icges(6,4) * t151 - Icges(6,2) * t148) * t152;
t188 = t148 * t96;
t187 = Icges(5,4) * t149;
t185 = t148 * t150;
t181 = t150 * t151;
t179 = t151 * t153;
t113 = -t148 * t182 - t181;
t114 = t149 * t179 - t185;
t176 = rSges(6,1) * t114 + rSges(6,2) * t113;
t175 = t138 * t182 + t153 * t177;
t174 = rSges(5,1) * t182 + rSges(5,2) * t178;
t173 = pkin(1) * t153 + qJ(2) * t150;
t66 = Icges(6,5) * t114 + Icges(6,6) * t113 - Icges(6,3) * t178;
t68 = Icges(6,4) * t114 + Icges(6,2) * t113 - Icges(6,6) * t178;
t70 = Icges(6,1) * t114 + Icges(6,4) * t113 - Icges(6,5) * t178;
t32 = t149 * t66 + (-t148 * t68 + t151 * t70) * t152;
t40 = t113 * t96 + t114 * t99 - t178 * t93;
t172 = -t32 / 0.2e1 - t40 / 0.2e1;
t111 = -t148 * t183 + t179;
t112 = t149 * t181 + t184;
t65 = Icges(6,5) * t112 + Icges(6,6) * t111 - Icges(6,3) * t180;
t67 = Icges(6,4) * t112 + Icges(6,2) * t111 - Icges(6,6) * t180;
t69 = Icges(6,1) * t112 + Icges(6,4) * t111 - Icges(6,5) * t180;
t31 = t149 * t65 + (-t148 * t67 + t151 * t69) * t152;
t39 = t111 * t96 + t112 * t99 - t180 * t93;
t171 = -t39 / 0.2e1 - t31 / 0.2e1;
t170 = qJ(3) * t153 + t173;
t169 = -pkin(5) * t148 - pkin(7);
t168 = -t180 / 0.2e1;
t167 = -t178 / 0.2e1;
t56 = Icges(7,5) * t104 + Icges(7,6) * t103 - Icges(7,3) * t180;
t58 = Icges(7,4) * t104 + Icges(7,2) * t103 - Icges(7,6) * t180;
t60 = Icges(7,1) * t104 + Icges(7,4) * t103 - Icges(7,5) * t180;
t17 = t103 * t58 + t104 * t60 - t180 * t56;
t57 = Icges(7,5) * t106 + Icges(7,6) * t105 - Icges(7,3) * t178;
t59 = Icges(7,4) * t106 + Icges(7,2) * t105 - Icges(7,6) * t178;
t61 = Icges(7,1) * t106 + Icges(7,4) * t105 - Icges(7,5) * t178;
t18 = t103 * t59 + t104 * t61 - t180 * t57;
t10 = -t150 * t18 + t153 * t17;
t19 = t105 * t58 + t106 * t60 - t178 * t56;
t20 = t105 * t59 + t106 * t61 - t178 * t57;
t11 = -t150 * t20 + t153 * t19;
t23 = t149 * t56 + (-t139 * t58 + t140 * t60) * t152;
t24 = t149 * t57 + (-t139 * t59 + t140 * t61) * t152;
t36 = t103 * t85 + t104 * t86 - t180 * t84;
t3 = t149 * t36 + (-t150 * t17 - t153 * t18) * t152;
t37 = t105 * t85 + t106 * t86 - t178 * t84;
t4 = t149 * t37 + (-t150 * t19 - t153 * t20) * t152;
t166 = t10 * t168 + t11 * t167 + t3 * t198 + t4 * t199 + (-t24 * t150 + t23 * t153) * t200;
t165 = t121 + (m(4) + m(6) + m(7)) * t133;
t44 = (-t152 * t189 + t192) * t149;
t164 = t44 + (t23 + t36) * t168 + (t24 + t37) * t167;
t137 = pkin(4) * t182;
t116 = -pkin(8) * t178 + t137;
t163 = -rSges(5,1) * t149 - rSges(5,2) * t152;
t162 = -rSges(6,1) * t112 - rSges(6,2) * t111;
t157 = Icges(5,2) * t152 + t187;
t156 = Icges(5,5) * t149 + Icges(5,6) * t152;
t5 = t149 * (t44 + (-t150 * t23 - t153 * t24) * t152);
t155 = t5 + (-t150 * t3 - t153 * t4) * t152;
t143 = t153 * qJ(2);
t131 = t152 * pkin(4) + t149 * pkin(8);
t130 = rSges(2,1) * t153 - rSges(2,2) * t150;
t128 = -rSges(2,1) * t150 - rSges(2,2) * t153;
t125 = -Icges(5,6) * t149 + t202;
t118 = t153 * t131;
t117 = t150 * t131;
t115 = pkin(4) * t183 - t136;
t108 = -rSges(3,2) * t153 + rSges(3,3) * t150 + t173;
t107 = rSges(3,3) * t153 + t143 + (rSges(3,2) - pkin(1)) * t150;
t102 = rSges(6,3) * t149 + (rSges(6,1) * t151 - rSges(6,2) * t148) * t152;
t95 = -Icges(5,3) * t150 + t153 * t156;
t94 = Icges(5,3) * t153 + t150 * t156;
t90 = rSges(4,2) * t150 + rSges(4,3) * t153 + t170;
t89 = rSges(4,2) * t153 + t143 + (-rSges(4,3) + t195) * t150;
t78 = t102 * t153 + t118;
t77 = t102 * t150 + t117;
t76 = t150 * t197 + t170 + t174;
t75 = t143 + t197 * t153 + (t163 + t195) * t150;
t74 = -pkin(5) * t185 - t116 + t175;
t72 = -rSges(6,3) * t178 + t176;
t71 = -rSges(6,3) * t180 - t162;
t64 = t145 * t163 - t153 * t174;
t54 = t63 * t180;
t53 = t153 * t190 + t118;
t52 = t150 * t190 + t117;
t51 = t102 * t178 + t149 * t72;
t50 = -t102 * t180 - t149 * t71;
t49 = -pkin(7) * t150 + t137 + (-rSges(6,3) - pkin(8)) * t178 + t170 + t176;
t48 = -t153 * pkin(7) + t136 + t143 + (rSges(6,3) * t152 - pkin(4) * t149 + t195) * t150 + t162;
t47 = (-t152 * t188 + t191) * t149;
t45 = -t149 * t62 - t180 * t87;
t43 = t150 * t169 + t170 + t175 + t63;
t42 = t143 + t169 * t153 + (-t138 * t149 + (rSges(7,3) - t154) * t152 + t195) * t150 + t161;
t41 = (t150 * t72 - t153 * t71) * t152;
t38 = -t178 * t62 + t54;
t33 = (-t116 - t72) * t153 + (-t115 - t71) * t150;
t30 = t149 * t74 + t178 * t83 + t46;
t29 = t149 * t193 - t180 * t190;
t28 = t113 * t68 + t114 * t70 - t178 * t66;
t27 = t113 * t67 + t114 * t69 - t178 * t65;
t26 = t111 * t68 + t112 * t70 - t180 * t66;
t25 = t111 * t67 + t112 * t69 - t180 * t65;
t16 = t54 + (t150 * t74 + t153 * t193) * t152;
t15 = (-t116 - t63 - t74) * t153 + (-t115 + t193) * t150;
t14 = -t150 * t28 + t153 * t27;
t13 = -t150 * t26 + t153 * t25;
t7 = t149 * t40 + (-t150 * t27 - t153 * t28) * t152;
t6 = t149 * t39 + (-t150 * t25 - t153 * t26) * t152;
t1 = [-t149 * (Icges(5,4) * t152 - Icges(5,2) * t149) + Icges(3,1) + Icges(4,1) + Icges(2,3) + (Icges(5,1) * t152 - t187 - t188 - t189) * t152 + m(7) * (t42 ^ 2 + t43 ^ 2) + m(6) * (t48 ^ 2 + t49 ^ 2) + m(5) * (t75 ^ 2 + t76 ^ 2) + m(3) * (t107 ^ 2 + t108 ^ 2) + m(4) * (t89 ^ 2 + t90 ^ 2) + m(2) * (t128 ^ 2 + t130 ^ 2) + t191 + t192; m(7) * (t150 * t42 - t153 * t43) + m(6) * (t150 * t48 - t153 * t49) + m(5) * (t150 * t75 - t153 * t76) + m(3) * (t107 * t150 - t108 * t153) + m(4) * (t150 * t89 - t153 * t90); m(3) * t133 + t165; m(7) * (t150 * t43 + t153 * t42) + m(6) * (t150 * t49 + t153 * t48) + m(5) * (t150 * t76 + t153 * t75) + m(4) * (t150 * t90 + t153 * t89); 0; t165; m(7) * (t42 * t53 + t43 * t52) + m(6) * (t48 * t78 + t49 * t77) + (t23 / 0.2e1 + t36 / 0.2e1 + t75 * t196 + t153 * t201 - t149 * (Icges(5,6) * t153 + t150 * t157) / 0.2e1 + t125 * t198 - t171) * t153 + (-t37 / 0.2e1 - t24 / 0.2e1 + t76 * t196 + t157 * t153 * t200 + t172 + (t201 - Icges(5,6) * t200 + t125 / 0.2e1) * t150) * t150; m(6) * (t150 * t78 - t153 * t77) + m(7) * (t150 * t53 - t153 * t52); m(6) * (t150 * t77 + t153 * t78) + m(7) * (t150 * t52 + t153 * t53) + t129 * t121; m(7) * (t15 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(6) * (t33 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(5) * (t129 ^ 2 * t133 + t64 ^ 2) + (-t145 * t95 - t11 - t14) * t150 + (t146 * t94 + t10 + t13 + (t150 * t94 - t153 * t95) * t150) * t153; t47 + m(7) * (t29 * t42 + t30 * t43) + m(6) * (t48 * t50 + t49 * t51) + (t150 * t171 + t153 * t172) * t152 + t164; m(6) * (t150 * t50 - t153 * t51) + m(7) * (t150 * t29 - t153 * t30); m(6) * (t150 * t51 + t153 * t50) + m(7) * (t150 * t30 + t153 * t29); (-t32 * t150 + t31 * t153) * t200 + t7 * t199 + t6 * t198 + (t13 * t199 - t153 * t14 / 0.2e1) * t152 + m(7) * (t15 * t16 + t29 * t53 + t30 * t52) + m(6) * (t33 * t41 + t50 * t78 + t51 * t77) + t166; t149 * t47 + t5 + m(7) * (t16 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(6) * (t41 ^ 2 + t50 ^ 2 + t51 ^ 2) + ((-t149 * t32 - t4 - t7) * t153 + (-t149 * t31 - t3 - t6) * t150) * t152; m(7) * (t42 * t45 + t43 * t46) + t164; m(7) * (t150 * t45 - t153 * t46); m(7) * (t150 * t46 + t153 * t45); m(7) * (t15 * t38 + t45 * t53 + t46 * t52) + t166; m(7) * (t16 * t38 + t29 * t45 + t30 * t46) + t155; m(7) * (t38 ^ 2 + t45 ^ 2 + t46 ^ 2) + t155;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
