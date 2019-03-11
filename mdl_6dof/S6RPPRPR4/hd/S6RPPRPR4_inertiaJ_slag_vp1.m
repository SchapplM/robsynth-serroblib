% Calculate joint inertia matrix for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:46:16
% EndTime: 2019-03-09 01:46:19
% DurationCPUTime: 1.62s
% Computational Cost: add. (3195->248), mult. (5097->367), div. (0->0), fcn. (6258->10), ass. (0->124)
t109 = qJ(4) + pkin(10);
t104 = cos(t109);
t113 = sin(qJ(1));
t116 = cos(qJ(1));
t148 = sin(pkin(9));
t149 = cos(pkin(9));
t88 = -t113 * t148 - t116 * t149;
t153 = t104 * t88;
t183 = Icges(6,6) * t104;
t115 = cos(qJ(4));
t179 = Icges(5,6) * t115;
t103 = sin(t109);
t180 = Icges(6,5) * t103;
t112 = sin(qJ(4));
t181 = Icges(5,5) * t112;
t182 = -t179 / 0.2e1 - t180 / 0.2e1 - t181 / 0.2e1;
t178 = Icges(5,3) + Icges(6,3);
t177 = -Icges(5,5) * t115 - Icges(6,5) * t104 + Icges(5,6) * t112 + Icges(6,6) * t103;
t137 = m(6) / 0.2e1 + m(7) / 0.2e1;
t173 = 0.2e1 * t137;
t89 = -t113 * t149 + t116 * t148;
t172 = t177 * t88 + t178 * t89;
t171 = -t177 * t89 + t178 * t88;
t154 = t103 * t88;
t114 = cos(qJ(6));
t111 = sin(qJ(6));
t140 = t104 * t111;
t58 = t114 * t89 + t88 * t140;
t139 = t104 * t114;
t59 = t111 * t89 - t88 * t139;
t31 = t59 * rSges(7,1) + t58 * rSges(7,2) - rSges(7,3) * t154;
t170 = -pkin(5) * t153 - pkin(8) * t154 + t31;
t86 = t88 ^ 2;
t85 = t89 ^ 2;
t169 = -t88 / 0.2e1;
t95 = -rSges(5,1) * t112 - rSges(5,2) * t115;
t168 = m(5) * t95;
t167 = pkin(7) * t89;
t166 = t104 / 0.2e1;
t165 = pkin(4) * t112;
t164 = pkin(5) * t104;
t163 = t88 * rSges(5,3);
t64 = Icges(7,3) * t104 + (-Icges(7,5) * t114 + Icges(7,6) * t111) * t103;
t65 = Icges(7,6) * t104 + (-Icges(7,4) * t114 + Icges(7,2) * t111) * t103;
t162 = t103 * t111 * t65 + t104 * t64;
t67 = rSges(7,3) * t104 + (-rSges(7,1) * t114 + rSges(7,2) * t111) * t103;
t161 = pkin(5) * t103 - pkin(8) * t104 - t67;
t102 = pkin(4) * t115 + pkin(3);
t110 = -qJ(5) - pkin(7);
t159 = -t88 * t102 - t89 * t110;
t156 = rSges(5,1) * t115;
t158 = t89 * rSges(5,3) - t88 * t156;
t157 = t85 + t86;
t155 = rSges(5,2) * t112;
t152 = t110 * t88;
t66 = Icges(7,5) * t104 + (-Icges(7,1) * t114 + Icges(7,4) * t111) * t103;
t151 = t114 * t66;
t150 = t89 * t103;
t144 = Icges(6,4) * t104;
t143 = Icges(7,5) * t103;
t142 = Icges(7,6) * t103;
t141 = Icges(7,3) * t103;
t138 = t116 * pkin(1) + t113 * qJ(2);
t56 = -t114 * t88 + t89 * t140;
t57 = -t111 * t88 - t89 * t139;
t24 = Icges(7,5) * t57 + Icges(7,6) * t56 - t89 * t141;
t26 = Icges(7,4) * t57 + Icges(7,2) * t56 - t89 * t142;
t28 = Icges(7,1) * t57 + Icges(7,4) * t56 - t89 * t143;
t10 = t104 * t24 + (t111 * t26 - t114 * t28) * t103;
t16 = -t64 * t150 + t56 * t65 + t57 * t66;
t136 = -t10 / 0.2e1 - t16 / 0.2e1;
t25 = Icges(7,5) * t59 + Icges(7,6) * t58 - t88 * t141;
t27 = Icges(7,4) * t59 + Icges(7,2) * t58 - t88 * t142;
t29 = Icges(7,1) * t59 + Icges(7,4) * t58 - t88 * t143;
t11 = t104 * t25 + (t111 * t27 - t114 * t29) * t103;
t17 = -t64 * t154 + t58 * t65 + t59 * t66;
t135 = -t17 / 0.2e1 - t11 / 0.2e1;
t134 = t116 * pkin(2) + t138;
t106 = t116 * qJ(2);
t132 = t106 + (-pkin(1) - pkin(2)) * t113;
t131 = -rSges(7,1) * t57 - rSges(7,2) * t56;
t130 = t134 + t159;
t129 = t155 - t156;
t128 = -rSges(6,1) * t104 + rSges(6,2) * t103;
t123 = -rSges(6,1) * t153 + rSges(6,2) * t154 + t89 * rSges(6,3);
t119 = Icges(6,2) * t103 - t144;
t97 = rSges(2,1) * t116 - t113 * rSges(2,2);
t96 = -t113 * rSges(2,1) - rSges(2,2) * t116;
t84 = -rSges(6,1) * t103 - rSges(6,2) * t104;
t80 = t88 * pkin(7);
t76 = t89 * t165;
t74 = rSges(3,1) * t116 + t113 * rSges(3,3) + t138;
t73 = rSges(3,3) * t116 + t106 + (-rSges(3,1) - pkin(1)) * t113;
t53 = -rSges(4,1) * t88 - rSges(4,2) * t89 + t134;
t52 = rSges(4,1) * t89 - rSges(4,2) * t88 + t132;
t51 = (-t84 + t165) * t88;
t50 = -t84 * t89 + t76;
t37 = t152 + t80 + (pkin(3) - t102) * t89;
t36 = t88 * (pkin(3) * t88 + t159 - t167);
t35 = (t161 + t165) * t88;
t34 = t161 * t89 + t76;
t33 = t167 + (-pkin(3) + t155) * t88 + t134 + t158;
t32 = t163 + t80 + (pkin(3) - t129) * t89 + t132;
t30 = -rSges(7,3) * t150 - t131;
t23 = t123 + t130;
t22 = (rSges(6,3) - t110) * t88 + (t102 - t128) * t89 + t132;
t21 = (-t103 * t151 + t162) * t104;
t20 = t88 * (t88 * t155 + t158) + (t129 * t89 - t163) * t89;
t19 = t104 * t31 + t67 * t154;
t18 = -t104 * t30 - t67 * t150;
t15 = t130 + t170;
t14 = -t152 + (t164 + t102 + (rSges(7,3) + pkin(8)) * t103) * t89 + t131 + t132;
t13 = (-t30 * t88 + t31 * t89) * t103;
t12 = t36 + t88 * t123 + (-t88 * rSges(6,3) + t128 * t89 + t37) * t89;
t9 = -t25 * t154 + t27 * t58 + t29 * t59;
t8 = -t24 * t154 + t26 * t58 + t28 * t59;
t7 = -t25 * t150 + t27 * t56 + t29 * t57;
t6 = -t24 * t150 + t26 * t56 + t28 * t57;
t5 = t36 + t170 * t88 + (t37 + t30 + (-pkin(8) * t103 - t164) * t89) * t89;
t4 = -t8 * t88 + t89 * t9;
t3 = -t6 * t88 + t7 * t89;
t2 = t104 * t17 + (-t8 * t89 - t88 * t9) * t103;
t1 = t104 * t16 + (-t6 * t89 - t7 * t88) * t103;
t38 = [-t104 * (-Icges(6,4) * t103 - Icges(6,2) * t104) + t115 ^ 2 * Icges(5,2) + Icges(3,2) + Icges(2,3) + Icges(4,3) + (Icges(6,1) * t103 + t144 - t151) * t103 + m(7) * (t14 ^ 2 + t15 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2) + m(6) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t52 ^ 2 + t53 ^ 2) + m(3) * (t73 ^ 2 + t74 ^ 2) + m(2) * (t96 ^ 2 + t97 ^ 2) + t162 + (Icges(5,1) * t112 + 0.2e1 * Icges(5,4) * t115) * t112; m(7) * (t113 * t14 - t116 * t15) + m(5) * (t113 * t32 - t116 * t33) + m(6) * (t113 * t22 - t116 * t23) + m(4) * (t113 * t52 - t116 * t53) + m(3) * (t113 * t73 - t116 * t74); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + t137) * (t113 ^ 2 + t116 ^ 2); 0; 0; m(4) + m(5) + m(6) + m(7); (-t119 * t153 / 0.2e1 - t135 + (-t183 / 0.2e1 + t182) * t89) * t89 + (t119 * t89 * t166 + t136 + (-Icges(6,6) * t166 + t182) * t88) * t88 + m(7) * (t14 * t35 + t15 * t34) + m(6) * (t22 * t51 + t23 * t50) + (-t32 * t88 - t33 * t89) * t168 + (-t179 - t180 - t181 - t183) * (t86 / 0.2e1 + t85 / 0.2e1); m(6) * (t51 * t113 - t116 * t50) + m(7) * (t35 * t113 - t116 * t34) + (-t113 * t88 + t116 * t89) * t168; -m(5) * t20 - m(6) * t12 - m(7) * t5; m(7) * (t34 ^ 2 + t35 ^ 2 + t5 ^ 2) + m(6) * (t12 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t157 * t95 ^ 2 + t20 ^ 2) + (t172 * t85 + t4) * t89 + (-t3 + t171 * t86 + (t171 * t89 + t172 * t88) * t89) * t88; m(7) * (t14 * t89 - t15 * t88) + m(6) * (t22 * t89 - t23 * t88); (t89 * t113 + t116 * t88) * t173; 0; m(7) * (-t34 * t88 + t35 * t89) + m(6) * (-t50 * t88 + t51 * t89); t157 * t173; m(7) * (t14 * t18 + t15 * t19) + t21 + (t135 * t88 + t136 * t89) * t103; m(7) * (t18 * t113 - t116 * t19); -m(7) * t13; t89 * t2 / 0.2e1 + t1 * t169 + (-t10 * t88 + t11 * t89) * t166 + m(7) * (t13 * t5 + t18 * t35 + t19 * t34) + (t4 * t169 - t89 * t3 / 0.2e1) * t103; m(7) * (t18 * t89 - t19 * t88); t104 * t21 + m(7) * (t13 ^ 2 + t18 ^ 2 + t19 ^ 2) + (-t88 * t2 - t89 * t1 + t104 * (-t10 * t89 - t11 * t88)) * t103;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t38(1) t38(2) t38(4) t38(7) t38(11) t38(16); t38(2) t38(3) t38(5) t38(8) t38(12) t38(17); t38(4) t38(5) t38(6) t38(9) t38(13) t38(18); t38(7) t38(8) t38(9) t38(10) t38(14) t38(19); t38(11) t38(12) t38(13) t38(14) t38(15) t38(20); t38(16) t38(17) t38(18) t38(19) t38(20) t38(21);];
Mq  = res;
