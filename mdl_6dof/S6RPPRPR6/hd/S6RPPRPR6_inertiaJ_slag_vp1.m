% Calculate joint inertia matrix for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:43
% EndTime: 2019-03-09 01:50:46
% DurationCPUTime: 1.32s
% Computational Cost: add. (1318->251), mult. (3063->371), div. (0->0), fcn. (3057->6), ass. (0->122)
t166 = Icges(6,1) + Icges(5,3);
t104 = sin(qJ(4));
t107 = cos(qJ(4));
t165 = (-Icges(6,5) + Icges(5,6)) * t107 + (-Icges(6,4) + Icges(5,5)) * t104;
t105 = sin(qJ(1));
t164 = -t105 / 0.2e1;
t163 = t105 / 0.2e1;
t108 = cos(qJ(1));
t162 = -t108 / 0.2e1;
t154 = t108 / 0.2e1;
t161 = t107 / 0.2e1;
t146 = rSges(6,3) * t107;
t139 = t105 * t107;
t90 = qJ(5) * t139;
t98 = t108 * qJ(2);
t150 = t90 + t98;
t151 = -pkin(1) - qJ(3);
t153 = -rSges(6,1) - pkin(7);
t21 = t153 * t108 + (t146 + (rSges(6,2) - pkin(4)) * t104 + t151) * t105 + t150;
t148 = t108 * pkin(1) + t105 * qJ(2);
t132 = t108 * qJ(3) + t148;
t147 = rSges(6,2) * t104;
t140 = t104 * t108;
t95 = pkin(4) * t140;
t22 = t95 + t153 * t105 + (-t147 + (-rSges(6,3) - qJ(5)) * t107) * t108 + t132;
t160 = m(6) * (t105 * t22 + t108 * t21);
t106 = cos(qJ(6));
t103 = sin(qJ(6));
t137 = t108 * t103;
t66 = -t106 * t139 - t137;
t136 = t108 * t106;
t67 = -t103 * t139 + t136;
t126 = -t67 * rSges(7,1) - t66 * rSges(7,2);
t156 = -pkin(5) - pkin(7);
t15 = t156 * t108 + ((-rSges(7,3) - pkin(4) - pkin(8)) * t104 + t151) * t105 + t126 + t150;
t64 = t105 * t103 - t107 * t136;
t65 = -t105 * t106 - t107 * t137;
t32 = t65 * rSges(7,1) + t64 * rSges(7,2) + rSges(7,3) * t140;
t138 = t107 * t108;
t69 = -qJ(5) * t138 + t95;
t134 = -pkin(8) * t140 - t32 - t69;
t16 = t156 * t105 + t132 - t134;
t159 = m(7) * (t105 * t16 + t108 * t15);
t158 = t165 * t105 + t166 * t108;
t157 = t166 * t105 - t165 * t108;
t100 = t105 ^ 2;
t102 = t108 ^ 2;
t87 = t107 * rSges(5,1) - t104 * rSges(5,2);
t155 = m(5) * t87;
t89 = t100 + t102;
t75 = m(5) * t89;
t152 = -rSges(5,3) - pkin(7);
t149 = rSges(5,1) * t140 + rSges(5,2) * t138;
t145 = Icges(5,4) * t104;
t144 = Icges(5,4) * t107;
t143 = Icges(6,6) * t104;
t142 = Icges(6,6) * t107;
t141 = t104 * t105;
t135 = m(6) / 0.2e1 + m(7) / 0.2e1;
t43 = Icges(7,3) * t107 + (Icges(7,5) * t103 + Icges(7,6) * t106) * t104;
t46 = Icges(7,6) * t107 + (Icges(7,4) * t103 + Icges(7,2) * t106) * t104;
t49 = Icges(7,5) * t107 + (Icges(7,1) * t103 + Icges(7,4) * t106) * t104;
t133 = t107 * t43 + (t103 * t49 + t106 * t46) * t104;
t26 = Icges(7,5) * t65 + Icges(7,6) * t64 + Icges(7,3) * t140;
t28 = Icges(7,4) * t65 + Icges(7,2) * t64 + Icges(7,6) * t140;
t30 = Icges(7,1) * t65 + Icges(7,4) * t64 + Icges(7,5) * t140;
t10 = t107 * t26 + (t103 * t30 + t106 * t28) * t104;
t12 = t43 * t140 + t64 * t46 + t65 * t49;
t131 = t12 / 0.2e1 + t10 / 0.2e1;
t27 = Icges(7,5) * t67 + Icges(7,6) * t66 + Icges(7,3) * t141;
t29 = Icges(7,4) * t67 + Icges(7,2) * t66 + Icges(7,6) * t141;
t31 = Icges(7,1) * t67 + Icges(7,4) * t66 + Icges(7,5) * t141;
t11 = t107 * t27 + (t103 * t31 + t106 * t29) * t104;
t13 = t43 * t141 + t46 * t66 + t49 * t67;
t130 = t13 / 0.2e1 + t11 / 0.2e1;
t129 = Icges(5,5) * t161 - Icges(6,4) * t107 / 0.2e1 + (-Icges(5,6) / 0.2e1 + Icges(6,5) / 0.2e1) * t104;
t58 = t107 * rSges(7,3) + (rSges(7,1) * t103 + rSges(7,2) * t106) * t104;
t128 = pkin(8) * t107 + t58;
t127 = t75 + (m(4) + m(6) + m(7)) * t89;
t125 = -rSges(5,1) * t104 - rSges(5,2) * t107;
t124 = t146 + t147;
t33 = rSges(7,3) * t141 - t126;
t19 = -t107 * t33 + t58 * t141;
t20 = t107 * t32 - t58 * t140;
t118 = t20 * t105 + t108 * t19;
t85 = t107 * pkin(4) + t104 * qJ(5);
t70 = t105 * t85;
t24 = t128 * t105 + t70;
t71 = t108 * t85;
t25 = t128 * t108 + t71;
t116 = t24 * t105 + t108 * t25;
t86 = -t107 * rSges(6,2) + t104 * rSges(6,3);
t36 = t105 * t86 + t70;
t37 = t108 * t86 + t71;
t115 = t36 * t105 + t108 * t37;
t114 = Icges(5,1) * t104 + t144;
t113 = Icges(5,2) * t107 + t145;
t110 = -Icges(6,2) * t104 - t142;
t109 = -Icges(6,3) * t107 - t143;
t88 = t108 * rSges(2,1) - t105 * rSges(2,2);
t84 = -t105 * rSges(2,1) - t108 * rSges(2,2);
t68 = pkin(4) * t141 - t90;
t60 = -t108 * rSges(3,2) + t105 * rSges(3,3) + t148;
t59 = t108 * rSges(3,3) + t98 + (rSges(3,2) - pkin(1)) * t105;
t42 = t105 * rSges(4,2) + t108 * rSges(4,3) + t132;
t41 = t108 * rSges(4,2) + t98 + (-rSges(4,3) + t151) * t105;
t35 = t152 * t105 + t132 + t149;
t34 = t98 + t152 * t108 + (t125 + t151) * t105;
t23 = t125 * t100 - t108 * t149;
t18 = t133 * t107;
t17 = (t124 * t108 - t69) * t108 + (t124 * t105 - t68) * t105;
t14 = (-t105 * t32 + t108 * t33) * t104;
t9 = t134 * t108 + (-pkin(8) * t141 - t33 - t68) * t105;
t8 = t27 * t141 + t29 * t66 + t31 * t67;
t7 = t26 * t141 + t28 * t66 + t30 * t67;
t6 = t27 * t140 + t64 * t29 + t65 * t31;
t5 = t26 * t140 + t64 * t28 + t65 * t30;
t4 = -t7 * t105 + t108 * t8;
t3 = -t5 * t105 + t108 * t6;
t2 = t13 * t107 + (t105 * t8 + t108 * t7) * t104;
t1 = t12 * t107 + (t105 * t6 + t108 * t5) * t104;
t38 = [Icges(3,1) + Icges(4,1) + Icges(2,3) + m(7) * (t15 ^ 2 + t16 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2) + m(5) * (t34 ^ 2 + t35 ^ 2) + m(4) * (t41 ^ 2 + t42 ^ 2) + m(3) * (t59 ^ 2 + t60 ^ 2) + m(2) * (t84 ^ 2 + t88 ^ 2) + t133 + (-t143 - t145 + (Icges(5,1) + Icges(6,2)) * t107) * t107 + (-t142 - t144 + (Icges(5,2) + Icges(6,3)) * t104) * t104; m(7) * (t105 * t15 - t108 * t16) + m(6) * (t105 * t21 - t108 * t22) + m(5) * (t105 * t34 - t108 * t35) + m(4) * (t105 * t41 - t108 * t42) + m(3) * (t105 * t59 - t108 * t60); m(3) * t89 + t127; t159 + t160 + m(5) * (t105 * t35 + t108 * t34) + m(4) * (t105 * t42 + t108 * t41); 0; t127; m(7) * (t15 * t25 + t16 * t24) + m(6) * (t21 * t37 + t22 * t36) + (t34 * t155 + t129 * t108 + (Icges(6,4) * t162 + Icges(5,5) * t154 + t110 * t164 + t114 * t163) * t107 + t130) * t108 + (t35 * t155 + t129 * t105 + (Icges(6,4) * t164 + Icges(5,5) * t163 + t110 * t154 + t114 * t162) * t107 - t131) * t105 + ((Icges(6,5) * t154 + Icges(5,6) * t162 + t109 * t163 + t113 * t164) * t108 + (Icges(6,5) * t163 + Icges(5,6) * t164 + t109 * t162 + t113 * t154) * t105) * t104; m(6) * (t37 * t105 - t108 * t36) + m(7) * (t25 * t105 - t108 * t24); m(6) * t115 + m(7) * t116 + t87 * t75; m(7) * (t24 ^ 2 + t25 ^ 2 + t9 ^ 2) + m(6) * (t17 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(5) * (t89 * t87 ^ 2 + t23 ^ 2) + (t158 * t102 + t4) * t108 + (-t3 + t157 * t100 + (t158 * t105 + t157 * t108) * t108) * t105; 0.2e1 * (-t159 / 0.2e1 - t160 / 0.2e1) * t107; 0; -0.2e1 * t135 * t89 * t107; m(7) * (t104 * t9 - t116 * t107) + m(6) * (t104 * t17 - t115 * t107); 0.2e1 * t135 * (t89 * t107 ^ 2 + t104 ^ 2); m(7) * (t15 * t19 + t16 * t20) + t18 + (t130 * t105 + t131 * t108) * t104; m(7) * (t19 * t105 - t108 * t20); m(7) * t118; t2 * t154 + (-t10 * t105 + t11 * t108) * t161 + m(7) * (t14 * t9 + t19 * t25 + t20 * t24) + t1 * t164 + (t3 * t154 + t4 * t163) * t104; m(7) * (t14 * t104 - t118 * t107); t107 * t18 + m(7) * (t14 ^ 2 + t19 ^ 2 + t20 ^ 2) + (t108 * t1 + t105 * t2 + t107 * (t10 * t108 + t105 * t11)) * t104;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t38(1) t38(2) t38(4) t38(7) t38(11) t38(16); t38(2) t38(3) t38(5) t38(8) t38(12) t38(17); t38(4) t38(5) t38(6) t38(9) t38(13) t38(18); t38(7) t38(8) t38(9) t38(10) t38(14) t38(19); t38(11) t38(12) t38(13) t38(14) t38(15) t38(20); t38(16) t38(17) t38(18) t38(19) t38(20) t38(21);];
Mq  = res;
