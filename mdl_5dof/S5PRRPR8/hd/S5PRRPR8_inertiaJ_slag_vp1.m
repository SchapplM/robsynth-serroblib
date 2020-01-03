% Calculate joint inertia matrix for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:19
% EndTime: 2019-12-31 17:42:22
% DurationCPUTime: 1.23s
% Computational Cost: add. (3943->184), mult. (3866->298), div. (0->0), fcn. (4038->10), ass. (0->100)
t166 = Icges(4,3) + Icges(5,3);
t105 = qJ(2) + qJ(3);
t100 = cos(t105);
t98 = pkin(9) + t105;
t95 = sin(t98);
t96 = cos(t98);
t99 = sin(t105);
t165 = Icges(4,5) * t100 + Icges(5,5) * t96 - Icges(4,6) * t99 - Icges(5,6) * t95;
t106 = sin(pkin(8));
t107 = cos(pkin(8));
t164 = -t165 * t106 + t166 * t107;
t163 = t166 * t106 + t165 * t107;
t103 = t106 ^ 2;
t104 = t107 ^ 2;
t162 = t103 + t104;
t120 = Icges(4,4) * t100 - Icges(4,2) * t99;
t121 = Icges(4,1) * t100 - Icges(4,4) * t99;
t125 = Icges(5,4) * t96 - Icges(5,2) * t95;
t126 = Icges(5,1) * t96 - Icges(5,4) * t95;
t161 = -t100 * (-Icges(4,5) * t107 + t121 * t106) + (-Icges(4,6) * t107 + t120 * t106) * t99 + (-Icges(5,6) * t107 + t125 * t106) * t95 - (-Icges(5,5) * t107 + t126 * t106) * t96;
t160 = (Icges(5,6) * t106 + t125 * t107) * t95 - (Icges(5,5) * t106 + t126 * t107) * t96 - t100 * (Icges(4,5) * t106 + t121 * t107) + (Icges(4,6) * t106 + t120 * t107) * t99;
t108 = sin(qJ(5));
t110 = cos(qJ(5));
t56 = -t96 * rSges(6,3) + (rSges(6,1) * t110 - rSges(6,2) * t108) * t95;
t158 = -t95 * pkin(4) + t96 * pkin(7) - t56;
t148 = t106 * t95;
t143 = t107 * t110;
t146 = t106 * t108;
t83 = -t96 * t146 - t143;
t144 = t107 * t108;
t145 = t106 * t110;
t84 = t96 * t145 - t144;
t36 = Icges(6,5) * t84 + Icges(6,6) * t83 + Icges(6,3) * t148;
t38 = Icges(6,4) * t84 + Icges(6,2) * t83 + Icges(6,6) * t148;
t40 = Icges(6,1) * t84 + Icges(6,4) * t83 + Icges(6,5) * t148;
t15 = t36 * t148 + t83 * t38 + t84 * t40;
t147 = t107 * t95;
t85 = -t96 * t144 + t145;
t86 = t96 * t143 + t146;
t37 = Icges(6,5) * t86 + Icges(6,6) * t85 + Icges(6,3) * t147;
t39 = Icges(6,4) * t86 + Icges(6,2) * t85 + Icges(6,6) * t147;
t41 = Icges(6,1) * t86 + Icges(6,4) * t85 + Icges(6,5) * t147;
t16 = t37 * t148 + t83 * t39 + t84 * t41;
t8 = t16 * t106 - t15 * t107;
t157 = -t8 + t164 * t104 + (t160 * t106 + (-t161 + t163) * t107) * t106;
t156 = pkin(3) * t99;
t109 = sin(qJ(2));
t154 = pkin(2) * t109;
t51 = -Icges(6,3) * t96 + (Icges(6,5) * t110 - Icges(6,6) * t108) * t95;
t153 = t96 * t51;
t151 = t162 * pkin(3) * t100;
t111 = cos(qJ(2));
t150 = t162 * t111 * pkin(2);
t35 = t162 * (rSges(4,1) * t100 - rSges(4,2) * t99);
t17 = t36 * t147 + t85 * t38 + t86 * t40;
t18 = t37 * t147 + t85 * t39 + t86 * t41;
t9 = t18 * t106 - t17 * t107;
t140 = (t9 + t163 * t103 + ((-t160 + t164) * t106 + t161 * t107) * t107) * t106;
t87 = t95 * rSges(5,1) + t96 * rSges(5,2);
t139 = -t87 - t156;
t90 = t99 * rSges(4,1) + t100 * rSges(4,2);
t138 = -t90 - t154;
t24 = t151 + t162 * (rSges(5,1) * t96 - rSges(5,2) * t95);
t22 = -t96 * t36 + (-t108 * t38 + t110 * t40) * t95;
t23 = -t96 * t37 + (-t108 * t39 + t110 * t41) * t95;
t52 = -Icges(6,6) * t96 + (Icges(6,4) * t110 - Icges(6,2) * t108) * t95;
t53 = -Icges(6,5) * t96 + (Icges(6,1) * t110 - Icges(6,4) * t108) * t95;
t3 = -(t83 * t52 + t84 * t53) * t96 + (t16 * t107 + (t15 - t153) * t106) * t95;
t4 = -(t85 * t52 + t86 * t53) * t96 + (t17 * t106 + (t18 - t153) * t107) * t95;
t137 = t106 * t4 / 0.2e1 - t96 * (t23 * t106 - t22 * t107) / 0.2e1 - t107 * t3 / 0.2e1 + t8 * t148 / 0.2e1 + t9 * t147 / 0.2e1;
t136 = -t156 + t158;
t134 = -t154 - t156;
t42 = t84 * rSges(6,1) + t83 * rSges(6,2) + rSges(6,3) * t148;
t43 = t86 * rSges(6,1) + t85 * rSges(6,2) + rSges(6,3) * t147;
t12 = t106 * t42 + t107 * t43 + t151 + t162 * (pkin(4) * t96 + pkin(7) * t95);
t116 = Icges(3,5) * t111 - Icges(3,6) * t109;
t115 = t134 - t87;
t114 = t134 + t158;
t113 = t157 * t107 + t140;
t93 = t109 * rSges(3,1) + t111 * rSges(3,2);
t78 = Icges(3,3) * t106 + t116 * t107;
t77 = -Icges(3,3) * t107 + t116 * t106;
t70 = t138 * t107;
t69 = t138 * t106;
t55 = t139 * t107;
t54 = t139 * t106;
t48 = t115 * t107;
t47 = t115 * t106;
t44 = t162 * (rSges(3,1) * t111 - rSges(3,2) * t109);
t32 = t136 * t107;
t31 = t136 * t106;
t30 = t114 * t107;
t29 = t114 * t106;
t28 = -t56 * t147 - t96 * t43;
t27 = t56 * t148 + t96 * t42;
t26 = t35 + t150;
t25 = (-t106 * t43 + t107 * t42) * t95;
t19 = t24 + t150;
t11 = t12 + t150;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t44 + m(4) * t26 + m(5) * t19 + m(6) * t11; m(6) * (t11 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(5) * (t19 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(4) * (t26 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(3) * (t162 * t93 ^ 2 + t44 ^ 2) + t106 * t103 * t78 + t140 + (-t104 * t77 + (-t106 * t77 + t107 * t78) * t106 + t157) * t107; m(4) * t35 + m(5) * t24 + m(6) * t12; m(6) * (t12 * t11 + t31 * t29 + t32 * t30) + m(5) * (t24 * t19 + t54 * t47 + t55 * t48) + m(4) * (t35 * t26 + (-t106 * t69 - t107 * t70) * t90) + t113; m(6) * (t12 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(5) * (t24 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(4) * (t162 * t90 ^ 2 + t35 ^ 2) + t113; 0; m(6) * (t106 * t30 - t107 * t29) + m(5) * (t106 * t48 - t107 * t47); m(6) * (t106 * t32 - t107 * t31) + m(5) * (t106 * t55 - t107 * t54); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t162; m(6) * t25; m(6) * (t25 * t11 + t27 * t30 + t28 * t29) + t137; m(6) * (t25 * t12 + t27 * t32 + t28 * t31) + t137; m(6) * (t27 * t106 - t28 * t107); m(6) * (t25 ^ 2 + t27 ^ 2 + t28 ^ 2) + t4 * t147 + t3 * t148 - t96 * (t96 ^ 2 * t51 + (t23 * t107 + t22 * t106 - (-t108 * t52 + t110 * t53) * t96) * t95);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
