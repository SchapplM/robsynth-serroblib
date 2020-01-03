% Calculate joint inertia matrix for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR15_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR15_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR15_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:38
% EndTime: 2019-12-31 18:36:42
% DurationCPUTime: 1.28s
% Computational Cost: add. (1990->262), mult. (3109->400), div. (0->0), fcn. (3225->8), ass. (0->129)
t112 = sin(qJ(1));
t165 = -t112 / 0.2e1;
t114 = cos(qJ(1));
t153 = t114 / 0.2e1;
t164 = rSges(5,3) + qJ(4);
t100 = t114 * qJ(2);
t111 = sin(qJ(3));
t108 = sin(pkin(8));
t138 = t114 * t108;
t109 = cos(pkin(8));
t141 = t112 * t109;
t75 = t111 * t138 + t141;
t137 = t114 * t109;
t142 = t112 * t108;
t76 = -t111 * t137 + t142;
t128 = -t76 * rSges(5,1) - t75 * rSges(5,2);
t113 = cos(qJ(3));
t129 = t164 * t113;
t155 = -pkin(1) - pkin(6);
t143 = t111 * t114;
t95 = pkin(3) * t143;
t21 = t155 * t112 - t114 * t129 + t100 + t128 + t95;
t147 = t114 * pkin(1) + t112 * qJ(2);
t130 = t114 * pkin(6) + t147;
t144 = t111 * t112;
t73 = -t111 * t142 + t137;
t74 = t111 * t141 + t138;
t161 = -t74 * rSges(5,1) - t73 * rSges(5,2) - pkin(3) * t144;
t22 = -t112 * t129 + t130 - t161;
t163 = m(5) * (t112 * t21 - t114 * t22);
t110 = -pkin(7) - qJ(4);
t104 = pkin(8) + qJ(5);
t97 = sin(t104);
t98 = cos(t104);
t59 = t112 * t98 + t97 * t143;
t60 = t112 * t97 - t98 * t143;
t127 = -t60 * rSges(6,1) - t59 * rSges(6,2);
t96 = t109 * pkin(4) + pkin(3);
t148 = t111 * t96;
t17 = t100 + (t148 + (-rSges(6,3) + t110) * t113) * t114 + (-pkin(4) * t108 + t155) * t112 + t127;
t140 = t112 * t113;
t57 = t114 * t98 - t97 * t144;
t58 = t114 * t97 + t98 * t144;
t31 = t58 * rSges(6,1) + t57 * rSges(6,2) - rSges(6,3) * t140;
t160 = -pkin(4) * t138 - t110 * t140 - t96 * t144 - t31;
t18 = t130 - t160;
t162 = m(6) * (t112 * t17 - t114 * t18);
t159 = (rSges(4,1) * t111 + rSges(4,2) * t113) * t114;
t158 = -qJ(4) - t110;
t105 = t112 ^ 2;
t107 = t114 ^ 2;
t54 = Icges(5,6) * t111 + (Icges(5,4) * t109 - Icges(5,2) * t108) * t113;
t157 = t54 / 0.2e1;
t55 = Icges(5,5) * t111 + (Icges(5,1) * t109 - Icges(5,4) * t108) * t113;
t156 = t55 / 0.2e1;
t154 = t112 / 0.2e1;
t48 = Icges(6,6) * t111 + (Icges(6,4) * t98 - Icges(6,2) * t97) * t113;
t152 = t97 * t48;
t47 = Icges(6,3) * t111 + (Icges(6,5) * t98 - Icges(6,6) * t97) * t113;
t49 = Icges(6,5) * t111 + (Icges(6,1) * t98 - Icges(6,4) * t97) * t113;
t151 = t113 * t98 * t49 + t111 * t47;
t50 = t111 * rSges(6,3) + (rSges(6,1) * t98 - rSges(6,2) * t97) * t113;
t150 = (-pkin(3) + t96) * t113 + t158 * t111 + t50;
t146 = Icges(4,4) * t111;
t145 = Icges(4,4) * t113;
t139 = t113 * t114;
t136 = t105 + t107;
t135 = m(5) / 0.2e1 + m(6) / 0.2e1;
t133 = rSges(4,1) * t144 + rSges(4,2) * t140 + t114 * rSges(4,3);
t25 = Icges(6,5) * t58 + Icges(6,6) * t57 - Icges(6,3) * t140;
t27 = Icges(6,4) * t58 + Icges(6,2) * t57 - Icges(6,6) * t140;
t29 = Icges(6,1) * t58 + Icges(6,4) * t57 - Icges(6,5) * t140;
t10 = t111 * t25 + (-t27 * t97 + t29 * t98) * t113;
t13 = -t47 * t140 + t57 * t48 + t58 * t49;
t132 = -t13 / 0.2e1 - t10 / 0.2e1;
t26 = Icges(6,5) * t60 + Icges(6,6) * t59 + Icges(6,3) * t139;
t28 = Icges(6,4) * t60 + Icges(6,2) * t59 + Icges(6,6) * t139;
t30 = Icges(6,1) * t60 + Icges(6,4) * t59 + Icges(6,5) * t139;
t11 = t111 * t26 + (-t28 * t97 + t30 * t98) * t113;
t14 = t47 * t139 + t59 * t48 + t60 * t49;
t131 = t14 / 0.2e1 + t11 / 0.2e1;
t19 = t111 * t31 + t50 * t140;
t32 = rSges(6,3) * t139 - t127;
t20 = -t111 * t32 + t50 * t139;
t122 = t20 * t112 - t19 * t114;
t86 = t113 * pkin(3) + t111 * qJ(4);
t78 = t112 * t86;
t23 = t150 * t112 + t78;
t24 = (-t86 - t150) * t114;
t120 = t23 * t112 - t24 * t114;
t56 = t111 * rSges(5,3) + (rSges(5,1) * t109 - rSges(5,2) * t108) * t113;
t40 = t112 * t56 + t78;
t41 = (-t56 - t86) * t114;
t119 = t40 * t112 - t41 * t114;
t118 = Icges(4,1) * t111 + t145;
t117 = Icges(4,2) * t113 + t146;
t116 = Icges(4,5) * t111 + Icges(4,6) * t113;
t43 = t100 + t159 + (-rSges(4,3) + t155) * t112;
t44 = t130 + t133;
t115 = m(4) * (t112 * t43 - t114 * t44);
t88 = t114 * rSges(2,1) - t112 * rSges(2,2);
t87 = t113 * rSges(4,1) - t111 * rSges(4,2);
t85 = -t112 * rSges(2,1) - t114 * rSges(2,2);
t82 = Icges(4,5) * t113 - Icges(4,6) * t111;
t69 = t114 * (qJ(4) * t139 - t95);
t68 = -t114 * rSges(3,2) + t112 * rSges(3,3) + t147;
t67 = t114 * rSges(3,3) + t100 + (rSges(3,2) - pkin(1)) * t112;
t62 = Icges(4,3) * t112 - t116 * t114;
t61 = Icges(4,3) * t114 + t116 * t112;
t39 = Icges(5,1) * t76 + Icges(5,4) * t75 + Icges(5,5) * t139;
t38 = Icges(5,1) * t74 + Icges(5,4) * t73 - Icges(5,5) * t140;
t37 = Icges(5,4) * t76 + Icges(5,2) * t75 + Icges(5,6) * t139;
t36 = Icges(5,4) * t74 + Icges(5,2) * t73 - Icges(5,6) * t140;
t35 = Icges(5,5) * t76 + Icges(5,6) * t75 + Icges(5,3) * t139;
t34 = Icges(5,5) * t74 + Icges(5,6) * t73 - Icges(5,3) * t140;
t33 = -t112 * t133 + (t112 * rSges(4,3) - t159) * t114;
t16 = (-t113 * t152 + t151) * t111;
t15 = (-t112 * t32 - t114 * t31) * t113;
t12 = t69 + t114 * (rSges(5,3) * t139 - t128) + (t164 * t140 + t161) * t112;
t9 = t26 * t139 + t59 * t28 + t60 * t30;
t8 = t25 * t139 + t59 * t27 + t60 * t29;
t7 = -t26 * t140 + t57 * t28 + t58 * t30;
t6 = -t25 * t140 + t57 * t27 + t58 * t29;
t5 = t69 + t160 * t112 + (pkin(4) * t142 + t32 + t95 + (t158 * t113 - t148) * t114) * t114;
t4 = t9 * t112 + t8 * t114;
t3 = t7 * t112 + t6 * t114;
t2 = t14 * t111 + (-t112 * t8 + t114 * t9) * t113;
t1 = t13 * t111 + (-t112 * t6 + t114 * t7) * t113;
t42 = [Icges(3,1) + Icges(2,3) + (-t145 + (Icges(5,5) * t109 - Icges(5,6) * t108) * t113 + (Icges(4,2) + Icges(5,3)) * t111) * t111 + (Icges(4,1) * t113 - t108 * t54 + t109 * t55 - t146 - t152) * t113 + m(6) * (t17 ^ 2 + t18 ^ 2) + m(5) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t43 ^ 2 + t44 ^ 2) + m(3) * (t67 ^ 2 + t68 ^ 2) + m(2) * (t85 ^ 2 + t88 ^ 2) + t151; t162 + t163 + t115 + m(3) * (t112 * t67 - t114 * t68); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t135) * t136; (t82 * t153 + t74 * t156 + t73 * t157 - t132) * t114 + (t82 * t154 + t76 * t156 + t75 * t157 + t131) * t112 + m(6) * (t23 * t17 + t24 * t18) + m(5) * (t40 * t21 + t41 * t22) + t87 * t115 + ((-Icges(4,6) * t114 / 0.2e1 + t117 * t165 + t34 / 0.2e1) * t114 + (Icges(4,6) * t165 + t117 * t153 + t35 / 0.2e1) * t112) * t111 + ((Icges(4,5) * t112 - t108 * t37 + t109 * t39 - t118 * t114) * t154 + (Icges(4,5) * t114 - t108 * t36 + t109 * t38 + t118 * t112) * t153) * t113; m(4) * t136 * t87 + m(5) * t119 + m(6) * t120; m(6) * (t23 ^ 2 + t24 ^ 2 + t5 ^ 2) + m(5) * (t12 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(4) * (t136 * t87 ^ 2 + t33 ^ 2) + (t107 * t61 + t3 + (-t34 * t140 + t73 * t36 + t74 * t38) * t114) * t114 + (t4 + t105 * t62 + (t35 * t139 + t75 * t37 + t76 * t39) * t112 + (t112 * t61 + t114 * t62 + t34 * t139 - t35 * t140 + t75 * t36 + t73 * t37 + t76 * t38 + t74 * t39) * t114) * t112; 0.2e1 * (-t162 / 0.2e1 - t163 / 0.2e1) * t113; -0.2e1 * t135 * t136 * t113; m(6) * (t111 * t5 - t120 * t113) + m(5) * (t111 * t12 - t119 * t113); 0.2e1 * t135 * (t136 * t113 ^ 2 + t111 ^ 2); m(6) * (t20 * t17 + t19 * t18) + t16 + (t132 * t112 + t131 * t114) * t113; m(6) * t122; t1 * t153 + t111 * (t10 * t114 + t11 * t112) / 0.2e1 + m(6) * (t15 * t5 + t19 * t24 + t20 * t23) + t2 * t154 + (t4 * t153 + t3 * t165) * t113; m(6) * (t15 * t111 - t122 * t113); m(6) * (t15 ^ 2 + t19 ^ 2 + t20 ^ 2) + t111 * t16 + (-t112 * t1 + t114 * t2 + t111 * (-t10 * t112 + t11 * t114)) * t113;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t42(1), t42(2), t42(4), t42(7), t42(11); t42(2), t42(3), t42(5), t42(8), t42(12); t42(4), t42(5), t42(6), t42(9), t42(13); t42(7), t42(8), t42(9), t42(10), t42(14); t42(11), t42(12), t42(13), t42(14), t42(15);];
Mq = res;
