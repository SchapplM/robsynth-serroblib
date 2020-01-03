% Calculate joint inertia matrix for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR11_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR11_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:11
% EndTime: 2019-12-31 18:27:15
% DurationCPUTime: 1.49s
% Computational Cost: add. (1902->184), mult. (2396->288), div. (0->0), fcn. (2469->8), ass. (0->92)
t84 = pkin(8) + qJ(3);
t81 = cos(t84);
t162 = t81 ^ 2;
t160 = Icges(5,4) + Icges(4,5);
t159 = Icges(4,6) - Icges(5,6);
t155 = Icges(5,2) + Icges(4,3);
t80 = sin(t84);
t154 = -t159 * t80 + t160 * t81;
t152 = (-Icges(5,6) / 0.2e1 + Icges(4,6) / 0.2e1) * t81 + (Icges(4,5) / 0.2e1 + Icges(5,4) / 0.2e1) * t80;
t90 = sin(qJ(5));
t92 = cos(qJ(5));
t57 = t80 * t92 - t81 * t90;
t151 = t57 / 0.2e1;
t105 = t80 * t90 + t81 * t92;
t144 = -t105 / 0.2e1;
t143 = Icges(6,5) * t151 + Icges(6,6) * t144;
t91 = sin(qJ(1));
t93 = cos(qJ(1));
t142 = -t154 * t91 + t155 * t93;
t141 = t154 * t93 + t155 * t91;
t85 = t91 ^ 2;
t86 = t93 ^ 2;
t140 = m(5) / 0.2e1;
t139 = t80 * t93;
t138 = t81 * t93;
t54 = t57 * t93;
t55 = t105 * t93;
t137 = t55 * rSges(6,1) + t54 * rSges(6,2);
t123 = qJ(4) * t80;
t134 = pkin(3) * t138 + t93 * t123;
t136 = t85 * (pkin(3) * t81 + t123) + t93 * t134;
t65 = t80 * pkin(3) - t81 * qJ(4);
t135 = -t80 * rSges(5,1) + t81 * rSges(5,3) - t65;
t133 = t85 + t86;
t53 = t105 * t91;
t132 = Icges(6,1) * t53;
t129 = Icges(6,4) * t53;
t126 = Icges(6,5) * t93;
t52 = t57 * t91;
t125 = Icges(6,2) * t52;
t124 = Icges(6,6) * t93;
t122 = rSges(3,3) + qJ(2);
t89 = -pkin(6) - qJ(2);
t121 = -rSges(6,3) - pkin(7) - t89;
t120 = t140 + m(6) / 0.2e1;
t28 = Icges(6,4) * t57 - Icges(6,2) * t105;
t29 = Icges(6,1) * t57 - Icges(6,4) * t105;
t94 = t124 + t125 + t129;
t95 = Icges(6,4) * t52 + t126 + t132;
t119 = t93 * t143 + t52 * t28 / 0.2e1 + t53 * t29 / 0.2e1 + t94 * t144 + t95 * t151;
t20 = Icges(6,4) * t55 + Icges(6,2) * t54 - Icges(6,6) * t91;
t21 = Icges(6,1) * t55 + Icges(6,4) * t54 - Icges(6,5) * t91;
t118 = t91 * t143 - t54 * t28 / 0.2e1 - t55 * t29 / 0.2e1 + t105 * t20 / 0.2e1 - t57 * t21 / 0.2e1;
t117 = rSges(5,1) * t138 + t91 * rSges(5,2) + rSges(5,3) * t139;
t88 = cos(pkin(8));
t78 = t88 * pkin(2) + pkin(1);
t71 = t93 * t78;
t116 = -t91 * t89 + t71;
t30 = t57 * rSges(6,1) - rSges(6,2) * t105;
t115 = -pkin(4) * t80 - t30 - t65;
t19 = Icges(6,5) * t55 + Icges(6,6) * t54 - Icges(6,3) * t91;
t113 = t93 * ((t93 * t19 + t52 * t20 + t53 * t21) * t91 - (Icges(6,3) * t86 + (0.2e1 * t126 + t132) * t53 + (0.2e1 * t124 + t125 + 0.2e1 * t129) * t52) * t93) - t91 * ((-t91 * t19 + t54 * t20 + t55 * t21) * t91 - (t55 * t95 + t54 * t94 - t91 * (Icges(6,5) * t53 + Icges(6,6) * t52 + Icges(6,3) * t93)) * t93);
t112 = rSges(4,1) * t81 - rSges(4,2) * t80;
t111 = -t53 * rSges(6,1) - t52 * rSges(6,2);
t14 = t115 * t91;
t15 = t115 * t93;
t110 = t14 * t91 + t15 * t93;
t98 = rSges(4,1) * t138 - rSges(4,2) * t139 + t91 * rSges(4,3);
t87 = sin(pkin(8));
t97 = rSges(3,1) * t88 - rSges(3,2) * t87 + pkin(1);
t11 = t121 * t93 + (-t123 - t78 + (-pkin(3) - pkin(4)) * t81) * t91 + t111;
t75 = pkin(4) * t138;
t12 = t121 * t91 + t134 + t137 + t71 + t75;
t96 = m(6) * (t11 * t93 + t12 * t91);
t69 = t93 * rSges(2,1) - t91 * rSges(2,2);
t68 = -t91 * rSges(2,1) - t93 * rSges(2,2);
t67 = t80 * rSges(4,1) + t81 * rSges(4,2);
t34 = t122 * t91 + t97 * t93;
t33 = t122 * t93 - t97 * t91;
t32 = t135 * t93;
t31 = t135 * t91;
t25 = t116 + t98;
t24 = (rSges(4,3) - t89) * t93 + (-t112 - t78) * t91;
t23 = -t91 * rSges(6,3) + t137;
t22 = t93 * rSges(6,3) - t111;
t18 = t93 * t98 + (-t93 * rSges(4,3) + t112 * t91) * t91;
t17 = t116 + t117 + t134;
t16 = (rSges(5,2) - t89) * t93 + (-t78 + (-rSges(5,1) - pkin(3)) * t81 + (-rSges(5,3) - qJ(4)) * t80) * t91;
t13 = t93 * t117 + (-t93 * rSges(5,2) + (rSges(5,1) * t81 + rSges(5,3) * t80) * t91) * t91 + t136;
t10 = -t91 * t22 - t93 * t23;
t5 = (t23 + t75) * t93 + (t91 * t81 * pkin(4) + t22) * t91 + t136;
t1 = [Icges(3,2) * t88 ^ 2 - t105 * t28 + t57 * t29 + Icges(2,3) + (Icges(3,1) * t87 + 0.2e1 * Icges(3,4) * t88) * t87 + m(6) * (t11 ^ 2 + t12 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(3) * (t33 ^ 2 + t34 ^ 2) + m(2) * (t68 ^ 2 + t69 ^ 2) + (Icges(4,2) + Icges(5,3)) * t162 + (0.2e1 * (Icges(4,4) - Icges(5,5)) * t81 + (Icges(4,1) + Icges(5,1)) * t80) * t80; m(6) * (t91 * t11 - t93 * t12) + m(5) * (t91 * t16 - t93 * t17) + m(4) * (t91 * t24 - t93 * t25) + m(3) * (t91 * t33 - t93 * t34); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t120) * t133; (t152 * t93 - t119) * t93 + (t152 * t91 - t118) * t91 + m(6) * (t15 * t11 + t14 * t12) + m(5) * (t32 * t16 + t31 * t17) + m(4) * (-t24 * t93 - t25 * t91) * t67 + (t159 * t81 + t160 * t80) * (t85 / 0.2e1 + t86 / 0.2e1); m(5) * (-t31 * t93 + t32 * t91) + m(6) * (-t14 * t93 + t15 * t91); m(6) * (t14 ^ 2 + t15 ^ 2 + t5 ^ 2) + m(5) * (t13 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(4) * (t133 * t67 ^ 2 + t18 ^ 2) - t113 + t142 * t93 * t86 + (t141 * t85 + (t141 * t93 + t142 * t91) * t93) * t91; 0.2e1 * (t96 / 0.2e1 + (t16 * t93 + t17 * t91) * t140) * t80; 0; m(6) * (t110 * t80 - t81 * t5) + m(5) * (-t81 * t13 + (t31 * t91 + t32 * t93) * t80); 0.2e1 * t120 * (t133 * t80 ^ 2 + t162); t118 * t91 + t119 * t93 + t30 * t96; 0; m(6) * (t10 * t5 + t110 * t30) + t113; m(6) * (t133 * t80 * t30 - t10 * t81); m(6) * (t133 * t30 ^ 2 + t10 ^ 2) - t113;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
