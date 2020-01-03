% Calculate joint inertia matrix for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR5_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR5_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:58
% EndTime: 2019-12-31 16:23:00
% DurationCPUTime: 0.61s
% Computational Cost: add. (1210->115), mult. (1696->198), div. (0->0), fcn. (1803->8), ass. (0->70)
t62 = sin(pkin(6));
t59 = t62 ^ 2;
t63 = cos(pkin(6));
t60 = t63 ^ 2;
t87 = t59 + t60;
t61 = qJ(2) + pkin(7);
t57 = sin(t61);
t58 = cos(t61);
t66 = sin(qJ(2));
t68 = cos(qJ(2));
t103 = Icges(3,5) * t68 + Icges(4,5) * t58 - Icges(3,6) * t66 - Icges(4,6) * t57;
t102 = Icges(3,3) + Icges(4,3);
t101 = t102 * t63 - t103 * t62;
t100 = t102 * t62 + t103 * t63;
t99 = t62 / 0.2e1;
t98 = pkin(2) * t66;
t95 = t57 * t62;
t94 = t57 * t63;
t65 = sin(qJ(4));
t67 = cos(qJ(4));
t29 = -Icges(5,3) * t58 + (Icges(5,5) * t67 - Icges(5,6) * t65) * t57;
t93 = t58 * t29;
t92 = t62 * t65;
t91 = t62 * t67;
t90 = t63 * t65;
t89 = t63 * t67;
t88 = t87 * t68 * pkin(2);
t86 = Icges(5,5) * t57;
t85 = Icges(5,6) * t57;
t84 = Icges(5,3) * t57;
t83 = -t57 * rSges(4,1) - t58 * rSges(4,2) - t98;
t32 = -t58 * rSges(5,3) + (rSges(5,1) * t67 - rSges(5,2) * t65) * t57;
t82 = -t57 * pkin(3) + t58 * pkin(5) - t32 - t98;
t81 = pkin(3) * t58 + pkin(5) * t57;
t54 = t66 * rSges(3,1) + t68 * rSges(3,2);
t50 = t58 * t89 + t92;
t49 = -t58 * t90 + t91;
t48 = t58 * t91 - t90;
t47 = -t58 * t92 - t89;
t34 = t83 * t63;
t33 = t83 * t62;
t31 = -Icges(5,5) * t58 + (Icges(5,1) * t67 - Icges(5,4) * t65) * t57;
t30 = -Icges(5,6) * t58 + (Icges(5,4) * t67 - Icges(5,2) * t65) * t57;
t26 = t87 * (rSges(3,1) * t68 - rSges(3,2) * t66);
t25 = t50 * rSges(5,1) + t49 * rSges(5,2) + rSges(5,3) * t94;
t24 = t48 * rSges(5,1) + t47 * rSges(5,2) + rSges(5,3) * t95;
t23 = Icges(5,1) * t50 + Icges(5,4) * t49 + t63 * t86;
t22 = Icges(5,1) * t48 + Icges(5,4) * t47 + t62 * t86;
t21 = Icges(5,4) * t50 + Icges(5,2) * t49 + t63 * t85;
t20 = Icges(5,4) * t48 + Icges(5,2) * t47 + t62 * t85;
t19 = Icges(5,5) * t50 + Icges(5,6) * t49 + t63 * t84;
t18 = Icges(5,5) * t48 + Icges(5,6) * t47 + t62 * t84;
t17 = t82 * t63;
t16 = t82 * t62;
t15 = -t58 * t25 - t32 * t94;
t14 = t58 * t24 + t32 * t95;
t13 = t88 + t87 * (rSges(4,1) * t58 - rSges(4,2) * t57);
t12 = (t24 * t63 - t25 * t62) * t57;
t11 = -t58 * t19 + (-t21 * t65 + t23 * t67) * t57;
t10 = -t58 * t18 + (-t20 * t65 + t22 * t67) * t57;
t9 = t19 * t94 + t49 * t21 + t50 * t23;
t8 = t18 * t94 + t49 * t20 + t50 * t22;
t7 = t19 * t95 + t47 * t21 + t48 * t23;
t6 = t18 * t95 + t47 * t20 + t48 * t22;
t5 = (t81 * t63 + t25) * t63 + (t81 * t62 + t24) * t62 + t88;
t4 = t9 * t62 - t8 * t63;
t3 = -t6 * t63 + t7 * t62;
t2 = -(t49 * t30 + t50 * t31) * t58 + (t8 * t62 + (t9 - t93) * t63) * t57;
t1 = -(t47 * t30 + t48 * t31) * t58 + (t7 * t63 + (t6 - t93) * t62) * t57;
t27 = [m(2) + m(3) + m(4) + m(5); m(3) * t26 + m(4) * t13 + m(5) * t5; m(3) * (t87 * t54 ^ 2 + t26 ^ 2) + m(4) * (t13 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2 + t5 ^ 2) + (t100 * t59 + t4) * t62 + (-t3 + t101 * t60 + (t100 * t63 + t101 * t62) * t62) * t63; 0; m(4) * (-t63 * t33 + t62 * t34) + m(5) * (-t63 * t16 + t62 * t17); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t87; m(5) * t12; -t58 * (-t10 * t63 + t11 * t62) / 0.2e1 + t2 * t99 - t63 * t1 / 0.2e1 + m(5) * (t12 * t5 + t14 * t17 + t15 * t16) + (t63 * t4 / 0.2e1 + t3 * t99) * t57; m(5) * (t14 * t62 - t15 * t63); m(5) * (t12 ^ 2 + t14 ^ 2 + t15 ^ 2) + t2 * t94 + t1 * t95 - t58 * (t58 ^ 2 * t29 + (t11 * t63 + t10 * t62 - (-t30 * t65 + t31 * t67) * t58) * t57);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t27(1), t27(2), t27(4), t27(7); t27(2), t27(3), t27(5), t27(8); t27(4), t27(5), t27(6), t27(9); t27(7), t27(8), t27(9), t27(10);];
Mq = res;
