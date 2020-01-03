% Calculate joint inertia matrix for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:03
% EndTime: 2019-12-31 17:38:05
% DurationCPUTime: 0.96s
% Computational Cost: add. (1441->167), mult. (3721->274), div. (0->0), fcn. (4334->8), ass. (0->77)
t82 = sin(pkin(7));
t78 = t82 ^ 2;
t83 = cos(pkin(7));
t79 = t83 ^ 2;
t109 = t78 + t79;
t85 = sin(qJ(2));
t87 = cos(qJ(2));
t117 = (Icges(4,4) + Icges(3,5)) * t87 + (-Icges(3,6) + Icges(4,6)) * t85;
t116 = Icges(4,2) + Icges(3,3);
t115 = t116 * t83 - t117 * t82;
t114 = t116 * t82 + t117 * t83;
t108 = cos(pkin(8));
t81 = sin(pkin(8));
t70 = t85 * t108 - t87 * t81;
t111 = t109 * (pkin(2) * t87 + qJ(3) * t85);
t72 = t85 * pkin(2) - t87 * qJ(3);
t110 = -t85 * rSges(4,1) + t87 * rSges(4,3) - t72;
t107 = -m(4) - m(5) - m(6);
t106 = m(5) / 0.2e1 + m(6) / 0.2e1;
t105 = -pkin(3) * t85 - t72;
t103 = t109 * pkin(3) * t87 + t111;
t69 = t87 * t108 + t85 * t81;
t102 = -t70 * rSges(5,1) + t69 * rSges(5,2) + t105;
t84 = sin(qJ(5));
t86 = cos(qJ(5));
t32 = t69 * rSges(6,3) + (rSges(6,1) * t86 - rSges(6,2) * t84) * t70;
t94 = -t70 * pkin(4) - t69 * pkin(6) + t105 - t32;
t74 = t85 * rSges(3,1) + t87 * rSges(3,2);
t63 = t69 * t83;
t62 = t70 * t83;
t61 = t69 * t82;
t60 = t70 * t82;
t47 = t110 * t83;
t46 = t110 * t82;
t45 = t63 * t86 - t82 * t84;
t44 = -t63 * t84 - t82 * t86;
t43 = t61 * t86 + t83 * t84;
t42 = -t61 * t84 + t83 * t86;
t39 = Icges(5,1) * t63 + Icges(5,4) * t62 - Icges(5,5) * t82;
t38 = Icges(5,1) * t61 + Icges(5,4) * t60 + Icges(5,5) * t83;
t37 = Icges(5,4) * t63 + Icges(5,2) * t62 - Icges(5,6) * t82;
t36 = Icges(5,4) * t61 + Icges(5,2) * t60 + Icges(5,6) * t83;
t35 = Icges(5,5) * t63 + Icges(5,6) * t62 - Icges(5,3) * t82;
t34 = Icges(5,5) * t61 + Icges(5,6) * t60 + Icges(5,3) * t83;
t33 = t109 * (rSges(3,1) * t87 - rSges(3,2) * t85);
t31 = Icges(6,5) * t69 + (Icges(6,1) * t86 - Icges(6,4) * t84) * t70;
t30 = Icges(6,6) * t69 + (Icges(6,4) * t86 - Icges(6,2) * t84) * t70;
t29 = Icges(6,3) * t69 + (Icges(6,5) * t86 - Icges(6,6) * t84) * t70;
t28 = t102 * t83;
t27 = t102 * t82;
t26 = t45 * rSges(6,1) + t44 * rSges(6,2) - t62 * rSges(6,3);
t25 = t43 * rSges(6,1) + t42 * rSges(6,2) - t60 * rSges(6,3);
t24 = Icges(6,1) * t45 + Icges(6,4) * t44 - Icges(6,5) * t62;
t23 = Icges(6,1) * t43 + Icges(6,4) * t42 - Icges(6,5) * t60;
t22 = Icges(6,4) * t45 + Icges(6,2) * t44 - Icges(6,6) * t62;
t21 = Icges(6,4) * t43 + Icges(6,2) * t42 - Icges(6,6) * t60;
t20 = Icges(6,5) * t45 + Icges(6,6) * t44 - Icges(6,3) * t62;
t19 = Icges(6,5) * t43 + Icges(6,6) * t42 - Icges(6,3) * t60;
t18 = t111 + t109 * (rSges(4,1) * t87 + rSges(4,3) * t85);
t17 = t94 * t83;
t16 = t94 * t82;
t15 = t69 * t26 + t62 * t32;
t14 = -t69 * t25 - t60 * t32;
t13 = t82 * (t61 * rSges(5,1) + t60 * rSges(5,2)) + t83 * (t63 * rSges(5,1) + t62 * rSges(5,2)) + t103;
t12 = -t62 * t25 + t60 * t26;
t11 = t69 * t20 + (-t22 * t84 + t24 * t86) * t70;
t10 = t69 * t19 + (-t21 * t84 + t23 * t86) * t70;
t9 = -t62 * t20 + t44 * t22 + t45 * t24;
t8 = -t62 * t19 + t44 * t21 + t45 * t23;
t7 = -t60 * t20 + t42 * t22 + t43 * t24;
t6 = -t60 * t19 + t42 * t21 + t43 * t23;
t5 = (t63 * pkin(4) - t62 * pkin(6) + t26) * t83 + (t61 * pkin(4) - t60 * pkin(6) + t25) * t82 + t103;
t4 = -t8 * t83 + t9 * t82;
t3 = -t6 * t83 + t7 * t82;
t2 = -t9 * t62 - t8 * t60 + (-t62 * t29 + t44 * t30 + t45 * t31) * t69;
t1 = -t7 * t62 - t6 * t60 + (-t60 * t29 + t42 * t30 + t43 * t31) * t69;
t40 = [m(2) + m(3) - t107; m(3) * t33 + m(4) * t18 + m(5) * t13 + m(6) * t5; m(6) * (t16 ^ 2 + t17 ^ 2 + t5 ^ 2) + m(5) * (t13 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(4) * (t18 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(3) * (t109 * t74 ^ 2 + t33 ^ 2) + (-t3 + (t83 * t34 + t60 * t36 + t61 * t38) * t83 + t115 * t79) * t83 + (t4 + (-t82 * t35 + t62 * t37 + t63 * t39) * t82 + t114 * t78 + (-t62 * t36 - t63 * t38 - t60 * t37 - t61 * t39 + (t34 + t115) * t82 + (-t35 + t114) * t83) * t83) * t82; t107 * t87; m(6) * (-t87 * t5 + (t16 * t82 + t17 * t83) * t85) + m(5) * (-t87 * t13 + (t27 * t82 + t28 * t83) * t85) + m(4) * (-t87 * t18 + (t46 * t82 + t47 * t83) * t85); 0.2e1 * (m(4) / 0.2e1 + t106) * (t109 * t85 ^ 2 + t87 ^ 2); 0; m(6) * (t83 * t16 - t82 * t17) + m(5) * (t83 * t27 - t82 * t28); 0; 0.2e1 * t106 * t109; m(6) * t12; t82 * t2 / 0.2e1 - t83 * t1 / 0.2e1 - t62 * t4 / 0.2e1 - t60 * t3 / 0.2e1 + t69 * (-t10 * t83 + t11 * t82) / 0.2e1 + m(6) * (t12 * t5 + t14 * t17 + t15 * t16); m(6) * (-t12 * t87 + (t14 * t83 + t15 * t82) * t85); m(6) * (-t14 * t82 + t15 * t83); m(6) * (t12 ^ 2 + t14 ^ 2 + t15 ^ 2) - t62 * t2 - t60 * t1 + t69 * (-t11 * t62 - t10 * t60 + (t69 * t29 + (-t30 * t84 + t31 * t86) * t70) * t69);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t40(1), t40(2), t40(4), t40(7), t40(11); t40(2), t40(3), t40(5), t40(8), t40(12); t40(4), t40(5), t40(6), t40(9), t40(13); t40(7), t40(8), t40(9), t40(10), t40(14); t40(11), t40(12), t40(13), t40(14), t40(15);];
Mq = res;
