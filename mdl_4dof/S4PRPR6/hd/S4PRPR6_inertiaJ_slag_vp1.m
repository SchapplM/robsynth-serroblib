% Calculate joint inertia matrix for
% S4PRPR6
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR6_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR6_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:20
% EndTime: 2019-12-31 16:24:22
% DurationCPUTime: 0.69s
% Computational Cost: add. (1304->156), mult. (2006->261), div. (0->0), fcn. (2156->8), ass. (0->85)
t69 = sin(pkin(6));
t64 = t69 ^ 2;
t71 = cos(pkin(6));
t65 = t71 ^ 2;
t90 = t64 + t65;
t74 = cos(qJ(2));
t100 = t74 ^ 2;
t99 = t69 / 0.2e1;
t98 = -m(4) - m(5);
t73 = sin(qJ(2));
t97 = t69 * t73;
t96 = t69 * t74;
t95 = t71 * t73;
t94 = t71 * t74;
t66 = pkin(7) + qJ(4);
t62 = sin(t66);
t63 = cos(t66);
t36 = -Icges(5,3) * t74 + (Icges(5,5) * t63 - Icges(5,6) * t62) * t73;
t93 = t74 * t36;
t59 = t73 * pkin(2) - t74 * qJ(3);
t68 = sin(pkin(7));
t70 = cos(pkin(7));
t92 = t74 * rSges(4,3) - (rSges(4,1) * t70 - rSges(4,2) * t68) * t73 - t59;
t81 = pkin(2) * t74 + qJ(3) * t73;
t91 = t90 * t81;
t89 = Icges(4,5) * t73;
t88 = Icges(5,5) * t73;
t87 = Icges(4,6) * t73;
t86 = Icges(5,6) * t73;
t85 = Icges(4,3) * t73;
t84 = Icges(5,3) * t73;
t39 = -t74 * rSges(5,3) + (rSges(5,1) * t63 - rSges(5,2) * t62) * t73;
t61 = t70 * pkin(3) + pkin(2);
t72 = -pkin(5) - qJ(3);
t83 = -(qJ(3) + t72) * t74 - (-pkin(2) + t61) * t73 - t39 - t59;
t76 = Icges(3,5) * t74 - Icges(3,6) * t73;
t75 = t61 * t74 - t72 * t73 - t81;
t60 = t73 * rSges(3,1) + t74 * rSges(3,2);
t57 = t69 * t68 + t70 * t94;
t56 = -t68 * t94 + t69 * t70;
t55 = -t71 * t68 + t70 * t96;
t54 = -t68 * t96 - t71 * t70;
t49 = t69 * t62 + t63 * t94;
t48 = -t62 * t94 + t69 * t63;
t47 = -t71 * t62 + t63 * t96;
t46 = -t62 * t96 - t71 * t63;
t41 = Icges(3,3) * t69 + t76 * t71;
t40 = -Icges(3,3) * t71 + t76 * t69;
t38 = -Icges(5,5) * t74 + (Icges(5,1) * t63 - Icges(5,4) * t62) * t73;
t37 = -Icges(5,6) * t74 + (Icges(5,4) * t63 - Icges(5,2) * t62) * t73;
t34 = t92 * t71;
t33 = t92 * t69;
t32 = Icges(4,1) * t57 + Icges(4,4) * t56 + t71 * t89;
t31 = Icges(4,1) * t55 + Icges(4,4) * t54 + t69 * t89;
t30 = Icges(4,4) * t57 + Icges(4,2) * t56 + t71 * t87;
t29 = Icges(4,4) * t55 + Icges(4,2) * t54 + t69 * t87;
t28 = Icges(4,5) * t57 + Icges(4,6) * t56 + t71 * t85;
t27 = Icges(4,5) * t55 + Icges(4,6) * t54 + t69 * t85;
t26 = t90 * (rSges(3,1) * t74 - rSges(3,2) * t73);
t25 = t49 * rSges(5,1) + t48 * rSges(5,2) + rSges(5,3) * t95;
t24 = t47 * rSges(5,1) + t46 * rSges(5,2) + rSges(5,3) * t97;
t23 = Icges(5,1) * t49 + Icges(5,4) * t48 + t71 * t88;
t22 = Icges(5,1) * t47 + Icges(5,4) * t46 + t69 * t88;
t21 = Icges(5,4) * t49 + Icges(5,2) * t48 + t71 * t86;
t20 = Icges(5,4) * t47 + Icges(5,2) * t46 + t69 * t86;
t19 = Icges(5,5) * t49 + Icges(5,6) * t48 + t71 * t84;
t18 = Icges(5,5) * t47 + Icges(5,6) * t46 + t69 * t84;
t17 = t83 * t71;
t16 = t83 * t69;
t15 = -t74 * t25 - t39 * t95;
t14 = t74 * t24 + t39 * t97;
t13 = (t24 * t71 - t25 * t69) * t73;
t12 = t69 * (t55 * rSges(4,1) + t54 * rSges(4,2) + rSges(4,3) * t97) + t71 * (t57 * rSges(4,1) + t56 * rSges(4,2) + rSges(4,3) * t95) + t91;
t11 = -t74 * t19 + (-t21 * t62 + t23 * t63) * t73;
t10 = -t74 * t18 + (-t20 * t62 + t22 * t63) * t73;
t9 = t19 * t95 + t48 * t21 + t49 * t23;
t8 = t18 * t95 + t48 * t20 + t49 * t22;
t7 = t19 * t97 + t46 * t21 + t47 * t23;
t6 = t18 * t97 + t46 * t20 + t47 * t22;
t5 = (t75 * t71 + t25) * t71 + (t75 * t69 + t24) * t69 + t91;
t4 = t9 * t69 - t8 * t71;
t3 = -t6 * t71 + t7 * t69;
t2 = -(t48 * t37 + t49 * t38) * t74 + (t8 * t69 + (t9 - t93) * t71) * t73;
t1 = -(t46 * t37 + t47 * t38) * t74 + (t7 * t71 + (t6 - t93) * t69) * t73;
t35 = [m(2) + m(3) - t98; m(3) * t26 + m(4) * t12 + m(5) * t5; m(3) * (t90 * t60 ^ 2 + t26 ^ 2) + m(4) * (t12 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2 + t5 ^ 2) + (-t65 * t40 - t3 + (t27 * t97 + t54 * t29 + t55 * t31) * t71) * t71 + (t64 * t41 + t4 + (t28 * t95 + t56 * t30 + t57 * t32) * t69 + (-t27 * t95 - t28 * t97 - t56 * t29 - t54 * t30 - t57 * t31 - t55 * t32 - t69 * t40 + t71 * t41) * t71) * t69; t98 * t74; m(4) * (-t74 * t12 + (t33 * t69 + t34 * t71) * t73) + m(5) * (-t74 * t5 + (t16 * t69 + t17 * t71) * t73); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t90 * t73 ^ 2 + t100); m(5) * t13; -t74 * (-t10 * t71 + t11 * t69) / 0.2e1 + t2 * t99 - t71 * t1 / 0.2e1 + m(5) * (t13 * t5 + t14 * t17 + t15 * t16) + (t71 * t4 / 0.2e1 + t3 * t99) * t73; m(5) * (-t13 * t74 + (t14 * t71 + t15 * t69) * t73); m(5) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) + t2 * t95 + t1 * t97 - t74 * (t100 * t36 + (t11 * t71 + t10 * t69 - (-t37 * t62 + t38 * t63) * t74) * t73);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t35(1), t35(2), t35(4), t35(7); t35(2), t35(3), t35(5), t35(8); t35(4), t35(5), t35(6), t35(9); t35(7), t35(8), t35(9), t35(10);];
Mq = res;
