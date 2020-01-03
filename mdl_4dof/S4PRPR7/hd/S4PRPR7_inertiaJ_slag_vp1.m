% Calculate joint inertia matrix for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR7_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR7_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:37
% EndTime: 2019-12-31 16:25:38
% DurationCPUTime: 0.59s
% Computational Cost: add. (714->117), mult. (1796->207), div. (0->0), fcn. (1902->6), ass. (0->66)
t62 = sin(qJ(2));
t64 = cos(qJ(2));
t96 = (-Icges(4,4) + Icges(3,5)) * t64 + (Icges(4,5) - Icges(3,6)) * t62;
t95 = Icges(4,1) + Icges(3,3);
t59 = sin(pkin(6));
t56 = t59 ^ 2;
t60 = cos(pkin(6));
t57 = t60 ^ 2;
t82 = t56 + t57;
t94 = -t96 * t59 + t95 * t60;
t93 = t95 * t59 + t96 * t60;
t92 = t62 ^ 2;
t91 = t59 / 0.2e1;
t90 = -m(4) - m(5);
t89 = t59 * t64;
t88 = t60 * t64;
t61 = sin(qJ(4));
t87 = t61 * t62;
t63 = cos(qJ(4));
t41 = Icges(5,3) * t62 + (-Icges(5,5) * t61 - Icges(5,6) * t63) * t64;
t86 = t62 * t41;
t85 = t62 * t63;
t84 = t82 * (pkin(2) * t64 + qJ(3) * t62);
t53 = t62 * pkin(2) - t64 * qJ(3);
t83 = t62 * rSges(4,2) + t64 * rSges(4,3) - t53;
t81 = Icges(5,5) * t64;
t80 = Icges(5,6) * t64;
t79 = Icges(5,3) * t64;
t44 = t62 * rSges(5,3) + (-rSges(5,1) * t61 - rSges(5,2) * t63) * t64;
t78 = -pkin(5) * t62 - t44 - t53;
t55 = t62 * rSges(3,1) + t64 * rSges(3,2);
t51 = t59 * t87 - t60 * t63;
t50 = t59 * t85 + t60 * t61;
t49 = t59 * t63 + t60 * t87;
t48 = -t59 * t61 + t60 * t85;
t43 = Icges(5,5) * t62 + (-Icges(5,1) * t61 - Icges(5,4) * t63) * t64;
t42 = Icges(5,6) * t62 + (-Icges(5,4) * t61 - Icges(5,2) * t63) * t64;
t28 = t83 * t60;
t27 = t83 * t59;
t26 = t51 * rSges(5,1) + t50 * rSges(5,2) + rSges(5,3) * t89;
t25 = t49 * rSges(5,1) + t48 * rSges(5,2) + rSges(5,3) * t88;
t24 = t78 * t60;
t23 = t78 * t59;
t22 = Icges(5,1) * t51 + Icges(5,4) * t50 + t59 * t81;
t21 = Icges(5,1) * t49 + Icges(5,4) * t48 + t60 * t81;
t20 = Icges(5,4) * t51 + Icges(5,2) * t50 + t59 * t80;
t19 = Icges(5,4) * t49 + Icges(5,2) * t48 + t60 * t80;
t18 = Icges(5,5) * t51 + Icges(5,6) * t50 + t59 * t79;
t17 = Icges(5,5) * t49 + Icges(5,6) * t48 + t60 * t79;
t16 = t82 * (rSges(3,1) * t64 - rSges(3,2) * t62);
t15 = t62 * t25 - t44 * t88;
t14 = -t62 * t26 + t44 * t89;
t13 = t84 + t82 * (-rSges(4,2) * t64 + rSges(4,3) * t62);
t12 = (-t25 * t59 + t26 * t60) * t64;
t11 = t62 * t18 + (-t20 * t63 - t22 * t61) * t64;
t10 = t62 * t17 + (-t19 * t63 - t21 * t61) * t64;
t9 = t82 * t64 * pkin(5) + t60 * t25 + t59 * t26 + t84;
t8 = t18 * t89 + t50 * t20 + t51 * t22;
t7 = t17 * t89 + t50 * t19 + t51 * t21;
t6 = t18 * t88 + t48 * t20 + t49 * t22;
t5 = t17 * t88 + t48 * t19 + t49 * t21;
t4 = t7 * t59 - t8 * t60;
t3 = t5 * t59 - t6 * t60;
t2 = (t50 * t42 + t51 * t43) * t62 + (t7 * t60 + (t8 + t86) * t59) * t64;
t1 = (t48 * t42 + t49 * t43) * t62 + (t6 * t59 + (t5 + t86) * t60) * t64;
t29 = [m(2) + m(3) - t90; m(3) * t16 + m(4) * t13 + m(5) * t9; m(3) * (t82 * t55 ^ 2 + t16 ^ 2) + m(4) * (t13 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(5) * (t23 ^ 2 + t24 ^ 2 + t9 ^ 2) + (t94 * t57 - t4) * t60 + (t3 + t93 * t56 + (t94 * t59 + t93 * t60) * t60) * t59; t90 * t64; m(4) * (-t64 * t13 + (t27 * t59 + t28 * t60) * t62) + m(5) * (-t64 * t9 + (t23 * t59 + t24 * t60) * t62); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t64 ^ 2 + t82 * t92); m(5) * t12; t62 * (t10 * t59 - t11 * t60) / 0.2e1 + t1 * t91 - t60 * t2 / 0.2e1 + m(5) * (t12 * t9 + t14 * t24 + t15 * t23) + (t60 * t3 / 0.2e1 + t4 * t91) * t64; m(5) * (-t12 * t64 + (t14 * t60 + t15 * t59) * t62); m(5) * (t12 ^ 2 + t14 ^ 2 + t15 ^ 2) + t1 * t88 + t2 * t89 + t62 * (t92 * t41 + (t10 * t60 + t11 * t59 + (-t42 * t63 - t43 * t61) * t62) * t64);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t29(1), t29(2), t29(4), t29(7); t29(2), t29(3), t29(5), t29(8); t29(4), t29(5), t29(6), t29(9); t29(7), t29(8), t29(9), t29(10);];
Mq = res;
