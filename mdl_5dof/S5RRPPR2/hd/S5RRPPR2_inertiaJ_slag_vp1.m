% Calculate joint inertia matrix for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:03
% EndTime: 2019-12-05 18:20:05
% DurationCPUTime: 0.43s
% Computational Cost: add. (1758->146), mult. (1474->215), div. (0->0), fcn. (1472->10), ass. (0->74)
t64 = qJ(1) + qJ(2);
t62 = sin(t64);
t91 = pkin(2) * t62;
t63 = cos(t64);
t90 = pkin(2) * t63;
t68 = sin(qJ(1));
t89 = t68 * pkin(1);
t70 = cos(qJ(1));
t88 = t70 * pkin(1);
t65 = sin(pkin(9));
t66 = cos(pkin(9));
t67 = sin(qJ(5));
t69 = cos(qJ(5));
t40 = -Icges(6,3) * t66 + (Icges(6,5) * t69 - Icges(6,6) * t67) * t65;
t41 = -Icges(6,6) * t66 + (Icges(6,4) * t69 - Icges(6,2) * t67) * t65;
t42 = -Icges(6,5) * t66 + (Icges(6,1) * t69 - Icges(6,4) * t67) * t65;
t16 = -t66 * t40 + (-t41 * t67 + t42 * t69) * t65;
t87 = t16 * t66;
t61 = pkin(8) + t64;
t59 = sin(t61);
t86 = t59 * t65;
t60 = cos(t61);
t85 = t60 * t65;
t84 = t66 * t67;
t83 = t66 * t69;
t36 = t59 * t84 + t60 * t69;
t37 = -t59 * t83 + t60 * t67;
t82 = t37 * rSges(6,1) + t36 * rSges(6,2);
t81 = Icges(6,5) * t65;
t80 = Icges(6,6) * t65;
t79 = Icges(6,3) * t65;
t76 = -rSges(5,1) * t66 - pkin(3);
t75 = t60 * qJ(4) - t91;
t50 = -t63 * rSges(3,1) + t62 * rSges(3,2);
t49 = -t62 * rSges(3,1) - t63 * rSges(3,2);
t38 = t59 * t69 - t60 * t84;
t39 = t59 * t67 + t60 * t83;
t74 = -t39 * rSges(6,1) - t38 * rSges(6,2);
t17 = Icges(6,5) * t37 + Icges(6,6) * t36 - t59 * t79;
t19 = Icges(6,4) * t37 + Icges(6,2) * t36 - t59 * t80;
t21 = Icges(6,1) * t37 + Icges(6,4) * t36 - t59 * t81;
t3 = -t66 * t17 + (-t19 * t67 + t21 * t69) * t65;
t18 = Icges(6,5) * t39 + Icges(6,6) * t38 + t60 * t79;
t20 = Icges(6,4) * t39 + Icges(6,2) * t38 + t60 * t80;
t22 = Icges(6,1) * t39 + Icges(6,4) * t38 + t60 * t81;
t4 = -t66 * t18 + (-t20 * t67 + t22 * t69) * t65;
t8 = t36 * t41 + t37 * t42 - t40 * t86;
t9 = t38 * t41 + t39 * t42 + t40 * t85;
t73 = -t87 - (t3 + t8) * t86 / 0.2e1 + (t4 + t9) * t85 / 0.2e1;
t35 = -t60 * rSges(4,1) + t59 * rSges(4,2) - t90;
t72 = -pkin(4) * t66 - pkin(3) + (-rSges(6,3) - pkin(7)) * t65;
t34 = -t59 * rSges(4,1) - t60 * rSges(4,2) - t91;
t27 = rSges(5,2) * t86 + t60 * rSges(5,3) + t76 * t59 + t75;
t71 = Icges(5,2) * t66 ^ 2 + Icges(3,3) + Icges(4,3) + t16 + (Icges(5,1) * t65 + 0.2e1 * Icges(5,4) * t66) * t65;
t28 = -t90 + rSges(5,2) * t85 + t76 * t60 + (-rSges(5,3) - qJ(4)) * t59;
t12 = t72 * t59 + t75 + t82;
t13 = -t59 * qJ(4) + t72 * t60 + t74 - t90;
t54 = -t70 * rSges(2,1) + t68 * rSges(2,2);
t53 = -t68 * rSges(2,1) - t70 * rSges(2,2);
t45 = t50 - t88;
t44 = t49 - t89;
t43 = -t66 * rSges(6,3) + (rSges(6,1) * t69 - rSges(6,2) * t67) * t65;
t31 = t35 - t88;
t30 = t34 - t89;
t26 = t28 - t88;
t25 = t27 - t89;
t24 = rSges(6,3) * t85 - t74;
t23 = -rSges(6,3) * t86 + t82;
t15 = t66 * t24 + t43 * t85;
t14 = -t66 * t23 + t43 * t86;
t11 = t13 - t88;
t10 = t12 - t89;
t5 = (-t23 * t60 - t24 * t59) * t65;
t1 = [Icges(2,3) + m(2) * (t53 ^ 2 + t54 ^ 2) + m(3) * (t44 ^ 2 + t45 ^ 2) + m(4) * (t30 ^ 2 + t31 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t71; m(3) * (t49 * t44 + t50 * t45) + m(4) * (t34 * t30 + t35 * t31) + m(5) * (t27 * t25 + t28 * t26) + m(6) * (t12 * t10 + t13 * t11) + t71; m(6) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(3) * (t49 ^ 2 + t50 ^ 2) + t71; 0; 0; m(4) + m(5) + m(6); m(5) * (t59 * t25 + t60 * t26) + m(6) * (t59 * t10 + t60 * t11); m(6) * (t59 * t12 + t60 * t13) + m(5) * (t59 * t27 + t60 * t28); 0; 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t59 ^ 2 + t60 ^ 2); m(6) * (t14 * t10 + t15 * t11) + t73; m(6) * (t14 * t12 + t15 * t13) + t73; m(6) * t5; m(6) * (t14 * t59 + t15 * t60); m(6) * (t14 ^ 2 + t15 ^ 2 + t5 ^ 2) - t66 * (-t87 + (-t3 * t59 + t4 * t60) * t65) - (-t8 * t66 - (-t17 * t86 + t36 * t19 + t37 * t21) * t86 + (-t18 * t86 + t36 * t20 + t37 * t22) * t85) * t86 + (-t9 * t66 - (t17 * t85 + t38 * t19 + t39 * t21) * t86 + (t18 * t85 + t38 * t20 + t39 * t22) * t85) * t85;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
