% Calculate joint inertia matrix for
% S5PRRPR2
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:19
% EndTime: 2019-12-05 16:17:20
% DurationCPUTime: 0.40s
% Computational Cost: add. (1673->134), mult. (1400->203), div. (0->0), fcn. (1416->8), ass. (0->65)
t66 = pkin(8) + qJ(2);
t63 = sin(t66);
t86 = pkin(2) * t63;
t67 = sin(pkin(9));
t68 = cos(pkin(9));
t69 = sin(qJ(5));
t70 = cos(qJ(5));
t38 = -Icges(6,3) * t68 + (Icges(6,5) * t70 - Icges(6,6) * t69) * t67;
t39 = -Icges(6,6) * t68 + (Icges(6,4) * t70 - Icges(6,2) * t69) * t67;
t40 = -Icges(6,5) * t68 + (Icges(6,1) * t70 - Icges(6,4) * t69) * t67;
t16 = -t68 * t38 + (-t39 * t69 + t40 * t70) * t67;
t85 = t16 * t68;
t65 = qJ(3) + t66;
t61 = sin(t65);
t84 = t61 * t67;
t62 = cos(t65);
t83 = t62 * t67;
t82 = t62 * t68;
t81 = t68 * t69;
t80 = t68 * t70;
t79 = t62 * pkin(3) + t61 * qJ(4);
t78 = Icges(6,5) * t67;
t77 = Icges(6,6) * t67;
t76 = Icges(6,3) * t67;
t36 = t61 * t70 - t62 * t81;
t37 = t61 * t69 + t62 * t80;
t24 = t37 * rSges(6,1) + t36 * rSges(6,2) + rSges(6,3) * t83;
t43 = t62 * rSges(4,1) - t61 * rSges(4,2);
t42 = -t61 * rSges(4,1) - t62 * rSges(4,2);
t34 = -t61 * t81 - t62 * t70;
t35 = t61 * t80 - t62 * t69;
t73 = -t35 * rSges(6,1) - t34 * rSges(6,2);
t17 = Icges(6,5) * t35 + Icges(6,6) * t34 + t61 * t76;
t19 = Icges(6,4) * t35 + Icges(6,2) * t34 + t61 * t77;
t21 = Icges(6,1) * t35 + Icges(6,4) * t34 + t61 * t78;
t3 = -t68 * t17 + (-t19 * t69 + t21 * t70) * t67;
t18 = Icges(6,5) * t37 + Icges(6,6) * t36 + t62 * t76;
t20 = Icges(6,4) * t37 + Icges(6,2) * t36 + t62 * t77;
t22 = Icges(6,1) * t37 + Icges(6,4) * t36 + t62 * t78;
t4 = -t68 * t18 + (-t20 * t69 + t22 * t70) * t67;
t8 = t34 * t39 + t35 * t40 + t38 * t84;
t9 = t36 * t39 + t37 * t40 + t38 * t83;
t72 = -t85 + (t3 + t8) * t84 / 0.2e1 + (t4 + t9) * t83 / 0.2e1;
t13 = pkin(4) * t82 + pkin(7) * t83 + t24 + t79;
t28 = rSges(5,1) * t82 - rSges(5,2) * t83 + t61 * rSges(5,3) + t79;
t55 = t62 * qJ(4);
t27 = rSges(5,2) * t84 + t55 + t62 * rSges(5,3) + (-rSges(5,1) * t68 - pkin(3)) * t61;
t71 = Icges(5,2) * t68 ^ 2 + Icges(4,3) + t16 + (Icges(5,1) * t67 + 0.2e1 * Icges(5,4) * t68) * t67;
t12 = t55 + (-pkin(4) * t68 - pkin(3) + (-rSges(6,3) - pkin(7)) * t67) * t61 + t73;
t64 = cos(t66);
t60 = pkin(2) * t64;
t48 = t64 * rSges(3,1) - t63 * rSges(3,2);
t47 = -t63 * rSges(3,1) - t64 * rSges(3,2);
t41 = -t68 * rSges(6,3) + (rSges(6,1) * t70 - rSges(6,2) * t69) * t67;
t33 = t43 + t60;
t32 = t42 - t86;
t26 = t28 + t60;
t25 = t27 - t86;
t23 = rSges(6,3) * t84 - t73;
t15 = -t68 * t24 - t41 * t83;
t14 = t68 * t23 + t41 * t84;
t11 = t60 + t13;
t10 = t12 - t86;
t5 = (t23 * t62 - t24 * t61) * t67;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; Icges(3,3) + m(6) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2) + m(4) * (t32 ^ 2 + t33 ^ 2) + m(3) * (t47 ^ 2 + t48 ^ 2) + t71; 0; m(6) * (t12 * t10 + t13 * t11) + m(5) * (t27 * t25 + t28 * t26) + m(4) * (t42 * t32 + t43 * t33) + t71; m(6) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(4) * (t42 ^ 2 + t43 ^ 2) + t71; 0; m(6) * (t61 * t10 - t62 * t11) + m(5) * (t61 * t25 - t62 * t26); m(6) * (t61 * t12 - t62 * t13) + m(5) * (t61 * t27 - t62 * t28); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t61 ^ 2 + t62 ^ 2); m(6) * t5; m(6) * (t14 * t10 + t15 * t11) + t72; m(6) * (t14 * t12 + t15 * t13) + t72; m(6) * (t14 * t61 - t15 * t62); m(6) * (t14 ^ 2 + t15 ^ 2 + t5 ^ 2) + ((t18 * t83 + t36 * t20 + t37 * t22) * t83 + (t17 * t83 + t36 * t19 + t37 * t21) * t84 - t9 * t68) * t83 + ((t18 * t84 + t34 * t20 + t35 * t22) * t83 + (t17 * t84 + t34 * t19 + t35 * t21) * t84 - t8 * t68) * t84 - t68 * (-t85 + (t3 * t61 + t4 * t62) * t67);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
