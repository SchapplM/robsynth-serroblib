% Calculate joint inertia matrix for
% S5PRRPR1
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:34
% EndTime: 2019-12-05 16:15:35
% DurationCPUTime: 0.38s
% Computational Cost: add. (1090->86), mult. (644->124), div. (0->0), fcn. (526->8), ass. (0->51)
t61 = pkin(8) + qJ(2);
t59 = qJ(3) + t61;
t52 = sin(t59);
t49 = t52 ^ 2;
t53 = cos(t59);
t50 = t53 ^ 2;
t60 = pkin(9) + qJ(5);
t55 = sin(t60);
t90 = Icges(6,5) * t55;
t57 = cos(t60);
t89 = Icges(6,6) * t57;
t27 = t89 + t90;
t88 = qJ(4) + rSges(5,3);
t62 = sin(pkin(9));
t63 = cos(pkin(9));
t87 = -rSges(5,1) * t63 + rSges(5,2) * t62 - pkin(3);
t86 = t52 * t53;
t33 = t55 * rSges(6,1) + t57 * rSges(6,2);
t83 = m(6) * t33;
t56 = sin(t61);
t82 = pkin(2) * t56;
t80 = rSges(6,1) * t57;
t78 = rSges(6,2) * t55;
t77 = t53 * rSges(6,3) + t52 * t78;
t76 = t49 + t50;
t73 = t27 * t50 + (t90 / 0.2e1 + t89 / 0.2e1 + t27 / 0.2e1) * t49;
t25 = t53 * rSges(4,1) - t52 * rSges(4,2);
t72 = Icges(5,2) * t63 ^ 2 + Icges(6,2) * t57 ^ 2 + Icges(4,3) + (Icges(5,1) * t62 + 0.2e1 * Icges(5,4) * t63) * t62 + (Icges(6,1) * t55 + 0.2e1 * Icges(6,4) * t57) * t55;
t24 = -t52 * rSges(4,1) - t53 * rSges(4,2);
t66 = Icges(6,5) * t57 - Icges(6,6) * t55;
t65 = t52 * rSges(6,3) + (-t78 + t80) * t53;
t13 = t88 * t52 - t87 * t53;
t12 = t87 * t52 + t53 * t88;
t54 = t63 * pkin(4) + pkin(3);
t64 = -pkin(7) - qJ(4);
t9 = -t52 * t64 + t53 * t54 + t65;
t8 = -t53 * t64 + (-t54 - t80) * t52 + t77;
t58 = cos(t61);
t51 = pkin(2) * t58;
t35 = t58 * rSges(3,1) - t56 * rSges(3,2);
t34 = -t56 * rSges(3,1) - t58 * rSges(3,2);
t21 = t25 + t51;
t20 = t24 - t82;
t15 = Icges(6,3) * t52 + t66 * t53;
t14 = -Icges(6,3) * t53 + t66 * t52;
t11 = t13 + t51;
t10 = t12 - t82;
t7 = t51 + t9;
t6 = t8 - t82;
t3 = t52 * (t52 * t80 - t77) + t53 * t65;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; Icges(3,3) + m(5) * (t10 ^ 2 + t11 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(3) * (t34 ^ 2 + t35 ^ 2) + t72; 0; m(5) * (t12 * t10 + t13 * t11) + m(6) * (t8 * t6 + t9 * t7) + m(4) * (t24 * t20 + t25 * t21) + t72; m(6) * (t8 ^ 2 + t9 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + t72; 0; m(5) * (t52 * t10 - t53 * t11) + m(6) * (t52 * t6 - t53 * t7); m(6) * (t52 * t8 - t53 * t9) + m(5) * (t52 * t12 - t53 * t13); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t76; m(6) * t3; (-t52 * t7 - t53 * t6) * t83 + t73; (-t52 * t9 - t53 * t8) * t83 + t73; 0; m(6) * (t76 * t33 ^ 2 + t3 ^ 2) + t52 * (-t14 * t86 + t49 * t15) - t53 * (t50 * t14 - t15 * t86);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
