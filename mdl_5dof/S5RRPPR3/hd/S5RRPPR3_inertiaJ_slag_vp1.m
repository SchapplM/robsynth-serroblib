% Calculate joint inertia matrix for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:28
% EndTime: 2019-12-31 19:26:29
% DurationCPUTime: 0.38s
% Computational Cost: add. (995->93), mult. (676->127), div. (0->0), fcn. (534->8), ass. (0->54)
t60 = qJ(1) + qJ(2);
t56 = pkin(8) + t60;
t53 = sin(t56);
t50 = t53 ^ 2;
t54 = cos(t56);
t51 = t54 ^ 2;
t63 = cos(qJ(5));
t91 = Icges(6,5) * t63;
t61 = sin(qJ(5));
t90 = Icges(6,6) * t61;
t35 = -t90 + t91;
t89 = rSges(6,1) * t61 + rSges(6,2) * t63;
t88 = t53 * t54;
t57 = sin(t60);
t85 = pkin(2) * t57;
t62 = sin(qJ(1));
t84 = t62 * pkin(1);
t81 = t89 * t54;
t80 = t50 + t51;
t77 = t54 * rSges(6,3) + t89 * t53;
t58 = cos(t60);
t55 = pkin(2) * t58;
t76 = t54 * pkin(3) + t53 * qJ(4) + t55;
t75 = t35 * t51 + (t91 / 0.2e1 - t90 / 0.2e1 + t35 / 0.2e1) * t50;
t74 = t54 * qJ(4) - t85;
t28 = t58 * rSges(3,1) - t57 * rSges(3,2);
t23 = t54 * rSges(4,1) - t53 * rSges(4,2) + t55;
t27 = -t57 * rSges(3,1) - t58 * rSges(3,2);
t68 = Icges(6,5) * t61 + Icges(6,6) * t63;
t8 = (-rSges(6,3) - pkin(3) - pkin(7)) * t53 + t74 + t81;
t6 = t8 - t84;
t64 = cos(qJ(1));
t59 = t64 * pkin(1);
t9 = t54 * pkin(7) + t76 + t77;
t7 = t59 + t9;
t67 = m(6) * (t53 * t6 - t54 * t7);
t66 = m(6) * (t53 * t8 - t54 * t9);
t13 = -t54 * rSges(5,2) + t53 * rSges(5,3) + t76;
t65 = Icges(6,1) * t63 ^ 2 + Icges(5,1) + Icges(3,3) + Icges(4,3) + (-0.2e1 * Icges(6,4) * t63 + Icges(6,2) * t61) * t61;
t22 = -t53 * rSges(4,1) - t54 * rSges(4,2) - t85;
t12 = t54 * rSges(5,3) + t74 + (rSges(5,2) - pkin(3)) * t53;
t40 = t64 * rSges(2,1) - t62 * rSges(2,2);
t39 = t63 * rSges(6,1) - t61 * rSges(6,2);
t38 = -t62 * rSges(2,1) - t64 * rSges(2,2);
t25 = t28 + t59;
t24 = t27 - t84;
t21 = t23 + t59;
t20 = t22 - t84;
t15 = Icges(6,3) * t53 - t68 * t54;
t14 = Icges(6,3) * t54 + t68 * t53;
t11 = t13 + t59;
t10 = t12 - t84;
t3 = t54 * (t53 * rSges(6,3) - t81) - t53 * t77;
t1 = [Icges(2,3) + m(6) * (t6 ^ 2 + t7 ^ 2) + m(5) * (t10 ^ 2 + t11 ^ 2) + m(3) * (t24 ^ 2 + t25 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(2) * (t38 ^ 2 + t40 ^ 2) + t65; m(6) * (t8 * t6 + t9 * t7) + m(5) * (t12 * t10 + t13 * t11) + m(3) * (t27 * t24 + t28 * t25) + m(4) * (t22 * t20 + t23 * t21) + t65; m(5) * (t12 ^ 2 + t13 ^ 2) + m(6) * (t8 ^ 2 + t9 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(3) * (t27 ^ 2 + t28 ^ 2) + t65; 0; 0; m(4) + m(5) + m(6); t67 + m(5) * (t53 * t10 - t54 * t11); m(5) * (t53 * t12 - t54 * t13) + t66; 0; 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t80; t39 * t67 + t75; t39 * t66 + t75; m(6) * t3; m(6) * t80 * t39; m(6) * (t80 * t39 ^ 2 + t3 ^ 2) + t54 * (t51 * t14 + t15 * t88) + t53 * (t14 * t88 + t50 * t15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
