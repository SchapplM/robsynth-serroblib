% Calculate joint inertia matrix for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:02
% EndTime: 2020-01-03 11:52:04
% DurationCPUTime: 0.39s
% Computational Cost: add. (1399->97), mult. (760->130), div. (0->0), fcn. (614->10), ass. (0->58)
t65 = qJ(1) + pkin(9);
t62 = qJ(3) + t65;
t59 = qJ(4) + t62;
t53 = sin(t59);
t91 = t53 ^ 2;
t54 = cos(t59);
t90 = t54 ^ 2;
t66 = sin(qJ(5));
t95 = Icges(6,5) * t66;
t68 = cos(qJ(5));
t94 = Icges(6,6) * t68;
t35 = t94 + t95;
t93 = -rSges(6,1) * t68 + rSges(6,2) * t66;
t92 = t53 * t54;
t38 = t66 * rSges(6,1) + t68 * rSges(6,2);
t87 = m(6) * t38;
t24 = t53 * rSges(5,1) + t54 * rSges(5,2);
t57 = sin(t62);
t58 = cos(t62);
t28 = t57 * rSges(4,1) + t58 * rSges(4,2);
t60 = sin(t65);
t67 = sin(qJ(1));
t63 = t67 * pkin(1);
t84 = pkin(2) * t60 + t63;
t61 = cos(t65);
t69 = cos(qJ(1));
t64 = t69 * pkin(1);
t83 = pkin(2) * t61 + t64;
t51 = pkin(3) * t57;
t22 = t51 + t24;
t80 = Icges(6,2) * t68 ^ 2 + Icges(5,3) + (Icges(6,1) * t66 + 0.2e1 * Icges(6,4) * t68) * t66;
t79 = t35 * t90 + (t95 / 0.2e1 + t94 / 0.2e1 + t35 / 0.2e1) * t91;
t29 = t58 * rSges(4,1) - t57 * rSges(4,2);
t25 = t54 * rSges(5,1) - t53 * rSges(5,2);
t78 = Icges(4,3) + t80;
t77 = t93 * t53;
t52 = pkin(3) * t58;
t23 = t25 + t52;
t71 = Icges(6,5) * t68 - Icges(6,6) * t66;
t70 = -t53 * rSges(6,3) + t93 * t54;
t11 = t54 * pkin(4) + t53 * pkin(8) - t70;
t9 = t11 + t52;
t10 = t53 * pkin(4) + (-rSges(6,3) - pkin(8)) * t54 - t77;
t8 = t51 + t10;
t40 = t69 * rSges(2,1) - t67 * rSges(2,2);
t39 = t67 * rSges(2,1) + t69 * rSges(2,2);
t27 = t61 * rSges(3,1) - t60 * rSges(3,2) + t64;
t26 = t60 * rSges(3,1) + t61 * rSges(3,2) + t63;
t21 = t29 + t83;
t20 = t84 + t28;
t15 = -Icges(6,3) * t53 - t71 * t54;
t14 = -Icges(6,3) * t54 + t71 * t53;
t13 = t23 + t83;
t12 = t22 + t84;
t5 = t9 + t83;
t4 = t8 + t84;
t3 = -t54 * t70 + t53 * (-t54 * rSges(6,3) - t77);
t1 = [Icges(2,3) + Icges(3,3) + m(2) * (t39 ^ 2 + t40 ^ 2) + m(3) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(6) * (t4 ^ 2 + t5 ^ 2) + t78; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t28 * t20 + t29 * t21) + m(5) * (t22 * t12 + t23 * t13) + m(6) * (t8 * t4 + t9 * t5) + t78; 0; m(4) * (t28 ^ 2 + t29 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(6) * (t8 ^ 2 + t9 ^ 2) + t78; m(5) * (t24 * t12 + t25 * t13) + m(6) * (t10 * t4 + t11 * t5) + t80; 0; m(5) * (t24 * t22 + t25 * t23) + m(6) * (t10 * t8 + t11 * t9) + t80; m(5) * (t24 ^ 2 + t25 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t80; (t4 * t54 - t5 * t53) * t87 + t79; m(6) * t3; (-t53 * t9 + t54 * t8) * t87 + t79; (t10 * t54 - t11 * t53) * t87 + t79; m(6) * (t3 ^ 2 + (t90 + t91) * t38 ^ 2) - t54 * (t90 * t14 + t15 * t92) - t53 * (t14 * t92 + t91 * t15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
