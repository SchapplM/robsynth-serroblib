% Calculate joint inertia matrix for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:20
% EndTime: 2019-12-05 17:49:22
% DurationCPUTime: 0.39s
% Computational Cost: add. (1114->94), mult. (674->134), div. (0->0), fcn. (550->10), ass. (0->57)
t58 = qJ(1) + pkin(8);
t56 = qJ(3) + t58;
t49 = sin(t56);
t47 = t49 ^ 2;
t50 = cos(t56);
t48 = t50 ^ 2;
t57 = pkin(9) + qJ(5);
t52 = sin(t57);
t93 = Icges(6,5) * t52;
t54 = cos(t57);
t92 = Icges(6,6) * t54;
t29 = t92 + t93;
t59 = sin(pkin(9));
t60 = cos(pkin(9));
t91 = -rSges(5,1) * t60 + rSges(5,2) * t59 - pkin(3);
t90 = rSges(5,3) + qJ(4);
t89 = t49 * t50;
t35 = t52 * rSges(6,1) + t54 * rSges(6,2);
t86 = m(6) * t35;
t62 = sin(qJ(1));
t85 = t62 * pkin(1);
t63 = cos(qJ(1));
t84 = t63 * pkin(1);
t83 = rSges(6,1) * t54;
t81 = rSges(6,2) * t52;
t80 = t50 * rSges(6,3) + t49 * t81;
t79 = t48 + t47;
t76 = t29 * t48 + (t93 / 0.2e1 + t92 / 0.2e1 + t29 / 0.2e1) * t47;
t27 = -t50 * rSges(4,1) + t49 * rSges(4,2);
t74 = -t60 * pkin(4) - pkin(3) - t83;
t73 = -t49 * rSges(6,3) + t50 * t81;
t53 = sin(t58);
t72 = -pkin(2) * t53 - t85;
t55 = cos(t58);
t71 = -pkin(2) * t55 - t84;
t70 = Icges(5,2) * t60 ^ 2 + Icges(6,2) * t54 ^ 2 + Icges(4,3) + (Icges(5,1) * t59 + 0.2e1 * Icges(5,4) * t60) * t59 + (Icges(6,1) * t52 + 0.2e1 * Icges(6,4) * t54) * t52;
t26 = -t49 * rSges(4,1) - t50 * rSges(4,2);
t64 = Icges(6,5) * t54 - Icges(6,6) * t52;
t12 = t91 * t49 + t90 * t50;
t61 = -pkin(7) - qJ(4);
t11 = t49 * t61 + t74 * t50 + t73;
t10 = t74 * t49 - t50 * t61 + t80;
t13 = -t90 * t49 + t91 * t50;
t42 = -t63 * rSges(2,1) + t62 * rSges(2,2);
t41 = -t62 * rSges(2,1) - t63 * rSges(2,2);
t25 = -t55 * rSges(3,1) + t53 * rSges(3,2) - t84;
t24 = -t53 * rSges(3,1) - t55 * rSges(3,2) - t85;
t21 = t27 + t71;
t20 = t26 + t72;
t15 = Icges(6,3) * t49 + t64 * t50;
t14 = Icges(6,3) * t50 - t64 * t49;
t9 = t13 + t71;
t8 = t12 + t72;
t7 = t11 + t71;
t6 = t10 + t72;
t3 = t50 * (t50 * t83 - t73) - t49 * (-t49 * t83 + t80);
t1 = [Icges(2,3) + Icges(3,3) + m(2) * (t41 ^ 2 + t42 ^ 2) + m(3) * (t24 ^ 2 + t25 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) + t70; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t26 * t20 + t27 * t21) + m(5) * (t12 * t8 + t13 * t9) + m(6) * (t10 * t6 + t11 * t7) + t70; 0; m(5) * (t12 ^ 2 + t13 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + m(4) * (t26 ^ 2 + t27 ^ 2) + t70; m(5) * (t49 * t8 + t50 * t9) + m(6) * (t49 * t6 + t50 * t7); 0; m(5) * (t49 * t12 + t50 * t13) + m(6) * (t49 * t10 + t50 * t11); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t79; (t49 * t7 - t50 * t6) * t86 + t76; m(6) * t3; (-t10 * t50 + t11 * t49) * t86 + t76; 0; m(6) * (t79 * t35 ^ 2 + t3 ^ 2) + t50 * (t48 * t14 + t15 * t89) + t49 * (t14 * t89 + t47 * t15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
