% Calculate joint inertia matrix for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:18
% EndTime: 2019-12-05 15:28:21
% DurationCPUTime: 0.64s
% Computational Cost: add. (970->92), mult. (761->135), div. (0->0), fcn. (645->6), ass. (0->47)
t52 = pkin(8) + qJ(4);
t50 = cos(t52);
t100 = t50 ^ 2;
t53 = pkin(7) + qJ(2);
t49 = sin(t53);
t46 = t49 ^ 2;
t51 = cos(t53);
t47 = t51 ^ 2;
t79 = t46 + t47;
t98 = Icges(6,4) + Icges(5,5);
t97 = Icges(5,6) - Icges(6,6);
t48 = sin(t52);
t95 = -t97 * t48 + t98 * t50;
t93 = Icges(6,2) + Icges(5,3);
t87 = rSges(6,1) + pkin(4);
t89 = t87 * t50;
t88 = rSges(6,3) + qJ(5);
t86 = -t95 * t49 + t93 * t51;
t85 = t93 * t49 + t95 * t51;
t82 = t48 * t51;
t81 = t50 * t51;
t80 = -t87 * t48 + t88 * t50;
t74 = qJ(5) * t48;
t73 = rSges(4,3) + qJ(3);
t55 = cos(pkin(8));
t44 = t55 * pkin(3) + pkin(2);
t56 = -pkin(6) - qJ(3);
t72 = t51 * t44 - t49 * t56;
t70 = t49 * rSges(6,2) + rSges(6,3) * t82 + t51 * t74 + t87 * t81;
t69 = rSges(5,1) * t50 - rSges(5,2) * t48;
t58 = rSges(5,1) * t81 - rSges(5,2) * t82 + t49 * rSges(5,3);
t54 = sin(pkin(8));
t57 = rSges(4,1) * t55 - rSges(4,2) * t54 + pkin(2);
t34 = t51 * rSges(3,1) - t49 * rSges(3,2);
t33 = -t49 * rSges(3,1) - t51 * rSges(3,2);
t32 = t48 * rSges(5,1) + t50 * rSges(5,2);
t10 = t80 * t51;
t9 = t80 * t49;
t8 = t73 * t49 + t57 * t51;
t7 = -t57 * t49 + t73 * t51;
t6 = t58 + t72;
t5 = (rSges(5,3) - t56) * t51 + (-t44 - t69) * t49;
t4 = t51 * t58 + (-t51 * rSges(5,3) + t69 * t49) * t49;
t3 = t70 + t72;
t2 = (rSges(6,2) - t56) * t51 + (-t88 * t48 - t44 - t89) * t49;
t1 = t70 * t51 + (-t51 * rSges(6,2) + (rSges(6,3) * t48 + t74 + t89) * t49) * t49;
t11 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; Icges(4,2) * t55 ^ 2 + Icges(3,3) + (Icges(4,1) * t54 + 0.2e1 * Icges(4,4) * t55) * t54 + m(3) * (t33 ^ 2 + t34 ^ 2) + m(4) * (t7 ^ 2 + t8 ^ 2) + m(5) * (t5 ^ 2 + t6 ^ 2) + m(6) * (t2 ^ 2 + t3 ^ 2) + (Icges(5,2) + Icges(6,3)) * t100 + (0.2e1 * (-Icges(6,5) + Icges(5,4)) * t50 + (Icges(5,1) + Icges(6,1)) * t48) * t48; 0; m(4) * (t49 * t7 - t51 * t8) + m(5) * (t49 * t5 - t51 * t6) + m(6) * (t49 * t2 - t51 * t3); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t79; m(5) * t4 + m(6) * t1; m(6) * (t10 * t2 + t9 * t3) + m(5) * (-t49 * t6 - t5 * t51) * t32 + t79 * (t98 * t48 + t97 * t50); m(6) * (t10 * t49 - t9 * t51); m(5) * (t79 * t32 ^ 2 + t4 ^ 2) + m(6) * (t1 ^ 2 + t10 ^ 2 + t9 ^ 2) + t85 * t49 * t46 + (t86 * t47 + (t86 * t49 + t85 * t51) * t49) * t51; -m(6) * t50; m(6) * (t2 * t51 + t3 * t49) * t48; 0; m(6) * (-t50 * t1 + (t10 * t51 + t49 * t9) * t48); m(6) * (t79 * t48 ^ 2 + t100);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t11(1), t11(2), t11(4), t11(7), t11(11); t11(2), t11(3), t11(5), t11(8), t11(12); t11(4), t11(5), t11(6), t11(9), t11(13); t11(7), t11(8), t11(9), t11(10), t11(14); t11(11), t11(12), t11(13), t11(14), t11(15);];
Mq = res;
