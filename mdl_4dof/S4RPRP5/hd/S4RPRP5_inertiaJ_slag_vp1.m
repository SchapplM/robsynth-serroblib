% Calculate joint inertia matrix for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP5_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP5_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:46
% EndTime: 2019-12-31 16:44:48
% DurationCPUTime: 0.63s
% Computational Cost: add. (577->85), mult. (716->132), div. (0->0), fcn. (610->6), ass. (0->46)
t48 = pkin(6) + qJ(3);
t45 = cos(t48);
t99 = t45 ^ 2;
t54 = sin(qJ(1));
t49 = t54 ^ 2;
t55 = cos(qJ(1));
t50 = t55 ^ 2;
t78 = t49 + t50;
t97 = Icges(5,4) + Icges(4,5);
t96 = Icges(4,6) - Icges(5,6);
t44 = sin(t48);
t94 = -t96 * t44 + t97 * t45;
t92 = Icges(5,2) + Icges(4,3);
t86 = rSges(5,1) + pkin(3);
t88 = t86 * t45;
t87 = rSges(5,3) + qJ(4);
t85 = -t94 * t54 + t92 * t55;
t84 = t92 * t54 + t94 * t55;
t81 = t44 * t55;
t80 = t45 * t55;
t79 = -t86 * t44 + t87 * t45;
t73 = qJ(4) * t44;
t72 = rSges(3,3) + qJ(2);
t52 = cos(pkin(6));
t42 = t52 * pkin(2) + pkin(1);
t53 = -pkin(5) - qJ(2);
t71 = t55 * t42 - t54 * t53;
t69 = t54 * rSges(5,2) + rSges(5,3) * t81 + t55 * t73 + t86 * t80;
t68 = rSges(4,1) * t45 - rSges(4,2) * t44;
t57 = rSges(4,1) * t80 - rSges(4,2) * t81 + t54 * rSges(4,3);
t51 = sin(pkin(6));
t56 = rSges(3,1) * t52 - rSges(3,2) * t51 + pkin(1);
t34 = t55 * rSges(2,1) - t54 * rSges(2,2);
t33 = -t54 * rSges(2,1) - t55 * rSges(2,2);
t32 = t44 * rSges(4,1) + t45 * rSges(4,2);
t10 = t72 * t54 + t56 * t55;
t9 = -t56 * t54 + t72 * t55;
t8 = t79 * t55;
t7 = t79 * t54;
t6 = t57 + t71;
t5 = (rSges(4,3) - t53) * t55 + (-t42 - t68) * t54;
t4 = t55 * t57 + (-t55 * rSges(4,3) + t68 * t54) * t54;
t3 = t69 + t71;
t2 = (rSges(5,2) - t53) * t55 + (-t87 * t44 - t42 - t88) * t54;
t1 = t69 * t55 + (-t55 * rSges(5,2) + (rSges(5,3) * t44 + t73 + t88) * t54) * t54;
t11 = [Icges(3,2) * t52 ^ 2 + Icges(2,3) + (Icges(3,1) * t51 + 0.2e1 * Icges(3,4) * t52) * t51 + m(2) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t10 ^ 2 + t9 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t2 ^ 2 + t3 ^ 2) + (Icges(4,2) + Icges(5,3)) * t99 + (0.2e1 * (-Icges(5,5) + Icges(4,4)) * t45 + (Icges(4,1) + Icges(5,1)) * t44) * t44; m(3) * (-t55 * t10 + t54 * t9) + m(4) * (t54 * t5 - t55 * t6) + m(5) * (t54 * t2 - t55 * t3); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * t78; m(5) * (t8 * t2 + t7 * t3) + m(4) * (-t5 * t55 - t54 * t6) * t32 + t78 * (t97 * t44 + t96 * t45); m(5) * (t8 * t54 - t7 * t55); m(4) * (t78 * t32 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t7 ^ 2 + t8 ^ 2) + t84 * t54 * t49 + (t85 * t50 + (t54 * t85 + t55 * t84) * t54) * t55; m(5) * (t2 * t55 + t3 * t54) * t44; 0; m(5) * (-t45 * t1 + (t54 * t7 + t55 * t8) * t44); m(5) * (t78 * t44 ^ 2 + t99);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1), t11(2), t11(4), t11(7); t11(2), t11(3), t11(5), t11(8); t11(4), t11(5), t11(6), t11(9); t11(7), t11(8), t11(9), t11(10);];
Mq = res;
