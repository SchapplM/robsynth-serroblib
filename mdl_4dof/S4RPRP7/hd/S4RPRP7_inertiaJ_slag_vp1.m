% Calculate joint inertia matrix for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP7_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP7_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:03
% EndTime: 2019-12-31 16:47:04
% DurationCPUTime: 0.55s
% Computational Cost: add. (317->80), mult. (718->117), div. (0->0), fcn. (602->4), ass. (0->46)
t51 = cos(qJ(3));
t100 = t51 ^ 2;
t50 = sin(qJ(1));
t46 = t50 ^ 2;
t52 = cos(qJ(1));
t48 = t52 ^ 2;
t77 = t46 + t48;
t98 = Icges(5,4) + Icges(4,5);
t97 = Icges(4,6) - Icges(5,6);
t49 = sin(qJ(3));
t72 = rSges(5,3) + qJ(4);
t68 = t72 * t51;
t89 = rSges(5,1) + pkin(3);
t96 = (t49 * t89 - t68) * t52;
t94 = t98 * t49 + t97 * t51;
t92 = Icges(5,2) + Icges(4,3);
t88 = (rSges(4,1) * t49 + rSges(4,2) * t51) * t52;
t87 = t94 * t50 + t92 * t52;
t86 = t92 * t50 - t94 * t52;
t83 = -pkin(1) - pkin(5);
t82 = m(5) * t51;
t81 = t49 * t50;
t80 = t50 * t51;
t79 = t72 * t49 + t89 * t51;
t78 = t52 * pkin(1) + t50 * qJ(2);
t71 = rSges(4,1) * t81 + rSges(4,2) * t80 + t52 * rSges(4,3);
t70 = -t52 * rSges(5,2) - t89 * t81;
t69 = t52 * pkin(5) + t78;
t41 = t52 * qJ(2);
t2 = t41 + (-rSges(5,2) + t83) * t50 + t96;
t3 = -t50 * t68 + t69 - t70;
t66 = t50 * t2 - t52 * t3;
t7 = t79 * t50;
t8 = t79 * t52;
t65 = t7 * t50 + t8 * t52;
t5 = t41 + t88 + (-rSges(4,3) + t83) * t50;
t6 = t69 + t71;
t53 = m(4) * (t50 * t5 - t52 * t6);
t34 = t52 * rSges(2,1) - t50 * rSges(2,2);
t33 = t51 * rSges(4,1) - t49 * rSges(4,2);
t30 = -t50 * rSges(2,1) - t52 * rSges(2,2);
t22 = -t52 * rSges(3,2) + t50 * rSges(3,3) + t78;
t21 = t52 * rSges(3,3) + t41 + (rSges(3,2) - pkin(1)) * t50;
t4 = -t50 * t71 + (t50 * rSges(4,3) - t88) * t52;
t1 = (t72 * t80 + t70) * t50 + (t50 * rSges(5,2) - t96) * t52;
t9 = [Icges(3,1) + Icges(2,3) + m(2) * (t30 ^ 2 + t34 ^ 2) + m(3) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t2 ^ 2 + t3 ^ 2) + (Icges(4,1) + Icges(5,1)) * t100 + (0.2e1 * (Icges(5,5) - Icges(4,4)) * t51 + (Icges(4,2) + Icges(5,3)) * t49) * t49; m(3) * (t50 * t21 - t52 * t22) + t53 + m(5) * t66; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * t77; m(5) * (t7 * t2 - t8 * t3) + t33 * t53 + t77 * (-t97 * t49 + t98 * t51); m(4) * t77 * t33 + m(5) * t65; m(4) * (t77 * t33 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t7 ^ 2 + t8 ^ 2) + t86 * t50 * t46 + (t87 * t48 + (t50 * t87 + t52 * t86) * t50) * t52; -t66 * t82; -t77 * t82; m(5) * (t49 * t1 - t65 * t51); m(5) * (t77 * t100 + t49 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t9(1), t9(2), t9(4), t9(7); t9(2), t9(3), t9(5), t9(8); t9(4), t9(5), t9(6), t9(9); t9(7), t9(8), t9(9), t9(10);];
Mq = res;
