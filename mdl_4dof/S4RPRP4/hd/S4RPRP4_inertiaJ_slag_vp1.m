% Calculate joint inertia matrix for
% S4RPRP4
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP4_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP4_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:38
% EndTime: 2019-12-31 16:43:39
% DurationCPUTime: 0.50s
% Computational Cost: add. (579->78), mult. (639->110), div. (0->0), fcn. (551->6), ass. (0->47)
t54 = cos(qJ(3));
t101 = t54 ^ 2;
t50 = qJ(1) + pkin(6);
t47 = sin(t50);
t45 = t47 ^ 2;
t48 = cos(t50);
t46 = t48 ^ 2;
t77 = t45 + t46;
t99 = Icges(5,4) + Icges(4,5);
t98 = Icges(4,6) - Icges(5,6);
t52 = sin(qJ(3));
t96 = -t98 * t52 + t99 * t54;
t94 = Icges(5,2) + Icges(4,3);
t88 = rSges(5,1) + pkin(3);
t90 = t88 * t54;
t89 = rSges(5,3) + qJ(4);
t87 = -t96 * t47 + t94 * t48;
t86 = t94 * t47 + t96 * t48;
t53 = sin(qJ(1));
t83 = t53 * pkin(1);
t82 = t48 * rSges(5,2);
t81 = t48 * rSges(4,3);
t80 = t48 * t52;
t79 = t48 * t54;
t78 = -t88 * t52 + t89 * t54;
t72 = qJ(4) * t52;
t55 = cos(qJ(1));
t49 = t55 * pkin(1);
t71 = t48 * pkin(2) + t47 * pkin(5) + t49;
t70 = t48 * pkin(5) - t83;
t68 = t47 * rSges(5,2) + rSges(5,3) * t80 + t48 * t72 + t88 * t79;
t66 = rSges(4,1) * t54 - rSges(4,2) * t52;
t56 = rSges(4,1) * t79 - rSges(4,2) * t80 + t47 * rSges(4,3);
t34 = t55 * rSges(2,1) - t53 * rSges(2,2);
t33 = -t53 * rSges(2,1) - t55 * rSges(2,2);
t32 = t52 * rSges(4,1) + t54 * rSges(4,2);
t22 = t48 * rSges(3,1) - t47 * rSges(3,2) + t49;
t21 = -t47 * rSges(3,1) - t48 * rSges(3,2) - t83;
t8 = t78 * t48;
t7 = t78 * t47;
t6 = t56 + t71;
t5 = t81 + (-pkin(2) - t66) * t47 + t70;
t4 = t68 + t71;
t3 = t82 + (-t89 * t52 - pkin(2) - t90) * t47 + t70;
t2 = t48 * t56 + (t66 * t47 - t81) * t47;
t1 = t68 * t48 + (-t82 + (rSges(5,3) * t52 + t72 + t90) * t47) * t47;
t9 = [Icges(2,3) + Icges(3,3) + m(2) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + (Icges(4,2) + Icges(5,3)) * t101 + (0.2e1 * (-Icges(5,5) + Icges(4,4)) * t54 + (Icges(4,1) + Icges(5,1)) * t52) * t52; 0; m(3) + m(4) + m(5); m(5) * (t8 * t3 + t7 * t4) + m(4) * (-t47 * t6 - t48 * t5) * t32 + t77 * (t99 * t52 + t98 * t54); m(4) * t2 + m(5) * t1; m(4) * (t77 * t32 ^ 2 + t2 ^ 2) + m(5) * (t1 ^ 2 + t7 ^ 2 + t8 ^ 2) + t86 * t47 * t45 + (t87 * t46 + (t87 * t47 + t86 * t48) * t47) * t48; m(5) * (t3 * t48 + t4 * t47) * t52; -m(5) * t54; m(5) * (-t54 * t1 + (t47 * t7 + t48 * t8) * t52); m(5) * (t77 * t52 ^ 2 + t101);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t9(1), t9(2), t9(4), t9(7); t9(2), t9(3), t9(5), t9(8); t9(4), t9(5), t9(6), t9(9); t9(7), t9(8), t9(9), t9(10);];
Mq = res;
