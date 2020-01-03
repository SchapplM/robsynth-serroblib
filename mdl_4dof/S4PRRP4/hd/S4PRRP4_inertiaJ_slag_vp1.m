% Calculate joint inertia matrix for
% S4PRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP4_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP4_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:46
% EndTime: 2019-12-31 16:27:47
% DurationCPUTime: 0.47s
% Computational Cost: add. (563->70), mult. (617->101), div. (0->0), fcn. (535->4), ass. (0->41)
t50 = cos(qJ(3));
t94 = t50 ^ 2;
t47 = pkin(6) + qJ(2);
t45 = sin(t47);
t43 = t45 ^ 2;
t46 = cos(t47);
t44 = t46 ^ 2;
t70 = t43 + t44;
t92 = Icges(5,4) + Icges(4,5);
t91 = Icges(4,6) - Icges(5,6);
t49 = sin(qJ(3));
t89 = -t91 * t49 + t92 * t50;
t87 = Icges(5,2) + Icges(4,3);
t81 = rSges(5,1) + pkin(3);
t83 = t81 * t50;
t82 = rSges(5,3) + qJ(4);
t80 = -t89 * t45 + t87 * t46;
t79 = t87 * t45 + t89 * t46;
t76 = t46 * rSges(5,2);
t75 = t46 * rSges(4,3);
t74 = t46 * t49;
t73 = t46 * t50;
t72 = -t49 * t81 + t82 * t50;
t71 = t46 * pkin(2) + t45 * pkin(5);
t65 = qJ(4) * t49;
t63 = t45 * rSges(5,2) + rSges(5,3) * t74 + t46 * t65 + t73 * t81;
t61 = rSges(4,1) * t50 - rSges(4,2) * t49;
t51 = rSges(4,1) * t73 - rSges(4,2) * t74 + t45 * rSges(4,3);
t41 = t46 * pkin(5);
t32 = t49 * rSges(4,1) + t50 * rSges(4,2);
t22 = t46 * rSges(3,1) - t45 * rSges(3,2);
t21 = -t45 * rSges(3,1) - t46 * rSges(3,2);
t8 = t72 * t46;
t7 = t72 * t45;
t6 = t51 + t71;
t5 = t75 + t41 + (-pkin(2) - t61) * t45;
t4 = t63 + t71;
t3 = t76 + t41 + (-t49 * t82 - pkin(2) - t83) * t45;
t2 = t46 * t51 + (t61 * t45 - t75) * t45;
t1 = t63 * t46 + (-t76 + (rSges(5,3) * t49 + t65 + t83) * t45) * t45;
t9 = [m(2) + m(3) + m(4) + m(5); 0; Icges(3,3) + m(3) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + (Icges(4,2) + Icges(5,3)) * t94 + (0.2e1 * (-Icges(5,5) + Icges(4,4)) * t50 + (Icges(4,1) + Icges(5,1)) * t49) * t49; m(4) * t2 + m(5) * t1; m(5) * (t8 * t3 + t7 * t4) + m(4) * (-t45 * t6 - t46 * t5) * t32 + t70 * (t92 * t49 + t91 * t50); m(4) * (t70 * t32 ^ 2 + t2 ^ 2) + m(5) * (t1 ^ 2 + t7 ^ 2 + t8 ^ 2) + t79 * t45 * t43 + (t80 * t44 + (t80 * t45 + t79 * t46) * t45) * t46; -m(5) * t50; m(5) * (t3 * t46 + t4 * t45) * t49; m(5) * (-t50 * t1 + (t45 * t7 + t46 * t8) * t49); m(5) * (t70 * t49 ^ 2 + t94);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t9(1), t9(2), t9(4), t9(7); t9(2), t9(3), t9(5), t9(8); t9(4), t9(5), t9(6), t9(9); t9(7), t9(8), t9(9), t9(10);];
Mq = res;
