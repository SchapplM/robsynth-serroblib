% Calculate joint inertia matrix for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:29
% EndTime: 2019-12-31 17:56:30
% DurationCPUTime: 0.43s
% Computational Cost: add. (1068->102), mult. (1068->152), div. (0->0), fcn. (1170->8), ass. (0->55)
t50 = sin(qJ(5));
t86 = t50 / 0.2e1;
t52 = cos(qJ(5));
t85 = t52 / 0.2e1;
t84 = Icges(6,5) * t86 + Icges(6,6) * t85;
t49 = qJ(1) + pkin(8);
t46 = sin(t49);
t47 = cos(t49);
t76 = sin(qJ(4));
t77 = cos(qJ(4));
t28 = -t46 * t76 - t47 * t77;
t29 = -t46 * t77 + t47 * t76;
t83 = t29 * t28;
t74 = rSges(6,2) * t50;
t66 = -pkin(4) + t74;
t75 = rSges(6,1) * t52;
t72 = t29 * rSges(6,3) - t28 * t75;
t7 = -t29 * pkin(7) - t66 * t28 - t72;
t73 = t28 * rSges(6,3) + t29 * t75;
t6 = -t28 * pkin(7) + t66 * t29 - t73;
t70 = Icges(6,4) * t52;
t57 = Icges(6,2) * t50 - t70;
t71 = Icges(6,4) * t50;
t58 = -Icges(6,1) * t52 + t71;
t82 = -(t28 * t84 - t52 * (-Icges(6,6) * t28 + t57 * t29) / 0.2e1 - t50 * (-Icges(6,5) * t28 + t58 * t29) / 0.2e1) * t28 - (t29 * t84 + (Icges(6,6) * t29 + t57 * t28) * t85 + (Icges(6,5) * t29 + t58 * t28) * t86) * t29;
t81 = t28 ^ 2;
t80 = t29 ^ 2;
t36 = -t50 * rSges(6,1) - t52 * rSges(6,2);
t79 = m(6) * t36;
t51 = sin(qJ(1));
t78 = t51 * pkin(1);
t53 = cos(qJ(1));
t48 = t53 * pkin(1);
t67 = t47 * pkin(2) + t46 * qJ(3) + t48;
t65 = t47 * qJ(3) - t78;
t64 = t47 * pkin(3) + t67;
t18 = -t29 * rSges(5,1) + t28 * rSges(5,2);
t19 = t28 * rSges(5,1) + t29 * rSges(5,2);
t56 = -Icges(6,5) * t52 + Icges(6,6) * t50;
t55 = t52 * (-Icges(6,2) * t52 - t71) + t50 * (-Icges(6,1) * t50 - t70) - Icges(5,3);
t54 = (-pkin(2) - pkin(3)) * t46 + t65;
t38 = t53 * rSges(2,1) - t51 * rSges(2,2);
t37 = -t51 * rSges(2,1) - t53 * rSges(2,2);
t27 = t47 * rSges(3,1) - t46 * rSges(3,2) + t48;
t26 = -t46 * rSges(3,1) - t47 * rSges(3,2) - t78;
t21 = t47 * rSges(4,1) + t46 * rSges(4,3) + t67;
t20 = t47 * rSges(4,3) + (-rSges(4,1) - pkin(2)) * t46 + t65;
t17 = -t19 + t64;
t16 = -t18 + t54;
t11 = Icges(6,3) * t29 + t56 * t28;
t10 = -Icges(6,3) * t28 + t56 * t29;
t5 = t64 - t7;
t4 = t54 - t6;
t1 = t29 * (t29 * t74 - t73) + t28 * (t28 * t74 + t72);
t2 = [Icges(4,2) + Icges(2,3) + Icges(3,3) + m(6) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(3) * (t26 ^ 2 + t27 ^ 2) + m(2) * (t37 ^ 2 + t38 ^ 2) - t55; 0; m(3) + m(4) + m(5) + m(6); m(6) * (t46 * t4 - t47 * t5) + m(5) * (t46 * t16 - t47 * t17) + m(4) * (t46 * t20 - t47 * t21); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t46 ^ 2 + t47 ^ 2); m(6) * (t6 * t4 + t7 * t5) + m(5) * (t18 * t16 + t19 * t17) + t55; 0; m(5) * (t18 * t46 - t19 * t47) + m(6) * (t6 * t46 - t7 * t47); m(5) * (t18 ^ 2 + t19 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) - t55; (-t28 * t4 - t29 * t5) * t79 + t82; m(6) * t1; (-t28 * t46 + t29 * t47) * t79; (-t28 * t6 - t29 * t7) * t79 - t82; m(6) * (t1 ^ 2 + (t80 + t81) * t36 ^ 2) + t29 * (-t10 * t83 + t80 * t11) - t28 * (t81 * t10 - t11 * t83);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
