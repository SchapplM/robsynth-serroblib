% Calculate joint inertia matrix for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:43
% EndTime: 2019-12-05 17:04:44
% DurationCPUTime: 0.36s
% Computational Cost: add. (937->85), mult. (706->119), div. (0->0), fcn. (566->8), ass. (0->51)
t52 = qJ(2) + qJ(3);
t50 = qJ(4) + t52;
t46 = sin(t50);
t78 = t46 ^ 2;
t47 = cos(t50);
t77 = t47 ^ 2;
t53 = sin(qJ(5));
t82 = Icges(6,5) * t53;
t55 = cos(qJ(5));
t81 = Icges(6,6) * t55;
t33 = t81 + t82;
t80 = rSges(6,1) * t55 - rSges(6,2) * t53;
t79 = t46 * t47;
t36 = t53 * rSges(6,1) + t55 * rSges(6,2);
t74 = m(6) * t36;
t48 = sin(t52);
t73 = pkin(3) * t48;
t54 = sin(qJ(2));
t72 = t54 * pkin(2);
t67 = Icges(6,2) * t55 ^ 2 + Icges(5,3) + (Icges(6,1) * t53 + 0.2e1 * Icges(6,4) * t55) * t53;
t66 = t33 * t77 + (t82 / 0.2e1 + t81 / 0.2e1 + t33 / 0.2e1) * t78;
t49 = cos(t52);
t27 = t49 * rSges(4,1) - t48 * rSges(4,2);
t25 = t47 * rSges(5,1) - t46 * rSges(5,2);
t65 = Icges(4,3) + t67;
t45 = pkin(3) * t49;
t21 = t25 + t45;
t26 = -t48 * rSges(4,1) - t49 * rSges(4,2);
t24 = -t46 * rSges(5,1) - t47 * rSges(5,2);
t59 = Icges(6,5) * t55 - Icges(6,6) * t53;
t58 = -t47 * rSges(6,3) + t80 * t46;
t57 = t46 * rSges(6,3) + t80 * t47;
t11 = t47 * pkin(6) - t58;
t10 = t46 * pkin(6) + t57;
t9 = t10 + t45;
t20 = t24 - t73;
t8 = t11 - t73;
t56 = cos(qJ(2));
t51 = t56 * pkin(2);
t38 = t56 * rSges(3,1) - t54 * rSges(3,2);
t37 = -t54 * rSges(3,1) - t56 * rSges(3,2);
t23 = t27 + t51;
t22 = t26 - t72;
t19 = t21 + t51;
t18 = t20 - t72;
t13 = Icges(6,3) * t46 + t59 * t47;
t12 = -Icges(6,3) * t47 + t59 * t46;
t7 = t51 + t9;
t6 = t8 - t72;
t3 = t46 * t58 + t47 * t57;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; Icges(3,3) + m(6) * (t6 ^ 2 + t7 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(3) * (t37 ^ 2 + t38 ^ 2) + t65; 0; m(6) * (t8 * t6 + t9 * t7) + m(5) * (t20 * t18 + t21 * t19) + m(4) * (t26 * t22 + t27 * t23) + t65; m(4) * (t26 ^ 2 + t27 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(6) * (t8 ^ 2 + t9 ^ 2) + t65; 0; m(6) * (t10 * t7 + t11 * t6) + m(5) * (t24 * t18 + t25 * t19) + t67; m(5) * (t24 * t20 + t25 * t21) + m(6) * (t10 * t9 + t11 * t8) + t67; m(5) * (t24 ^ 2 + t25 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t67; m(6) * t3; (-t46 * t7 - t47 * t6) * t74 + t66; (-t46 * t9 - t47 * t8) * t74 + t66; (-t10 * t46 - t11 * t47) * t74 + t66; m(6) * (t3 ^ 2 + (t77 + t78) * t36 ^ 2) + t46 * (-t12 * t79 + t78 * t13) - t47 * (t77 * t12 - t13 * t79);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
