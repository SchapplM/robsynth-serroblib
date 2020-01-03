% Calculate joint inertia matrix for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:08
% EndTime: 2019-12-31 17:22:09
% DurationCPUTime: 0.35s
% Computational Cost: add. (984->85), mult. (714->121), div. (0->0), fcn. (578->8), ass. (0->52)
t53 = qJ(1) + qJ(2);
t51 = qJ(3) + t53;
t47 = sin(t51);
t79 = t47 ^ 2;
t48 = cos(t51);
t78 = t48 ^ 2;
t54 = sin(qJ(4));
t82 = Icges(5,5) * t54;
t56 = cos(qJ(4));
t81 = Icges(5,6) * t56;
t33 = t81 + t82;
t80 = t47 * t48;
t36 = t54 * rSges(5,1) + t56 * rSges(5,2);
t75 = m(5) * t36;
t49 = sin(t53);
t74 = pkin(2) * t49;
t55 = sin(qJ(1));
t73 = t55 * pkin(1);
t72 = rSges(5,1) * t56;
t71 = rSges(5,2) * t54;
t70 = t48 * rSges(5,3) + t47 * t71;
t67 = Icges(5,2) * t56 ^ 2 + Icges(4,3) + (Icges(5,1) * t54 + 0.2e1 * Icges(5,4) * t56) * t54;
t66 = t33 * t78 + (t82 / 0.2e1 + t81 / 0.2e1 + t33 / 0.2e1) * t79;
t50 = cos(t53);
t27 = t50 * rSges(3,1) - t49 * rSges(3,2);
t25 = t48 * rSges(4,1) - t47 * rSges(4,2);
t65 = Icges(3,3) + t67;
t46 = pkin(2) * t50;
t21 = t25 + t46;
t26 = -t49 * rSges(3,1) - t50 * rSges(3,2);
t24 = -t47 * rSges(4,1) - t48 * rSges(4,2);
t59 = Icges(5,5) * t56 - Icges(5,6) * t54;
t58 = t47 * rSges(5,3) + (-t71 + t72) * t48;
t11 = t48 * pkin(3) + t47 * pkin(7) + t58;
t20 = t24 - t74;
t9 = t11 + t46;
t10 = t48 * pkin(7) + (-pkin(3) - t72) * t47 + t70;
t8 = t10 - t74;
t57 = cos(qJ(1));
t52 = t57 * pkin(1);
t38 = t57 * rSges(2,1) - t55 * rSges(2,2);
t37 = -t55 * rSges(2,1) - t57 * rSges(2,2);
t23 = t27 + t52;
t22 = t26 - t73;
t19 = t21 + t52;
t18 = t20 - t73;
t13 = Icges(5,3) * t47 + t59 * t48;
t12 = -Icges(5,3) * t48 + t59 * t47;
t7 = t52 + t9;
t6 = t8 - t73;
t3 = t47 * (t47 * t72 - t70) + t48 * t58;
t1 = [Icges(2,3) + m(2) * (t37 ^ 2 + t38 ^ 2) + m(3) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2) + t65; m(3) * (t26 * t22 + t27 * t23) + m(4) * (t20 * t18 + t21 * t19) + m(5) * (t8 * t6 + t9 * t7) + t65; m(3) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2) + t65; m(4) * (t24 * t18 + t25 * t19) + m(5) * (t10 * t6 + t11 * t7) + t67; m(4) * (t24 * t20 + t25 * t21) + m(5) * (t10 * t8 + t11 * t9) + t67; m(4) * (t24 ^ 2 + t25 ^ 2) + m(5) * (t10 ^ 2 + t11 ^ 2) + t67; (-t47 * t7 - t48 * t6) * t75 + t66; (-t47 * t9 - t48 * t8) * t75 + t66; (-t10 * t48 - t11 * t47) * t75 + t66; m(5) * (t3 ^ 2 + (t78 + t79) * t36 ^ 2) + t47 * (-t12 * t80 + t79 * t13) - t48 * (t78 * t12 - t13 * t80);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
