% Calculate joint inertia matrix for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:04
% EndTime: 2019-12-05 17:06:06
% DurationCPUTime: 0.40s
% Computational Cost: add. (1375->90), mult. (730->122), div. (0->0), fcn. (590->8), ass. (0->53)
t56 = pkin(9) + qJ(2);
t55 = qJ(3) + t56;
t52 = qJ(4) + t55;
t47 = sin(t52);
t80 = t47 ^ 2;
t48 = cos(t52);
t79 = t48 ^ 2;
t57 = sin(qJ(5));
t83 = Icges(6,5) * t57;
t58 = cos(qJ(5));
t82 = Icges(6,6) * t58;
t35 = t82 + t83;
t81 = t47 * t48;
t38 = t57 * rSges(6,1) + t58 * rSges(6,2);
t76 = m(6) * t38;
t53 = sin(t56);
t75 = pkin(2) * t53;
t50 = sin(t55);
t74 = pkin(3) * t50;
t73 = rSges(6,1) * t58;
t72 = rSges(6,2) * t57;
t71 = t48 * rSges(6,3) + t47 * t72;
t68 = Icges(6,2) * t58 ^ 2 + Icges(5,3) + (Icges(6,1) * t57 + 0.2e1 * Icges(6,4) * t58) * t57;
t67 = t35 * t79 + (t83 / 0.2e1 + t82 / 0.2e1 + t35 / 0.2e1) * t80;
t51 = cos(t55);
t27 = t51 * rSges(4,1) - t50 * rSges(4,2);
t25 = t48 * rSges(5,1) - t47 * rSges(5,2);
t66 = Icges(4,3) + t68;
t46 = pkin(3) * t51;
t21 = t25 + t46;
t26 = -t50 * rSges(4,1) - t51 * rSges(4,2);
t24 = -t47 * rSges(5,1) - t48 * rSges(5,2);
t60 = Icges(6,5) * t58 - Icges(6,6) * t57;
t59 = t47 * rSges(6,3) + (-t72 + t73) * t48;
t11 = t48 * pkin(4) + t47 * pkin(8) + t59;
t20 = t24 - t74;
t9 = t11 + t46;
t10 = t48 * pkin(8) + (-pkin(4) - t73) * t47 + t71;
t8 = t10 - t74;
t54 = cos(t56);
t49 = pkin(2) * t54;
t29 = t54 * rSges(3,1) - t53 * rSges(3,2);
t28 = -t53 * rSges(3,1) - t54 * rSges(3,2);
t23 = t27 + t49;
t22 = t26 - t75;
t19 = t21 + t49;
t18 = t20 - t75;
t13 = Icges(6,3) * t47 + t60 * t48;
t12 = -Icges(6,3) * t48 + t60 * t47;
t7 = t49 + t9;
t6 = t8 - t75;
t3 = t47 * (t47 * t73 - t71) + t48 * t59;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; Icges(3,3) + m(6) * (t6 ^ 2 + t7 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(3) * (t28 ^ 2 + t29 ^ 2) + t66; 0; m(6) * (t8 * t6 + t9 * t7) + m(5) * (t20 * t18 + t21 * t19) + m(4) * (t26 * t22 + t27 * t23) + t66; m(4) * (t26 ^ 2 + t27 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(6) * (t8 ^ 2 + t9 ^ 2) + t66; 0; m(6) * (t10 * t6 + t11 * t7) + m(5) * (t24 * t18 + t25 * t19) + t68; m(5) * (t24 * t20 + t25 * t21) + m(6) * (t10 * t8 + t11 * t9) + t68; m(5) * (t24 ^ 2 + t25 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t68; m(6) * t3; (-t47 * t7 - t48 * t6) * t76 + t67; (-t47 * t9 - t48 * t8) * t76 + t67; (-t10 * t48 - t11 * t47) * t76 + t67; m(6) * (t3 ^ 2 + (t79 + t80) * t38 ^ 2) + t47 * (-t12 * t81 + t80 * t13) - t48 * (t79 * t12 - t13 * t81);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
