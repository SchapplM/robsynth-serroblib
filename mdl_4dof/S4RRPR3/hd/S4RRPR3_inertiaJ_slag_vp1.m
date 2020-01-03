% Calculate joint inertia matrix for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:28
% EndTime: 2019-12-31 17:01:29
% DurationCPUTime: 0.32s
% Computational Cost: add. (669->72), mult. (494->101), div. (0->0), fcn. (400->8), ass. (0->47)
t49 = qJ(1) + qJ(2);
t45 = pkin(7) + t49;
t42 = sin(t45);
t74 = t42 ^ 2;
t43 = cos(t45);
t73 = t43 ^ 2;
t50 = sin(qJ(4));
t77 = Icges(5,5) * t50;
t52 = cos(qJ(4));
t76 = Icges(5,6) * t52;
t29 = t76 + t77;
t75 = t42 * t43;
t32 = t50 * rSges(5,1) + t52 * rSges(5,2);
t70 = m(5) * t32;
t46 = sin(t49);
t69 = pkin(2) * t46;
t51 = sin(qJ(1));
t68 = t51 * pkin(1);
t67 = rSges(5,1) * t52;
t66 = rSges(5,2) * t50;
t65 = t43 * rSges(5,3) + t42 * t66;
t62 = t29 * t73 + (t77 / 0.2e1 + t76 / 0.2e1 + t29 / 0.2e1) * t74;
t47 = cos(t49);
t23 = t47 * rSges(3,1) - t46 * rSges(3,2);
t61 = Icges(5,2) * t52 ^ 2 + Icges(3,3) + Icges(4,3) + (Icges(5,1) * t50 + 0.2e1 * Icges(5,4) * t52) * t50;
t44 = pkin(2) * t47;
t19 = t43 * rSges(4,1) - t42 * rSges(4,2) + t44;
t22 = -t46 * rSges(3,1) - t47 * rSges(3,2);
t55 = Icges(5,5) * t52 - Icges(5,6) * t50;
t54 = t42 * rSges(5,3) + (-t66 + t67) * t43;
t18 = -t42 * rSges(4,1) - t43 * rSges(4,2) - t69;
t9 = t43 * pkin(3) + t42 * pkin(6) + t44 + t54;
t8 = -t69 + t43 * pkin(6) + (-pkin(3) - t67) * t42 + t65;
t53 = cos(qJ(1));
t48 = t53 * pkin(1);
t34 = t53 * rSges(2,1) - t51 * rSges(2,2);
t33 = -t51 * rSges(2,1) - t53 * rSges(2,2);
t21 = t23 + t48;
t20 = t22 - t68;
t17 = t19 + t48;
t16 = t18 - t68;
t11 = Icges(5,3) * t42 + t55 * t43;
t10 = -Icges(5,3) * t43 + t55 * t42;
t7 = t48 + t9;
t6 = t8 - t68;
t3 = t42 * (t42 * t67 - t65) + t43 * t54;
t1 = [Icges(2,3) + m(2) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2) + t61; m(3) * (t22 * t20 + t23 * t21) + m(4) * (t18 * t16 + t19 * t17) + m(5) * (t8 * t6 + t9 * t7) + t61; m(3) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2) + t61; 0; 0; m(4) + m(5); (-t42 * t7 - t43 * t6) * t70 + t62; (-t42 * t9 - t43 * t8) * t70 + t62; m(5) * t3; m(5) * (t3 ^ 2 + (t73 + t74) * t32 ^ 2) + t42 * (-t10 * t75 + t74 * t11) - t43 * (t73 * t10 - t11 * t75);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
