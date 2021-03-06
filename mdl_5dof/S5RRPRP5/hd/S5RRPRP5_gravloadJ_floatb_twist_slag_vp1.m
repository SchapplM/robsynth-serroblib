% Calculate Gravitation load on the joints for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:01
% EndTime: 2019-12-31 19:54:03
% DurationCPUTime: 0.40s
% Computational Cost: add. (258->75), mult. (233->92), div. (0->0), fcn. (190->8), ass. (0->43)
t68 = rSges(6,3) + qJ(5);
t67 = rSges(6,1) + pkin(4);
t25 = qJ(2) + pkin(8);
t22 = qJ(4) + t25;
t18 = cos(t22);
t66 = t68 * t18;
t30 = cos(qJ(1));
t65 = g(2) * t30;
t17 = sin(t22);
t64 = t18 * rSges(5,1) - t17 * rSges(5,2);
t20 = sin(t25);
t21 = cos(t25);
t29 = cos(qJ(2));
t23 = t29 * pkin(2);
t63 = t21 * rSges(4,1) - t20 * rSges(4,2) + t23;
t28 = sin(qJ(1));
t62 = g(1) * t30 + g(2) * t28;
t61 = t62 * t17;
t60 = t68 * t17 + t67 * t18;
t46 = pkin(3) * t21 + t23;
t4 = pkin(1) + t46;
t59 = t4 * t65;
t58 = t66 * t28;
t57 = t66 * t30;
t27 = sin(qJ(2));
t54 = t27 * pkin(2);
t52 = rSges(3,3) + pkin(6);
t26 = -qJ(3) - pkin(6);
t24 = -pkin(7) + t26;
t49 = rSges(6,2) - t24;
t48 = rSges(4,3) - t26;
t47 = rSges(5,3) - t24;
t42 = t29 * rSges(3,1) - t27 * rSges(3,2);
t38 = -rSges(5,1) * t17 - rSges(5,2) * t18;
t37 = pkin(1) + t42;
t36 = pkin(1) + t63;
t35 = t38 * t28;
t34 = t38 * t30;
t31 = t67 * t61;
t5 = -pkin(3) * t20 - t54;
t3 = t30 * t5;
t2 = t28 * t5;
t1 = [-m(2) * (g(1) * (-t28 * rSges(2,1) - t30 * rSges(2,2)) + g(2) * (t30 * rSges(2,1) - t28 * rSges(2,2))) - m(3) * ((g(1) * t52 + g(2) * t37) * t30 + (-g(1) * t37 + g(2) * t52) * t28) - m(4) * ((g(1) * t48 + g(2) * t36) * t30 + (-g(1) * t36 + g(2) * t48) * t28) - m(5) * (t59 + (g(1) * t47 + g(2) * t64) * t30 + (g(1) * (-t64 - t4) + g(2) * t47) * t28) - m(6) * (t59 + (g(1) * t49 + g(2) * t60) * t30 + (g(1) * (-t60 - t4) + g(2) * t49) * t28), -m(3) * (g(3) * t42 + t62 * (-rSges(3,1) * t27 - rSges(3,2) * t29)) - m(4) * (g(3) * t63 + t62 * (-rSges(4,1) * t20 - rSges(4,2) * t21 - t54)) - m(5) * (g(1) * (t3 + t34) + g(2) * (t2 + t35) + g(3) * (t64 + t46)) - m(6) * (g(1) * (t3 + t57) + g(2) * (t2 + t58) + g(3) * (t60 + t46) - t31), (-m(4) - m(5) - m(6)) * (g(1) * t28 - t65), -m(5) * (g(1) * t34 + g(2) * t35 + g(3) * t64) - m(6) * (g(1) * t57 + g(2) * t58 + g(3) * t60 - t31), -m(6) * (-g(3) * t18 + t61)];
taug = t1(:);
