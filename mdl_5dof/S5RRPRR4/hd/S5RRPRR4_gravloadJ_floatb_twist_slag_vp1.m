% Calculate Gravitation load on the joints for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:01:48
% EndTime: 2020-01-03 12:01:49
% DurationCPUTime: 0.21s
% Computational Cost: add. (272->70), mult. (177->84), div. (0->0), fcn. (135->10), ass. (0->41)
t63 = rSges(5,3) + pkin(7);
t62 = rSges(6,3) + pkin(8) + pkin(7);
t35 = qJ(4) + qJ(5);
t30 = cos(t35);
t19 = t30 * rSges(6,1);
t28 = sin(t35);
t49 = -rSges(6,2) * t28 + t19;
t37 = sin(qJ(4));
t39 = cos(qJ(4));
t58 = rSges(5,1) * t39;
t61 = -rSges(5,2) * t37 + t58;
t36 = qJ(1) + qJ(2);
t27 = pkin(9) + t36;
t23 = cos(t27);
t54 = t23 * t28;
t55 = rSges(6,2) * t30;
t60 = rSges(6,1) * t54 + t23 * t55;
t22 = sin(t27);
t59 = g(2) * t22;
t53 = t23 * t37;
t29 = sin(t36);
t31 = cos(t36);
t52 = t29 * rSges(3,1) + t31 * rSges(3,2);
t24 = pkin(2) * t29;
t51 = t22 * rSges(4,1) + t23 * rSges(4,2) + t24;
t50 = t31 * rSges(3,1) - rSges(3,2) * t29;
t25 = pkin(2) * t31;
t48 = t23 * rSges(4,1) - rSges(4,2) * t22 + t25;
t47 = rSges(5,1) * t37 + rSges(5,2) * t39;
t46 = -rSges(6,1) * t28 - t55;
t45 = -rSges(5,2) * t53 + t25 + (pkin(3) + t58) * t23 + t63 * t22;
t33 = t39 * pkin(4);
t26 = t33 + pkin(3);
t44 = t24 - t62 * t23 + (t26 + t49) * t22;
t43 = -rSges(6,2) * t54 + t25 + (t19 + t26) * t23 + t62 * t22;
t42 = t24 - t63 * t23 + (pkin(3) + t61) * t22;
t40 = cos(qJ(1));
t38 = sin(qJ(1));
t34 = t40 * pkin(1);
t32 = t38 * pkin(1);
t1 = [-m(2) * (g(2) * (rSges(2,1) * t40 - rSges(2,2) * t38) + g(3) * (rSges(2,1) * t38 + rSges(2,2) * t40)) - m(3) * (g(2) * (t34 + t50) + g(3) * (t32 + t52)) - m(4) * (g(2) * (t34 + t48) + g(3) * (t32 + t51)) - m(5) * (g(2) * (t34 + t45) + g(3) * (t32 + t42)) - m(6) * (g(2) * (t34 + t43) + g(3) * (t32 + t44)), -m(3) * (g(2) * t50 + g(3) * t52) - m(4) * (g(2) * t48 + g(3) * t51) - m(5) * (g(2) * t45 + g(3) * t42) - m(6) * (g(2) * t43 + g(3) * t44), (-m(4) - m(5) - m(6)) * g(1), -m(5) * (g(3) * t47 * t23 + g(1) * t61) - m(6) * (g(1) * (t33 + t49) + g(3) * (pkin(4) * t53 + t60)) + (m(5) * t47 - m(6) * (-pkin(4) * t37 + t46)) * t59, -m(6) * (g(1) * t49 + g(3) * t60 + t46 * t59)];
taug = t1(:);
