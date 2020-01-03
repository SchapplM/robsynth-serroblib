% Calculate Gravitation load on the joints for
% S5RRPRR3
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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:21
% EndTime: 2020-01-03 12:00:22
% DurationCPUTime: 0.19s
% Computational Cost: add. (285->57), mult. (152->69), div. (0->0), fcn. (112->10), ass. (0->35)
t54 = rSges(6,3) + pkin(8);
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t53 = rSges(6,1) * t33 - rSges(6,2) * t31;
t52 = pkin(4) + t53;
t30 = qJ(1) + qJ(2);
t25 = pkin(9) + t30;
t24 = qJ(4) + t25;
t15 = sin(t24);
t16 = cos(t24);
t51 = t15 * rSges(5,1) + t16 * rSges(5,2);
t20 = sin(t25);
t26 = sin(t30);
t22 = pkin(2) * t26;
t48 = pkin(3) * t20 + t22;
t21 = cos(t25);
t27 = cos(t30);
t23 = pkin(2) * t27;
t47 = pkin(3) * t21 + t23;
t46 = t26 * rSges(3,1) + t27 * rSges(3,2);
t45 = t20 * rSges(4,1) + t21 * rSges(4,2) + t22;
t44 = t16 * rSges(5,1) - rSges(5,2) * t15;
t43 = t48 + t51;
t42 = t27 * rSges(3,1) - rSges(3,2) * t26;
t41 = t21 * rSges(4,1) - rSges(4,2) * t20 + t23;
t39 = t44 + t47;
t38 = t54 * t15 + t52 * t16;
t37 = t38 + t47;
t36 = t52 * t15 - t54 * t16;
t35 = t36 + t48;
t34 = cos(qJ(1));
t32 = sin(qJ(1));
t29 = t34 * pkin(1);
t28 = t32 * pkin(1);
t1 = [-m(2) * (g(2) * (rSges(2,1) * t34 - t32 * rSges(2,2)) + g(3) * (t32 * rSges(2,1) + rSges(2,2) * t34)) - m(3) * (g(2) * (t29 + t42) + g(3) * (t28 + t46)) - m(4) * (g(2) * (t29 + t41) + g(3) * (t28 + t45)) - m(5) * (g(2) * (t29 + t39) + g(3) * (t28 + t43)) - m(6) * (g(2) * (t29 + t37) + g(3) * (t28 + t35)), -m(3) * (g(2) * t42 + g(3) * t46) - m(4) * (g(2) * t41 + g(3) * t45) - m(5) * (g(2) * t39 + g(3) * t43) - m(6) * (g(2) * t37 + g(3) * t35), (-m(4) - m(5) - m(6)) * g(1), -m(5) * (g(2) * t44 + g(3) * t51) - m(6) * (g(2) * t38 + g(3) * t36), -m(6) * (g(1) * t53 + (-g(2) * t15 + g(3) * t16) * (rSges(6,1) * t31 + rSges(6,2) * t33))];
taug = t1(:);
