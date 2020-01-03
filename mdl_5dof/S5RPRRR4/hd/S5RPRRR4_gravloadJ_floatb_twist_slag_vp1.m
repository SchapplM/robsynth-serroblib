% Calculate Gravitation load on the joints for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:02
% EndTime: 2020-01-03 11:52:03
% DurationCPUTime: 0.19s
% Computational Cost: add. (266->54), mult. (140->66), div. (0->0), fcn. (102->10), ass. (0->33)
t49 = rSges(6,3) + pkin(8);
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t48 = rSges(6,1) * t30 - rSges(6,2) * t28;
t47 = pkin(4) + t48;
t27 = qJ(1) + pkin(9);
t24 = qJ(3) + t27;
t21 = qJ(4) + t24;
t15 = sin(t21);
t16 = cos(t21);
t46 = t15 * rSges(5,1) + t16 * rSges(5,2);
t19 = sin(t24);
t20 = cos(t24);
t43 = t19 * rSges(4,1) + t20 * rSges(4,2);
t22 = sin(t27);
t29 = sin(qJ(1));
t25 = t29 * pkin(1);
t42 = pkin(2) * t22 + t25;
t23 = cos(t27);
t31 = cos(qJ(1));
t26 = t31 * pkin(1);
t41 = pkin(2) * t23 + t26;
t13 = pkin(3) * t19;
t40 = t13 + t46;
t39 = t16 * rSges(5,1) - rSges(5,2) * t15;
t38 = t20 * rSges(4,1) - rSges(4,2) * t19;
t14 = pkin(3) * t20;
t37 = t14 + t39;
t35 = t49 * t15 + t47 * t16;
t34 = t14 + t35;
t33 = t47 * t15 - t49 * t16;
t32 = t13 + t33;
t1 = [-m(2) * (g(2) * (rSges(2,1) * t31 - t29 * rSges(2,2)) + g(3) * (t29 * rSges(2,1) + rSges(2,2) * t31)) - m(3) * (g(2) * (rSges(3,1) * t23 - rSges(3,2) * t22 + t26) + g(3) * (rSges(3,1) * t22 + rSges(3,2) * t23 + t25)) - m(4) * (g(2) * (t38 + t41) + g(3) * (t42 + t43)) - m(5) * (g(2) * (t37 + t41) + g(3) * (t40 + t42)) - m(6) * (g(2) * (t34 + t41) + g(3) * (t32 + t42)), (-m(3) - m(4) - m(5) - m(6)) * g(1), -m(4) * (g(2) * t38 + g(3) * t43) - m(5) * (g(2) * t37 + g(3) * t40) - m(6) * (g(2) * t34 + g(3) * t32), -m(5) * (g(2) * t39 + g(3) * t46) - m(6) * (g(2) * t35 + g(3) * t33), -m(6) * (g(1) * t48 + (-g(2) * t15 + g(3) * t16) * (rSges(6,1) * t28 + rSges(6,2) * t30))];
taug = t1(:);
