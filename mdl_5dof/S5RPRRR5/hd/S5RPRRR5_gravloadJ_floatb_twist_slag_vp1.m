% Calculate Gravitation load on the joints for
% S5RPRRR5
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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:53:47
% EndTime: 2020-01-03 11:53:49
% DurationCPUTime: 0.32s
% Computational Cost: add. (253->65), mult. (165->81), div. (0->0), fcn. (125->10), ass. (0->39)
t60 = rSges(5,3) + pkin(7);
t59 = rSges(6,3) + pkin(8) + pkin(7);
t33 = qJ(4) + qJ(5);
t28 = cos(t33);
t20 = t28 * rSges(6,1);
t27 = sin(t33);
t45 = -rSges(6,2) * t27 + t20;
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t55 = rSges(5,1) * t36;
t58 = -rSges(5,2) * t34 + t55;
t32 = qJ(1) + pkin(9);
t26 = qJ(3) + t32;
t22 = cos(t26);
t51 = t22 * t27;
t52 = rSges(6,2) * t28;
t57 = rSges(6,1) * t51 + t22 * t52;
t21 = sin(t26);
t56 = g(2) * t21;
t50 = t22 * t34;
t49 = t21 * rSges(4,1) + t22 * rSges(4,2);
t24 = sin(t32);
t35 = sin(qJ(1));
t29 = t35 * pkin(1);
t48 = pkin(2) * t24 + t29;
t25 = cos(t32);
t37 = cos(qJ(1));
t31 = t37 * pkin(1);
t47 = pkin(2) * t25 + t31;
t46 = t22 * rSges(4,1) - rSges(4,2) * t21;
t44 = rSges(5,1) * t34 + rSges(5,2) * t36;
t43 = -rSges(6,1) * t27 - t52;
t42 = -rSges(5,2) * t50 + (pkin(3) + t55) * t22 + t60 * t21;
t30 = t36 * pkin(4);
t23 = t30 + pkin(3);
t41 = -t59 * t22 + (t23 + t45) * t21;
t40 = -rSges(6,2) * t51 + (t20 + t23) * t22 + t59 * t21;
t39 = -t60 * t22 + (pkin(3) + t58) * t21;
t1 = [-m(2) * (g(2) * (rSges(2,1) * t37 - rSges(2,2) * t35) + g(3) * (rSges(2,1) * t35 + rSges(2,2) * t37)) - m(3) * (g(2) * (rSges(3,1) * t25 - rSges(3,2) * t24 + t31) + g(3) * (rSges(3,1) * t24 + rSges(3,2) * t25 + t29)) - m(4) * (g(2) * (t46 + t47) + g(3) * (t48 + t49)) - m(5) * (g(2) * (t42 + t47) + g(3) * (t39 + t48)) - m(6) * (g(2) * (t40 + t47) + g(3) * (t41 + t48)), (-m(3) - m(4) - m(5) - m(6)) * g(1), -m(4) * (g(2) * t46 + g(3) * t49) - m(5) * (g(2) * t42 + g(3) * t39) - m(6) * (g(2) * t40 + g(3) * t41), -m(5) * (g(3) * t44 * t22 + g(1) * t58) - m(6) * (g(1) * (t30 + t45) + g(3) * (pkin(4) * t50 + t57)) + (m(5) * t44 - m(6) * (-pkin(4) * t34 + t43)) * t56, -m(6) * (g(1) * t45 + g(3) * t57 + t43 * t56)];
taug = t1(:);
