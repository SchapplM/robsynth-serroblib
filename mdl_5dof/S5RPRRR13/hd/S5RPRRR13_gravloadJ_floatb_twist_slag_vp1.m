% Calculate Gravitation load on the joints for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR13_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR13_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:33
% EndTime: 2019-12-31 19:14:34
% DurationCPUTime: 0.50s
% Computational Cost: add. (185->88), mult. (300->120), div. (0->0), fcn. (277->8), ass. (0->44)
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t64 = -g(1) * t24 + g(2) * t27;
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t48 = rSges(5,3) + pkin(7);
t62 = t23 * pkin(3) - t48 * t26;
t37 = rSges(6,3) + pkin(8) + pkin(7);
t60 = t23 * rSges(4,1) + t26 * rSges(4,2);
t25 = cos(qJ(4));
t14 = t25 * pkin(4) + pkin(3);
t21 = qJ(4) + qJ(5);
t15 = sin(t21);
t16 = cos(t21);
t22 = sin(qJ(4));
t29 = m(5) * (rSges(5,1) * t25 - rSges(5,2) * t22 + pkin(3)) + m(6) * (rSges(6,1) * t16 - rSges(6,2) * t15 + t14) + m(4) * rSges(4,1);
t59 = -m(4) * rSges(4,2) + m(5) * t48 + m(6) * t37;
t58 = -pkin(1) - pkin(6);
t40 = t27 * t16;
t46 = t24 * t15;
t5 = -t23 * t46 + t40;
t41 = t27 * t15;
t45 = t24 * t16;
t6 = t23 * t45 + t41;
t57 = t5 * rSges(6,1) - t6 * rSges(6,2);
t7 = t23 * t41 + t45;
t8 = t23 * t40 - t46;
t56 = t7 * rSges(6,1) + t8 * rSges(6,2);
t53 = pkin(4) * t22;
t50 = g(3) * t26;
t44 = t24 * t22;
t43 = t24 * t25;
t39 = t27 * t22;
t38 = t27 * t25;
t36 = t27 * pkin(1) + t24 * qJ(2);
t35 = t27 * pkin(6) + t36;
t34 = -rSges(6,1) * t15 - rSges(6,2) * t16;
t11 = t23 * t39 + t43;
t9 = -t23 * t44 + t38;
t31 = t23 * t14 - t37 * t26;
t18 = t27 * qJ(2);
t12 = t23 * t38 - t44;
t10 = t23 * t43 + t39;
t1 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - t27 * rSges(2,2)) + g(2) * (t27 * rSges(2,1) - t24 * rSges(2,2))) - m(3) * (g(1) * (t27 * rSges(3,3) + t18 + (rSges(3,2) - pkin(1)) * t24) + g(2) * (-t27 * rSges(3,2) + t24 * rSges(3,3) + t36)) - m(4) * (g(1) * (t60 * t27 + t18) + g(2) * (t27 * rSges(4,3) + t35) + (g(1) * (-rSges(4,3) + t58) + g(2) * t60) * t24) - m(5) * ((t10 * rSges(5,1) + t9 * rSges(5,2) + t62 * t24 + t35) * g(2) + (t12 * rSges(5,1) - t11 * rSges(5,2) + t58 * t24 + t62 * t27 + t18) * g(1)) - m(6) * (g(1) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t18) + g(2) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t35) + (g(1) * t31 + g(2) * t53) * t27 + (g(1) * (-t53 + t58) + g(2) * t31) * t24), -(-m(3) - m(4) - m(5) - m(6)) * t64, (t29 * t23 - t26 * t59) * g(3) + t64 * (t59 * t23 + t29 * t26), -m(5) * (g(1) * (t9 * rSges(5,1) - t10 * rSges(5,2)) + g(2) * (t11 * rSges(5,1) + t12 * rSges(5,2))) - m(6) * (g(1) * (t9 * pkin(4) + t57) + g(2) * (t11 * pkin(4) + t56)) + (-m(5) * (-rSges(5,1) * t22 - rSges(5,2) * t25) - m(6) * (t34 - t53)) * t50, -m(6) * (g(1) * t57 + g(2) * t56 + t34 * t50)];
taug = t1(:);
