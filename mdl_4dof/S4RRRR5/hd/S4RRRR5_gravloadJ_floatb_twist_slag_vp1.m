% Calculate Gravitation load on the joints for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:22
% EndTime: 2019-12-31 17:27:24
% DurationCPUTime: 0.36s
% Computational Cost: add. (164->75), mult. (271->111), div. (0->0), fcn. (255->8), ass. (0->38)
t51 = rSges(5,3) + pkin(7) + pkin(6);
t21 = sin(qJ(2));
t24 = cos(qJ(2));
t50 = t24 * rSges(3,1) - t21 * rSges(3,2);
t41 = rSges(4,3) + pkin(6);
t49 = t24 * pkin(2) + t21 * t41;
t23 = cos(qJ(3));
t13 = pkin(3) * t23 + pkin(2);
t19 = qJ(3) + qJ(4);
t14 = sin(t19);
t15 = cos(t19);
t20 = sin(qJ(3));
t48 = m(3) * rSges(3,1) + m(4) * (rSges(4,1) * t23 - rSges(4,2) * t20 + pkin(2)) + m(5) * (rSges(5,1) * t15 - rSges(5,2) * t14 + t13);
t25 = cos(qJ(1));
t36 = t25 * t15;
t22 = sin(qJ(1));
t39 = t22 * t24;
t5 = t14 * t39 + t36;
t37 = t25 * t14;
t6 = -t15 * t39 + t37;
t47 = -t5 * rSges(5,1) + t6 * rSges(5,2);
t7 = t15 * t22 - t24 * t37;
t8 = t14 * t22 + t24 * t36;
t46 = t7 * rSges(5,1) - t8 * rSges(5,2);
t44 = pkin(3) * t20;
t43 = g(3) * t21;
t35 = t25 * t20;
t34 = t25 * t23;
t33 = t25 * pkin(1) + t22 * pkin(5);
t32 = -rSges(5,1) * t14 - rSges(5,2) * t15;
t11 = t22 * t23 - t24 * t35;
t9 = t20 * t39 + t34;
t29 = t24 * t13 + t21 * t51;
t28 = m(3) * rSges(3,2) - m(4) * t41 - m(5) * t51;
t17 = t25 * pkin(5);
t12 = t20 * t22 + t24 * t34;
t10 = -t23 * t39 + t35;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t22 - rSges(2,2) * t25) + g(2) * (rSges(2,1) * t25 - rSges(2,2) * t22)) - m(3) * (g(1) * (t25 * rSges(3,3) + t17) + g(2) * (t25 * t50 + t33) + (g(1) * (-pkin(1) - t50) + g(2) * rSges(3,3)) * t22) - m(4) * ((t12 * rSges(4,1) + t11 * rSges(4,2) + t25 * t49 + t33) * g(2) + (t10 * rSges(4,1) + t9 * rSges(4,2) + t17 + (-pkin(1) - t49) * t22) * g(1)) - m(5) * (g(1) * (t6 * rSges(5,1) + t5 * rSges(5,2) + t17) + g(2) * (t8 * rSges(5,1) + t7 * rSges(5,2) + t33) + (g(1) * t44 + g(2) * t29) * t25 + (g(1) * (-pkin(1) - t29) + g(2) * t44) * t22), (t28 * t21 - t24 * t48) * g(3) + (g(1) * t25 + g(2) * t22) * (t21 * t48 + t28 * t24), -m(4) * (g(1) * (rSges(4,1) * t11 - rSges(4,2) * t12) + g(2) * (-rSges(4,1) * t9 + rSges(4,2) * t10)) - m(5) * (g(1) * (pkin(3) * t11 + t46) + g(2) * (-pkin(3) * t9 + t47)) + (-m(4) * (-rSges(4,1) * t20 - rSges(4,2) * t23) - m(5) * (t32 - t44)) * t43, -m(5) * (g(1) * t46 + g(2) * t47 + t32 * t43)];
taug = t1(:);
