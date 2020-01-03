% Calculate Gravitation load on the joints for
% S5RRPRP2
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:26
% EndTime: 2019-12-31 19:49:27
% DurationCPUTime: 0.25s
% Computational Cost: add. (268->64), mult. (186->79), div. (0->0), fcn. (147->8), ass. (0->33)
t45 = rSges(6,1) + pkin(4);
t24 = sin(qJ(4));
t26 = cos(qJ(4));
t52 = rSges(5,1) * t26 - rSges(5,2) * t24;
t40 = rSges(6,3) + qJ(5);
t23 = qJ(1) + qJ(2);
t19 = pkin(8) + t23;
t16 = sin(t19);
t17 = cos(t19);
t51 = g(1) * t17 + g(2) * t16;
t50 = t40 * t24 + t45 * t26;
t25 = sin(qJ(1));
t49 = pkin(1) * t25;
t20 = sin(t23);
t48 = pkin(2) * t20;
t42 = t17 * t24;
t41 = t17 * t26;
t21 = cos(t23);
t18 = pkin(2) * t21;
t39 = t17 * pkin(3) + t16 * pkin(7) + t18;
t38 = t17 * pkin(7) - t48;
t37 = t21 * rSges(3,1) - rSges(3,2) * t20;
t36 = t17 * rSges(6,2) + t38;
t35 = t17 * rSges(4,1) - rSges(4,2) * t16 + t18;
t34 = -rSges(3,1) * t20 - rSges(3,2) * t21;
t33 = t16 * rSges(6,2) + t40 * t42 + t45 * t41 + t39;
t32 = -rSges(4,1) * t16 - rSges(4,2) * t17 - t48;
t31 = rSges(5,1) * t41 - rSges(5,2) * t42 + t16 * rSges(5,3) + t39;
t30 = t17 * rSges(5,3) + t38 + (-pkin(3) - t52) * t16;
t29 = g(1) * (-pkin(3) - t50) * t16;
t27 = cos(qJ(1));
t22 = t27 * pkin(1);
t1 = [-m(2) * (g(1) * (-t25 * rSges(2,1) - rSges(2,2) * t27) + g(2) * (rSges(2,1) * t27 - t25 * rSges(2,2))) - m(3) * (g(1) * (t34 - t49) + g(2) * (t22 + t37)) - m(4) * (g(1) * (t32 - t49) + g(2) * (t22 + t35)) - m(5) * (g(1) * (t30 - t49) + g(2) * (t22 + t31)) - m(6) * (g(1) * (t36 - t49) + g(2) * (t22 + t33) + t29), -m(3) * (g(1) * t34 + g(2) * t37) - m(4) * (g(1) * t32 + g(2) * t35) - m(5) * (g(1) * t30 + g(2) * t31) - m(6) * (g(1) * t36 + g(2) * t33 + t29), (-m(4) - m(5) - m(6)) * g(3), (-m(5) * t52 - m(6) * t50) * g(3) + t51 * (-m(5) * (-rSges(5,1) * t24 - rSges(5,2) * t26) - m(6) * (-t45 * t24 + t40 * t26)), -m(6) * (-g(3) * t26 + t51 * t24)];
taug = t1(:);
