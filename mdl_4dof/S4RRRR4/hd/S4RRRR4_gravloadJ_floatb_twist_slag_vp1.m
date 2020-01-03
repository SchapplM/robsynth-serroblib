% Calculate Gravitation load on the joints for
% S4RRRR4
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
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:41
% EndTime: 2019-12-31 17:25:42
% DurationCPUTime: 0.38s
% Computational Cost: add. (173->64), mult. (226->91), div. (0->0), fcn. (198->8), ass. (0->37)
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t55 = g(1) * t27 + g(2) * t24;
t59 = rSges(5,3) + pkin(7);
t21 = qJ(2) + qJ(3);
t18 = sin(t21);
t19 = cos(t21);
t56 = t19 * rSges(4,1) - rSges(4,2) * t18;
t31 = t19 * pkin(3) + t59 * t18;
t25 = cos(qJ(4));
t49 = rSges(5,1) * t25;
t54 = (-pkin(3) - t49) * t18;
t23 = sin(qJ(2));
t53 = pkin(2) * t23;
t50 = rSges(3,3) + pkin(5);
t22 = sin(qJ(4));
t47 = rSges(5,2) * t22;
t44 = t22 * t24;
t43 = t22 * t27;
t42 = t24 * t25;
t41 = t25 * t27;
t28 = -pkin(6) - pkin(5);
t40 = rSges(4,3) - t28;
t26 = cos(qJ(2));
t36 = rSges(3,1) * t26 - rSges(3,2) * t23;
t34 = -rSges(4,1) * t18 - rSges(4,2) * t19;
t33 = pkin(1) + t36;
t32 = t31 + (-t47 + t49) * t19;
t30 = t55 * (t18 * t47 + t19 * t59);
t20 = t26 * pkin(2);
t17 = t20 + pkin(1);
t7 = t27 * t17;
t4 = t19 * t41 + t44;
t3 = -t19 * t43 + t42;
t2 = -t19 * t42 + t43;
t1 = t19 * t44 + t41;
t5 = [-m(2) * (g(1) * (-rSges(2,1) * t24 - rSges(2,2) * t27) + g(2) * (rSges(2,1) * t27 - rSges(2,2) * t24)) - m(3) * ((g(1) * t50 + g(2) * t33) * t27 + (-g(1) * t33 + g(2) * t50) * t24) - m(4) * (g(2) * t7 + (g(1) * t40 + g(2) * t56) * t27 + (g(1) * (-t17 - t56) + g(2) * t40) * t24) - m(5) * (g(1) * (t2 * rSges(5,1) + t1 * rSges(5,2)) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t7) + (-g(1) * t28 + g(2) * t31) * t27 + (g(1) * (-t17 - t31) - g(2) * t28) * t24), -m(5) * t30 + (-m(3) * t36 - m(4) * (t20 + t56) - m(5) * (t20 + t32)) * g(3) + t55 * (-m(3) * (-rSges(3,1) * t23 - rSges(3,2) * t26) - m(4) * (t34 - t53) - m(5) * (-t53 + t54)), -m(4) * (g(3) * t56 + t55 * t34) - m(5) * (g(3) * t32 + t55 * t54 + t30), -m(5) * (g(1) * (rSges(5,1) * t3 - rSges(5,2) * t4) + g(2) * (-rSges(5,1) * t1 + rSges(5,2) * t2) + g(3) * (-rSges(5,1) * t22 - rSges(5,2) * t25) * t18)];
taug = t5(:);
