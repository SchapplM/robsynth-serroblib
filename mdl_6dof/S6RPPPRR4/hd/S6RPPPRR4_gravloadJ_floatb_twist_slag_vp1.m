% Calculate Gravitation load on the joints for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:18
% EndTime: 2019-03-09 01:35:19
% DurationCPUTime: 0.40s
% Computational Cost: add. (198->75), mult. (361->104), div. (0->0), fcn. (399->8), ass. (0->34)
t34 = sin(pkin(9));
t35 = cos(pkin(9));
t40 = sin(qJ(1));
t41 = cos(qJ(1));
t6 = -t40 * t34 - t41 * t35;
t7 = t41 * t34 - t40 * t35;
t51 = -g(1) * t7 + g(2) * t6;
t42 = rSges(7,3) + pkin(8);
t16 = sin(qJ(6));
t18 = cos(qJ(6));
t27 = rSges(7,1) * t18 - rSges(7,2) * t16 + pkin(5);
t48 = m(6) * rSges(6,1) + m(7) * t27;
t25 = -m(6) * rSges(6,2) + m(7) * t42;
t43 = rSges(6,3) + pkin(7);
t17 = sin(qJ(5));
t39 = t16 * t17;
t38 = t17 * t18;
t37 = t41 * pkin(1) + t40 * qJ(2);
t36 = rSges(5,3) + qJ(4);
t33 = -m(5) - m(6) - m(7);
t32 = t41 * pkin(2) + t37;
t31 = m(4) - t33;
t30 = -t6 * pkin(3) + t32;
t29 = -t40 * pkin(1) + t41 * qJ(2);
t28 = t16 * rSges(7,1) + t18 * rSges(7,2);
t19 = cos(qJ(5));
t26 = -t42 * t19 + qJ(4);
t24 = t17 * rSges(6,1) + t19 * rSges(6,2) + qJ(4);
t23 = -t40 * pkin(2) + t29;
t22 = t7 * pkin(3) + t23;
t20 = g(1) * t22;
t2 = -t6 * t16 + t7 * t38;
t1 = -t6 * t18 - t7 * t39;
t3 = [-m(2) * (g(1) * (-t40 * rSges(2,1) - t41 * rSges(2,2)) + g(2) * (t41 * rSges(2,1) - t40 * rSges(2,2))) - m(3) * (g(1) * (-t40 * rSges(3,1) + t41 * rSges(3,3) + t29) + g(2) * (t41 * rSges(3,1) + t40 * rSges(3,3) + t37)) - m(4) * (g(1) * (t7 * rSges(4,1) - t6 * rSges(4,2) + t23) + g(2) * (-t6 * rSges(4,1) - t7 * rSges(4,2) + t32)) - m(5) * (g(1) * (-t7 * rSges(5,2) + t36 * t6 + t22) + g(2) * (t6 * rSges(5,2) + t36 * t7 + t30)) - m(6) * (t20 + g(2) * t30 + (g(1) * t43 + g(2) * t24) * t7 + (g(1) * t24 - g(2) * t43) * t6) - m(7) * (t20 + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t30) + (g(1) * (pkin(7) + t28) + g(2) * (t17 * pkin(5) + t26)) * t7 + (-g(2) * pkin(7) + (t17 * t27 + t26) * g(1)) * t6) (-m(3) - t31) * (g(1) * t40 - g(2) * t41) t31 * g(3), -t33 * t51 (-t17 * t48 + t25 * t19) * g(3) + t51 * (t25 * t17 + t48 * t19) -m(7) * (g(1) * (t1 * rSges(7,1) - t2 * rSges(7,2)) + g(2) * ((-t7 * t18 + t6 * t39) * rSges(7,1) + (t7 * t16 + t6 * t38) * rSges(7,2)) + g(3) * t28 * t19)];
taug  = t3(:);
