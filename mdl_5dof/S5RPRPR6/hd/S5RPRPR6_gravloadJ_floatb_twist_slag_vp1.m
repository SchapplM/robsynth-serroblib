% Calculate Gravitation load on the joints for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:42
% EndTime: 2019-12-31 18:17:42
% DurationCPUTime: 0.22s
% Computational Cost: add. (202->54), mult. (126->62), div. (0->0), fcn. (92->8), ass. (0->30)
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t47 = rSges(6,1) * t23 + rSges(6,2) * t25;
t48 = rSges(6,3) + pkin(7);
t22 = qJ(1) + pkin(8);
t20 = qJ(3) + t22;
t17 = cos(t20);
t46 = t47 * t17;
t16 = sin(t20);
t40 = g(1) * t16;
t45 = -g(2) * t17 + t40;
t44 = t17 * rSges(5,3) + (rSges(5,2) - pkin(3)) * t16;
t43 = -m(5) - m(6);
t24 = sin(qJ(1));
t42 = pkin(1) * t24;
t38 = t17 * pkin(3) + t16 * qJ(4);
t19 = cos(t22);
t26 = cos(qJ(1));
t21 = t26 * pkin(1);
t35 = pkin(2) * t19 + t21;
t34 = t17 * rSges(4,1) - rSges(4,2) * t16;
t18 = sin(t22);
t33 = -pkin(2) * t18 - t42;
t32 = (-pkin(3) - t48) * t40;
t31 = -rSges(4,1) * t16 - rSges(4,2) * t17;
t29 = -rSges(5,2) * t17 + t16 * rSges(5,3) + t38;
t28 = t47 * t16 + t48 * t17 + t38;
t7 = t17 * qJ(4);
t27 = t33 + t7;
t1 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - rSges(2,2) * t26) + g(2) * (rSges(2,1) * t26 - t24 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t18 - rSges(3,2) * t19 - t42) + g(2) * (rSges(3,1) * t19 - rSges(3,2) * t18 + t21)) - m(4) * (g(1) * (t31 + t33) + g(2) * (t34 + t35)) - m(5) * (g(1) * (t27 + t44) + g(2) * (t29 + t35)) - m(6) * (g(1) * (t27 + t46) + g(2) * (t28 + t35) + t32), (-m(3) - m(4) + t43) * g(3), -m(4) * (g(1) * t31 + g(2) * t34) - m(5) * (g(1) * (t7 + t44) + g(2) * t29) - m(6) * (g(1) * (t7 + t46) + g(2) * t28 + t32), t43 * t45, -m(6) * (-g(3) * t47 + t45 * (rSges(6,1) * t25 - rSges(6,2) * t23))];
taug = t1(:);
