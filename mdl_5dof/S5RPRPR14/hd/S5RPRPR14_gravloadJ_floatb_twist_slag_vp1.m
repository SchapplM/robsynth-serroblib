% Calculate Gravitation load on the joints for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR14_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:31
% EndTime: 2019-12-31 18:34:32
% DurationCPUTime: 0.36s
% Computational Cost: add. (156->81), mult. (216->108), div. (0->0), fcn. (187->8), ass. (0->38)
t16 = qJ(3) + pkin(8);
t12 = cos(t16);
t19 = sin(qJ(3));
t22 = cos(qJ(3));
t18 = sin(qJ(5));
t21 = cos(qJ(5));
t24 = rSges(6,1) * t21 - rSges(6,2) * t18 + pkin(4);
t51 = -(m(5) * rSges(5,1) + m(6) * t24) * t12 - m(4) * (rSges(4,1) * t22 - rSges(4,2) * t19);
t20 = sin(qJ(1));
t40 = g(1) * t20;
t35 = rSges(6,3) + pkin(7);
t49 = t35 * t12;
t23 = cos(qJ(1));
t48 = g(1) * t23 + g(2) * t20;
t42 = pkin(3) * t22;
t46 = t42 * t40;
t45 = -m(5) - m(6);
t43 = pkin(3) * t19;
t11 = sin(t16);
t41 = pkin(4) * t11;
t37 = g(2) * t23;
t36 = rSges(4,3) + pkin(6);
t34 = t18 * t23;
t33 = t20 * t18;
t32 = t20 * t21;
t31 = t21 * t23;
t30 = t23 * pkin(1) + t20 * qJ(2);
t14 = t23 * qJ(2);
t17 = -qJ(4) - pkin(6);
t29 = t20 * t17 + t23 * t43 + t14;
t28 = t20 * t43 + t30;
t26 = rSges(4,1) * t19 + rSges(4,2) * t22;
t25 = rSges(5,1) * t11 + rSges(5,2) * t12;
t4 = t11 * t31 - t33;
t3 = t11 * t34 + t32;
t2 = t11 * t32 + t34;
t1 = -t11 * t33 + t31;
t5 = [-m(2) * (g(1) * (-t20 * rSges(2,1) - rSges(2,2) * t23) + g(2) * (rSges(2,1) * t23 - t20 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t23 + t14 + (rSges(3,2) - pkin(1)) * t20) + g(2) * (-rSges(3,2) * t23 + t20 * rSges(3,3) + t30)) - m(4) * (g(1) * t14 + g(2) * t30 + (g(1) * t26 + g(2) * t36) * t23 + (g(1) * (-pkin(1) - t36) + g(2) * t26) * t20) - m(5) * (g(1) * t29 + g(2) * t28 + (g(1) * t25 + g(2) * (rSges(5,3) - t17)) * t23 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t25) * t20) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) - t20 * pkin(1) + t23 * t41 + t29) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) - t17 * t23 + t20 * t41 + t28) - t48 * t49), (-m(3) - m(4) + t45) * (-t37 + t40), m(4) * g(3) * t26 - m(5) * (t46 + g(3) * (-t25 - t43)) - m(6) * (t46 + g(3) * (-t24 * t11 - t43 + t49)) + ((m(5) * rSges(5,2) - m(6) * t35) * t11 + t51) * t40 + (-m(5) * (rSges(5,2) * t11 - t42) - m(6) * (-t35 * t11 - t42) - t51) * t37, t45 * t48, -m(6) * (g(1) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (rSges(6,1) * t3 + rSges(6,2) * t4) + g(3) * (-rSges(6,1) * t18 - rSges(6,2) * t21) * t12)];
taug = t5(:);
