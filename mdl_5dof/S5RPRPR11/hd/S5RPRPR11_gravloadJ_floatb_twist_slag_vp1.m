% Calculate Gravitation load on the joints for
% S5RPRPR11
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:11
% EndTime: 2019-12-31 18:27:12
% DurationCPUTime: 0.34s
% Computational Cost: add. (216->77), mult. (252->105), div. (0->0), fcn. (231->8), ass. (0->40)
t20 = pkin(8) + qJ(3);
t18 = sin(t20);
t15 = t18 * qJ(4);
t19 = cos(t20);
t41 = t19 * pkin(3) + t15;
t25 = sin(qJ(1));
t27 = cos(qJ(1));
t50 = g(1) * t27 + g(2) * t25;
t49 = -m(5) - m(6);
t46 = t19 * pkin(4);
t26 = cos(qJ(5));
t45 = t18 * t26;
t44 = t19 * t27;
t23 = -pkin(6) - qJ(2);
t43 = rSges(5,2) - t23;
t42 = rSges(4,3) - t23;
t40 = qJ(4) * t19;
t39 = rSges(3,3) + qJ(2);
t38 = -rSges(6,3) - pkin(7) - t23;
t22 = cos(pkin(8));
t17 = t22 * pkin(2) + pkin(1);
t12 = t27 * t17;
t37 = pkin(3) * t44 + t27 * t15 + t12;
t24 = sin(qJ(5));
t7 = -t19 * t24 + t45;
t2 = t7 * t25;
t31 = t18 * t24 + t19 * t26;
t3 = t31 * t25;
t36 = t2 * rSges(6,1) - t3 * rSges(6,2);
t4 = t24 * t44 - t27 * t45;
t5 = t31 * t27;
t35 = -t4 * rSges(6,1) - t5 * rSges(6,2);
t34 = -rSges(6,1) * t31 - t7 * rSges(6,2);
t33 = t19 * rSges(4,1) - t18 * rSges(4,2);
t32 = t19 * rSges(5,1) + t18 * rSges(5,3);
t30 = rSges(3,1) * t22 - rSges(3,2) * sin(pkin(8)) + pkin(1);
t29 = -t17 - t41;
t11 = t27 * t40;
t9 = t25 * t40;
t1 = [-m(2) * (g(1) * (-t25 * rSges(2,1) - t27 * rSges(2,2)) + g(2) * (t27 * rSges(2,1) - t25 * rSges(2,2))) - m(3) * ((g(1) * t39 + g(2) * t30) * t27 + (-g(1) * t30 + g(2) * t39) * t25) - m(4) * (g(2) * t12 + (g(1) * t42 + g(2) * t33) * t27 + (g(1) * (-t17 - t33) + g(2) * t42) * t25) - m(5) * (g(2) * t37 + (g(1) * t43 + g(2) * t32) * t27 + (g(1) * (t29 - t32) + g(2) * t43) * t25) - m(6) * (g(1) * (-t3 * rSges(6,1) - t2 * rSges(6,2)) + g(2) * (t5 * rSges(6,1) - t4 * rSges(6,2) + t37) + (g(1) * t38 + g(2) * t46) * t27 + (g(1) * (t29 - t46) + g(2) * t38) * t25), (-m(3) - m(4) + t49) * (g(1) * t25 - g(2) * t27), -m(4) * g(3) * t33 - m(5) * (g(1) * t11 + g(2) * t9 + g(3) * (t32 + t41)) - m(6) * (g(1) * (t11 - t35) + g(2) * (-t36 + t9) + g(3) * (-t34 + t41 + t46)) + t50 * ((m(4) * rSges(4,2) - m(5) * rSges(5,3)) * t19 + (m(4) * rSges(4,1) - m(5) * (-rSges(5,1) - pkin(3)) - m(6) * (-pkin(3) - pkin(4))) * t18), t49 * (-g(3) * t19 + t50 * t18), -m(6) * (g(1) * t35 + g(2) * t36 + g(3) * t34)];
taug = t1(:);
