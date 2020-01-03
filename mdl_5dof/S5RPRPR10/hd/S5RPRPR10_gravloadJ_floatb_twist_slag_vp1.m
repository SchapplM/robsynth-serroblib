% Calculate Gravitation load on the joints for
% S5RPRPR10
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:55
% EndTime: 2019-12-31 18:25:55
% DurationCPUTime: 0.26s
% Computational Cost: add. (198->66), mult. (234->79), div. (0->0), fcn. (240->8), ass. (0->36)
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t36 = qJ(3) + pkin(8);
t33 = sin(t36);
t34 = cos(t36);
t1 = -t22 * t33 - t25 * t34;
t2 = -t22 * t34 + t25 * t33;
t20 = sin(qJ(5));
t23 = cos(qJ(5));
t30 = -rSges(6,1) * t23 + rSges(6,2) * t20;
t26 = pkin(4) - t30;
t42 = rSges(6,3) + pkin(7);
t48 = -(g(1) * t42 - g(2) * t26) * t1 - (g(1) * t26 + g(2) * t42) * t2;
t24 = cos(qJ(3));
t16 = pkin(3) * t24 + pkin(2);
t37 = t25 * pkin(1) + t22 * qJ(2);
t47 = t25 * t16 + t37;
t18 = t25 * qJ(2);
t44 = t18 + (-pkin(1) - t16) * t22;
t43 = m(5) + m(6);
t21 = sin(qJ(3));
t41 = t21 * t25;
t40 = t22 * t21;
t39 = t22 * t24;
t38 = t25 * t24;
t35 = pkin(3) * t38;
t3 = -t38 - t40;
t4 = -t39 + t41;
t32 = -rSges(4,1) * t4 + rSges(4,2) * t3;
t31 = t3 * rSges(4,1) + t4 * rSges(4,2);
t10 = pkin(3) * t40;
t28 = t1 * rSges(5,1) + t2 * rSges(5,2) - t10;
t12 = pkin(3) * t41;
t27 = rSges(5,1) * t2 - rSges(5,2) * t1 + t12;
t11 = pkin(3) * t39;
t5 = [-m(2) * (g(1) * (-t22 * rSges(2,1) - t25 * rSges(2,2)) + g(2) * (rSges(2,1) * t25 - t22 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t25 + t18 + (-rSges(3,1) - pkin(1)) * t22) + g(2) * (rSges(3,1) * t25 + t22 * rSges(3,3) + t37)) - m(4) * (g(1) * (t18 + (-pkin(1) - pkin(2)) * t22 - t32) + g(2) * (pkin(2) * t25 - t31 + t37)) - m(5) * (g(1) * (t27 + t44) + g(2) * (-t28 + t47)) - m(6) * (g(1) * (t12 + t44) + g(2) * (t10 + t47) - t48), (-m(3) - m(4) - t43) * (g(1) * t22 - g(2) * t25), -m(4) * (g(1) * t32 + g(2) * t31) - m(5) * (g(1) * (t11 - t27) + g(2) * (t28 - t35)) - m(6) * (g(1) * (t11 - t12) + g(2) * (-t10 - t35) + t48), t43 * g(3), -m(6) * (g(3) * t30 + (g(1) * t1 + g(2) * t2) * (rSges(6,1) * t20 + rSges(6,2) * t23))];
taug = t5(:);
