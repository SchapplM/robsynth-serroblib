% Calculate Gravitation load on the joints for
% S5RPRPR1
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:03
% EndTime: 2019-12-05 17:47:04
% DurationCPUTime: 0.25s
% Computational Cost: add. (148->67), mult. (172->83), div. (0->0), fcn. (135->8), ass. (0->33)
t16 = qJ(3) + pkin(8);
t11 = qJ(5) + t16;
t7 = sin(t11);
t8 = cos(t11);
t28 = rSges(6,1) * t7 + rSges(6,2) * t8;
t18 = sin(qJ(3));
t36 = pkin(3) * t18;
t9 = sin(t16);
t25 = pkin(4) * t9 + t28 + t36;
t10 = cos(t16);
t20 = cos(qJ(3));
t35 = pkin(3) * t20;
t40 = -m(4) * (rSges(4,1) * t20 - rSges(4,2) * t18) - m(5) * (rSges(5,1) * t10 - rSges(5,2) * t9 + t35);
t39 = -m(5) - m(6);
t38 = rSges(6,1) * t8;
t37 = rSges(6,2) * t7;
t19 = sin(qJ(1));
t34 = g(1) * t19;
t21 = cos(qJ(1));
t33 = g(2) * t21;
t32 = rSges(4,3) + pkin(6);
t17 = -qJ(4) - pkin(6);
t31 = rSges(5,3) - t17;
t30 = rSges(6,3) + pkin(7) - t17;
t29 = t21 * pkin(1) + t19 * qJ(2);
t26 = rSges(4,1) * t18 + rSges(4,2) * t20;
t13 = t21 * qJ(2);
t24 = g(1) * t13 + g(2) * t29;
t23 = rSges(5,1) * t9 + rSges(5,2) * t10 + t36;
t4 = t21 * t37;
t3 = t19 * t38;
t2 = pkin(4) * t10 + t35;
t1 = [-m(2) * (g(1) * (-t19 * rSges(2,1) - rSges(2,2) * t21) + g(2) * (rSges(2,1) * t21 - t19 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t21 + t13 + (rSges(3,2) - pkin(1)) * t19) + g(2) * (-rSges(3,2) * t21 + t19 * rSges(3,3) + t29)) - m(4) * ((g(1) * t26 + g(2) * t32) * t21 + (g(1) * (-pkin(1) - t32) + g(2) * t26) * t19 + t24) - m(5) * ((g(1) * t23 + g(2) * t31) * t21 + (g(1) * (-pkin(1) - t31) + g(2) * t23) * t19 + t24) - m(6) * ((g(1) * t25 + g(2) * t30) * t21 + (g(1) * (-pkin(1) - t30) + g(2) * t25) * t19 + t24), (-m(3) - m(4) + t39) * (-t33 + t34), -m(6) * (g(1) * t3 + g(2) * t4) + (m(4) * t26 + m(5) * t23 + m(6) * t25) * g(3) + (-m(6) * (-t2 - t38) - t40) * t33 + (-m(6) * (t2 - t37) + t40) * t34, t39 * (g(1) * t21 + g(2) * t19), -m(6) * (g(1) * (-t19 * t37 + t3) + g(2) * (-t21 * t38 + t4) - g(3) * t28)];
taug = t1(:);
