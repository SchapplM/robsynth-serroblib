% Calculate Gravitation load on the joints for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:05
% EndTime: 2019-12-31 16:49:06
% DurationCPUTime: 0.18s
% Computational Cost: add. (113->40), mult. (106->54), div. (0->0), fcn. (79->8), ass. (0->23)
t13 = cos(qJ(3));
t10 = qJ(3) + qJ(4);
t5 = sin(t10);
t6 = cos(t10);
t22 = t6 * rSges(5,1) - rSges(5,2) * t5;
t30 = t13 * pkin(3) + t22;
t9 = qJ(1) + pkin(7);
t3 = sin(t9);
t4 = cos(t9);
t29 = g(1) * t4 + g(2) * t3;
t12 = sin(qJ(1));
t25 = t12 * pkin(1);
t24 = rSges(4,3) + pkin(5);
t23 = rSges(5,3) + pkin(6) + pkin(5);
t21 = -rSges(5,1) * t5 - rSges(5,2) * t6;
t11 = sin(qJ(3));
t20 = rSges(4,1) * t13 - rSges(4,2) * t11;
t14 = cos(qJ(1));
t8 = t14 * pkin(1);
t19 = -g(1) * t25 + g(2) * t8;
t18 = pkin(2) + t30;
t17 = pkin(2) + t20;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t12 - rSges(2,2) * t14) + g(2) * (rSges(2,1) * t14 - rSges(2,2) * t12)) - m(3) * (g(1) * (-rSges(3,1) * t3 - rSges(3,2) * t4 - t25) + g(2) * (rSges(3,1) * t4 - rSges(3,2) * t3 + t8)) - m(4) * ((g(1) * t24 + g(2) * t17) * t4 + (-g(1) * t17 + g(2) * t24) * t3 + t19) - m(5) * ((g(1) * t23 + g(2) * t18) * t4 + (-g(1) * t18 + g(2) * t23) * t3 + t19), (-m(3) - m(4) - m(5)) * g(3), (-m(4) * t20 - m(5) * t30) * g(3) + t29 * (-m(4) * (-rSges(4,1) * t11 - rSges(4,2) * t13) - m(5) * (-pkin(3) * t11 + t21)), -m(5) * (g(3) * t22 + t29 * t21)];
taug = t1(:);
