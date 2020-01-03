% Calculate Gravitation load on the joints for
% S4RPRR6
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:29
% EndTime: 2019-12-31 16:52:29
% DurationCPUTime: 0.20s
% Computational Cost: add. (117->41), mult. (120->55), div. (0->0), fcn. (93->8), ass. (0->23)
t13 = pkin(7) + qJ(3);
t10 = qJ(4) + t13;
t5 = sin(t10);
t6 = cos(t10);
t25 = t6 * rSges(5,1) - rSges(5,2) * t5;
t9 = cos(t13);
t33 = pkin(3) * t9 + t25;
t17 = sin(qJ(1));
t18 = cos(qJ(1));
t32 = g(1) * t18 + g(2) * t17;
t15 = cos(pkin(7));
t7 = t15 * pkin(2) + pkin(1);
t16 = -pkin(5) - qJ(2);
t28 = rSges(4,3) - t16;
t27 = rSges(5,3) + pkin(6) - t16;
t26 = rSges(3,3) + qJ(2);
t8 = sin(t13);
t24 = rSges(4,1) * t9 - rSges(4,2) * t8;
t23 = -rSges(5,1) * t5 - rSges(5,2) * t6;
t22 = t24 + t7;
t21 = t7 + t33;
t20 = rSges(3,1) * t15 - rSges(3,2) * sin(pkin(7)) + pkin(1);
t1 = [-m(2) * (g(1) * (-t17 * rSges(2,1) - rSges(2,2) * t18) + g(2) * (rSges(2,1) * t18 - t17 * rSges(2,2))) - m(3) * ((g(1) * t26 + g(2) * t20) * t18 + (-g(1) * t20 + g(2) * t26) * t17) - m(4) * ((g(1) * t28 + g(2) * t22) * t18 + (-g(1) * t22 + g(2) * t28) * t17) - m(5) * ((g(1) * t27 + g(2) * t21) * t18 + (-g(1) * t21 + g(2) * t27) * t17), (-m(3) - m(4) - m(5)) * (g(1) * t17 - g(2) * t18), (-m(4) * t24 - m(5) * t33) * g(3) + t32 * (-m(4) * (-rSges(4,1) * t8 - rSges(4,2) * t9) - m(5) * (-pkin(3) * t8 + t23)), -m(5) * (g(3) * t25 + t32 * t23)];
taug = t1(:);
