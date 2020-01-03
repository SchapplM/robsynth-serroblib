% Calculate Gravitation load on the joints for
% S4RPRR4
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:14
% EndTime: 2019-12-31 16:50:15
% DurationCPUTime: 0.24s
% Computational Cost: add. (127->52), mult. (148->76), div. (0->0), fcn. (129->8), ass. (0->25)
t13 = sin(qJ(3));
t16 = cos(qJ(3));
t33 = rSges(4,1) * t16 - rSges(4,2) * t13;
t12 = sin(qJ(4));
t15 = cos(qJ(4));
t32 = m(4) * rSges(4,1) + m(5) * (rSges(5,1) * t15 - rSges(5,2) * t12 + pkin(3));
t27 = rSges(5,3) + pkin(6);
t31 = pkin(3) * t16 + t27 * t13;
t14 = sin(qJ(1));
t29 = pkin(1) * t14;
t24 = t12 * t16;
t23 = t15 * t16;
t17 = cos(qJ(1));
t10 = t17 * pkin(1);
t11 = qJ(1) + pkin(7);
t8 = sin(t11);
t9 = cos(t11);
t22 = t9 * pkin(2) + t8 * pkin(5) + t10;
t21 = t9 * pkin(5) - t29;
t19 = m(4) * rSges(4,2) - m(5) * t27;
t4 = t12 * t8 + t9 * t23;
t3 = t15 * t8 - t9 * t24;
t2 = t12 * t9 - t8 * t23;
t1 = t15 * t9 + t8 * t24;
t5 = [-m(2) * (g(1) * (-t14 * rSges(2,1) - rSges(2,2) * t17) + g(2) * (rSges(2,1) * t17 - t14 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t8 - rSges(3,2) * t9 - t29) + g(2) * (rSges(3,1) * t9 - rSges(3,2) * t8 + t10)) - m(4) * (g(1) * (rSges(4,3) * t9 + t21) + g(2) * (t33 * t9 + t22) + (g(1) * (-pkin(2) - t33) + g(2) * rSges(4,3)) * t8) - m(5) * ((rSges(5,1) * t4 + rSges(5,2) * t3 + t31 * t9 + t22) * g(2) + (rSges(5,1) * t2 + rSges(5,2) * t1 + t21 + (-pkin(2) - t31) * t8) * g(1)), (-m(3) - m(4) - m(5)) * g(3), (t19 * t13 - t32 * t16) * g(3) + (g(1) * t9 + g(2) * t8) * (t32 * t13 + t19 * t16), -m(5) * (g(1) * (rSges(5,1) * t3 - rSges(5,2) * t4) + g(2) * (-rSges(5,1) * t1 + rSges(5,2) * t2) + g(3) * (-rSges(5,1) * t12 - rSges(5,2) * t15) * t13)];
taug = t5(:);
