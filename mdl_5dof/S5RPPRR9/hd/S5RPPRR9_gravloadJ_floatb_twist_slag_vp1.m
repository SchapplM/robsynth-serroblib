% Calculate Gravitation load on the joints for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:19
% EndTime: 2019-12-31 18:02:20
% DurationCPUTime: 0.29s
% Computational Cost: add. (154->65), mult. (291->99), div. (0->0), fcn. (319->8), ass. (0->30)
t15 = sin(qJ(5));
t17 = cos(qJ(5));
t42 = m(5) * rSges(5,1) + m(6) * (rSges(6,1) * t17 - rSges(6,2) * t15 + pkin(4));
t18 = cos(qJ(4));
t40 = t18 * pkin(4);
t39 = rSges(5,3) + pkin(6);
t38 = rSges(6,3) + pkin(7);
t37 = cos(qJ(1));
t36 = sin(qJ(1));
t35 = t15 * t18;
t34 = t17 * t18;
t33 = t37 * pkin(1) + t36 * qJ(2);
t32 = cos(pkin(8));
t31 = sin(pkin(8));
t30 = m(4) + m(5) + m(6);
t29 = t37 * pkin(2) + t33;
t5 = -t36 * t31 - t37 * t32;
t28 = -t5 * pkin(3) + t29;
t27 = -t36 * pkin(1) + t37 * qJ(2);
t16 = sin(qJ(4));
t26 = t18 * rSges(5,1) - t16 * rSges(5,2);
t6 = t37 * t31 - t36 * t32;
t25 = t5 * t15 + t6 * t34;
t24 = -t5 * t17 + t6 * t35;
t22 = -m(5) * rSges(5,2) + m(6) * t38;
t21 = -t36 * pkin(2) + t27;
t20 = t6 * pkin(3) + t21;
t2 = t6 * t15 - t5 * t34;
t1 = t6 * t17 + t5 * t35;
t3 = [-m(2) * (g(1) * (-t36 * rSges(2,1) - t37 * rSges(2,2)) + g(2) * (t37 * rSges(2,1) - t36 * rSges(2,2))) - m(3) * (g(1) * (-t36 * rSges(3,1) + t37 * rSges(3,3) + t27) + g(2) * (t37 * rSges(3,1) + t36 * rSges(3,3) + t33)) - m(4) * (g(1) * (t6 * rSges(4,1) - t5 * rSges(4,2) + t21) + g(2) * (-t5 * rSges(4,1) - t6 * rSges(4,2) + t29)) - m(5) * (g(1) * t20 + g(2) * t28 + (g(1) * t26 + g(2) * t39) * t6 + (g(1) * t39 - g(2) * t26) * t5) - m(6) * (g(1) * (t25 * rSges(6,1) - t24 * rSges(6,2) + t5 * pkin(6) + t6 * t40 + t20) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t6 * pkin(6) - t5 * t40 + t28) + (g(1) * t6 - g(2) * t5) * t16 * t38), (-m(3) - t30) * (g(1) * t36 - g(2) * t37), t30 * g(3), (t22 * t16 + t42 * t18) * g(3) + (g(1) * t5 + g(2) * t6) * (-t42 * t16 + t22 * t18), -m(6) * (g(1) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (t24 * rSges(6,1) + t25 * rSges(6,2)) + g(3) * (t15 * rSges(6,1) + t17 * rSges(6,2)) * t16)];
taug = t3(:);
