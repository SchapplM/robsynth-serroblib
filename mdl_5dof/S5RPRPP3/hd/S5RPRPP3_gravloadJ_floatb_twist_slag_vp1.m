% Calculate Gravitation load on the joints for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:20
% EndTime: 2019-12-31 18:12:21
% DurationCPUTime: 0.32s
% Computational Cost: add. (178->65), mult. (200->78), div. (0->0), fcn. (166->6), ass. (0->30)
t13 = pkin(7) + qJ(3);
t11 = sin(t13);
t8 = t11 * qJ(4);
t12 = cos(t13);
t9 = t12 * pkin(3);
t37 = t8 + t9;
t17 = sin(qJ(1));
t18 = cos(qJ(1));
t23 = g(1) * t18 + g(2) * t17;
t41 = t23 * t12;
t29 = rSges(6,3) + qJ(5);
t34 = rSges(6,2) * t11;
t40 = t29 * t12 + t34;
t38 = -m(5) - m(6);
t16 = -pkin(6) - qJ(2);
t33 = rSges(5,1) - t16;
t32 = rSges(4,3) - t16;
t30 = rSges(3,3) + qJ(2);
t28 = rSges(6,1) + pkin(4) - t16;
t27 = -pkin(3) - t29;
t15 = cos(pkin(7));
t10 = pkin(2) * t15 + pkin(1);
t26 = -t10 - t8;
t5 = t18 * t10;
t25 = g(2) * (t37 * t18 + t5);
t24 = qJ(4) * t41;
t22 = rSges(4,1) * t12 - rSges(4,2) * t11;
t21 = -rSges(5,2) * t12 + rSges(5,3) * t11;
t20 = rSges(3,1) * t15 - rSges(3,2) * sin(pkin(7)) + pkin(1);
t1 = [-m(2) * (g(1) * (-t17 * rSges(2,1) - rSges(2,2) * t18) + g(2) * (rSges(2,1) * t18 - t17 * rSges(2,2))) - m(3) * ((g(1) * t30 + g(2) * t20) * t18 + (-g(1) * t20 + g(2) * t30) * t17) - m(4) * (g(2) * t5 + (g(1) * t32 + g(2) * t22) * t18 + (g(1) * (-t10 - t22) + g(2) * t32) * t17) - m(5) * (t25 + (g(1) * t33 + g(2) * t21) * t18 + (g(1) * (-t21 + t26 - t9) + g(2) * t33) * t17) - m(6) * (t25 + (g(1) * t28 + g(2) * t40) * t18 + (g(2) * t28 + (t27 * t12 + t26 - t34) * g(1)) * t17), (-m(3) - m(4) + t38) * (g(1) * t17 - g(2) * t18), -m(4) * g(3) * t22 - m(5) * (g(3) * (t21 + t37) + t24) - m(6) * (g(3) * (t37 + t40) + t24) + t23 * ((m(4) * rSges(4,2) - m(5) * rSges(5,3) - m(6) * rSges(6,2)) * t12 + (m(4) * rSges(4,1) - m(5) * (rSges(5,2) - pkin(3)) - m(6) * t27) * t11), t38 * (-g(3) * t12 + t23 * t11), -m(6) * (g(3) * t11 + t41)];
taug = t1(:);
