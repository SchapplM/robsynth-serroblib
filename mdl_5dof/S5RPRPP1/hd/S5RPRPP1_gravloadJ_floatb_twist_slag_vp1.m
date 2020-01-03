% Calculate Gravitation load on the joints for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:38
% EndTime: 2019-12-31 18:08:39
% DurationCPUTime: 0.25s
% Computational Cost: add. (192->61), mult. (163->72), div. (0->0), fcn. (129->8), ass. (0->30)
t34 = rSges(6,1) + pkin(4);
t11 = qJ(1) + pkin(7);
t5 = sin(t11);
t7 = cos(t11);
t33 = g(1) * t7 + g(2) * t5;
t32 = rSges(6,3) + qJ(5);
t10 = qJ(3) + pkin(8);
t4 = sin(t10);
t6 = cos(t10);
t18 = t32 * t4 + t34 * t6;
t29 = -m(5) - m(6);
t14 = sin(qJ(1));
t28 = pkin(1) * t14;
t13 = sin(qJ(3));
t27 = pkin(3) * t13;
t26 = rSges(4,3) + pkin(6);
t12 = -qJ(4) - pkin(6);
t25 = rSges(6,2) - t12;
t24 = rSges(5,3) - t12;
t23 = g(1) * t28;
t22 = rSges(5,1) * t6 - rSges(5,2) * t4;
t15 = cos(qJ(3));
t21 = rSges(4,1) * t15 - rSges(4,2) * t13;
t20 = pkin(2) + t21;
t8 = t15 * pkin(3);
t3 = t8 + pkin(2);
t16 = cos(qJ(1));
t9 = t16 * pkin(1);
t19 = -t23 + g(2) * (t7 * t3 + t9);
t1 = [-m(2) * (g(1) * (-t14 * rSges(2,1) - rSges(2,2) * t16) + g(2) * (rSges(2,1) * t16 - t14 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t5 - rSges(3,2) * t7 - t28) + g(2) * (rSges(3,1) * t7 - rSges(3,2) * t5 + t9)) - m(4) * (-t23 + g(2) * t9 + (g(1) * t26 + g(2) * t20) * t7 + (-g(1) * t20 + g(2) * t26) * t5) - m(5) * ((g(1) * t24 + g(2) * t22) * t7 + (g(1) * (-t22 - t3) + g(2) * t24) * t5 + t19) - m(6) * ((g(1) * t25 + g(2) * t18) * t7 + (g(1) * (-t18 - t3) + g(2) * t25) * t5 + t19), (-m(3) - m(4) + t29) * g(3), (-m(4) * t21 - m(5) * (t22 + t8) - m(6) * (t8 + t18)) * g(3) + t33 * (-m(4) * (-rSges(4,1) * t13 - rSges(4,2) * t15) - m(5) * (-rSges(5,1) * t4 - rSges(5,2) * t6 - t27) - m(6) * (t32 * t6 - t34 * t4 - t27)), t29 * (g(1) * t5 - g(2) * t7), -m(6) * (-g(3) * t6 + t33 * t4)];
taug = t1(:);
