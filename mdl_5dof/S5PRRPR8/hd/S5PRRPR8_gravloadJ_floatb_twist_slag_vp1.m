% Calculate Gravitation load on the joints for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:19
% EndTime: 2019-12-31 17:42:20
% DurationCPUTime: 0.31s
% Computational Cost: add. (216->59), mult. (206->80), div. (0->0), fcn. (171->10), ass. (0->38)
t24 = qJ(2) + qJ(3);
t20 = pkin(9) + t24;
t17 = sin(t20);
t18 = cos(t20);
t27 = sin(qJ(5));
t44 = t27 * rSges(6,2);
t60 = rSges(6,3) + pkin(7);
t61 = t17 * t44 + t18 * t60;
t29 = cos(qJ(5));
t59 = -t29 * rSges(6,1) - pkin(4);
t25 = sin(pkin(8));
t26 = cos(pkin(8));
t57 = g(1) * t26 + g(2) * t25;
t56 = t59 * t17;
t55 = -m(5) - m(6);
t28 = sin(qJ(2));
t54 = pkin(2) * t28;
t21 = sin(t24);
t53 = pkin(3) * t21;
t48 = t25 * t27;
t47 = t25 * t29;
t46 = t26 * t27;
t45 = t26 * t29;
t42 = t61 * t25;
t41 = t61 * t26;
t22 = cos(t24);
t38 = t22 * rSges(4,1) - rSges(4,2) * t21;
t19 = pkin(3) * t22;
t37 = t18 * rSges(5,1) - rSges(5,2) * t17 + t19;
t35 = -rSges(4,1) * t21 - rSges(4,2) * t22;
t34 = -rSges(5,1) * t17 - rSges(5,2) * t18;
t32 = t19 + t60 * t17 + (-t44 - t59) * t18;
t30 = cos(qJ(2));
t23 = t30 * pkin(2);
t5 = -t53 - t54;
t2 = t26 * t5;
t1 = t25 * t5;
t3 = [(-m(2) - m(3) - m(4) + t55) * g(3), -m(3) * (g(3) * (rSges(3,1) * t30 - t28 * rSges(3,2)) + t57 * (-rSges(3,1) * t28 - rSges(3,2) * t30)) - m(4) * (g(3) * (t23 + t38) + t57 * (t35 - t54)) - m(5) * (g(1) * (t34 * t26 + t2) + g(2) * (t34 * t25 + t1) + g(3) * (t23 + t37)) - m(6) * (g(1) * (t2 + t41) + g(2) * (t1 + t42) + g(3) * (t23 + t32) + t57 * t56), -m(6) * (g(1) * t41 + g(2) * t42) + (-m(4) * t38 - m(5) * t37 - m(6) * t32) * g(3) + t57 * (-m(4) * t35 - m(5) * (t34 - t53) - m(6) * (-t53 + t56)), t55 * (g(1) * t25 - g(2) * t26), -m(6) * (g(1) * ((-t18 * t46 + t47) * rSges(6,1) + (-t18 * t45 - t48) * rSges(6,2)) + g(2) * ((-t18 * t48 - t45) * rSges(6,1) + (-t18 * t47 + t46) * rSges(6,2)) + g(3) * (-rSges(6,1) * t27 - rSges(6,2) * t29) * t17)];
taug = t3(:);
