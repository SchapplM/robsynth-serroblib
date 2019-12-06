% Calculate Gravitation load on the joints for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:23:58
% EndTime: 2019-12-05 15:24:00
% DurationCPUTime: 0.26s
% Computational Cost: add. (144->47), mult. (166->65), div. (0->0), fcn. (140->10), ass. (0->23)
t12 = sin(pkin(7));
t14 = cos(pkin(7));
t33 = g(1) * t14 + g(2) * t12;
t32 = -m(5) - m(6);
t16 = sin(qJ(2));
t31 = pkin(2) * t16;
t10 = qJ(2) + pkin(8);
t7 = cos(t10);
t28 = t12 * t7;
t27 = t14 * t7;
t26 = rSges(6,3) + pkin(6) + qJ(4);
t25 = rSges(5,3) + qJ(4);
t24 = -m(4) + t32;
t13 = cos(pkin(9));
t9 = pkin(9) + qJ(5);
t4 = sin(t9);
t6 = cos(t9);
t22 = t6 * rSges(6,1) - t4 * rSges(6,2) + pkin(4) * t13 + pkin(3);
t21 = rSges(5,1) * t13 - rSges(5,2) * sin(pkin(9)) + pkin(3);
t17 = cos(qJ(2));
t8 = t17 * pkin(2);
t5 = sin(t10);
t1 = [(-m(2) - m(3) + t24) * g(3), -m(3) * (g(3) * (rSges(3,1) * t17 - t16 * rSges(3,2)) + t33 * (-rSges(3,1) * t16 - rSges(3,2) * t17)) - m(4) * (g(3) * (rSges(4,1) * t7 - rSges(4,2) * t5 + t8) + t33 * (-rSges(4,1) * t5 - rSges(4,2) * t7 - t31)) - m(5) * (g(3) * (t21 * t7 + t25 * t5 + t8) + t33 * (-t21 * t5 + t25 * t7 - t31)) - m(6) * (g(3) * (t22 * t7 + t26 * t5 + t8) + t33 * (-t22 * t5 + t26 * t7 - t31)), t24 * (g(1) * t12 - g(2) * t14), t32 * (-g(3) * t7 + t33 * t5), -m(6) * (g(1) * ((t12 * t6 - t4 * t27) * rSges(6,1) + (-t12 * t4 - t6 * t27) * rSges(6,2)) + g(2) * ((-t14 * t6 - t4 * t28) * rSges(6,1) + (t14 * t4 - t6 * t28) * rSges(6,2)) + g(3) * (-rSges(6,1) * t4 - rSges(6,2) * t6) * t5)];
taug = t1(:);
