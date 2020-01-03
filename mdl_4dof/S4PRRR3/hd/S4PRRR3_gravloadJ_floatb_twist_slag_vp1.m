% Calculate Gravitation load on the joints for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:34
% EndTime: 2019-12-31 16:31:34
% DurationCPUTime: 0.12s
% Computational Cost: add. (115->30), mult. (77->39), div. (0->0), fcn. (54->6), ass. (0->18)
t28 = rSges(5,3) + pkin(6);
t16 = sin(qJ(4));
t17 = cos(qJ(4));
t27 = rSges(5,1) * t17 - t16 * rSges(5,2);
t26 = -pkin(3) - t27;
t15 = pkin(7) + qJ(2);
t12 = sin(t15);
t25 = pkin(2) * t12;
t14 = qJ(3) + t15;
t10 = sin(t14);
t11 = cos(t14);
t22 = t11 * rSges(4,1) - rSges(4,2) * t10;
t21 = -rSges(4,1) * t10 - rSges(4,2) * t11;
t19 = t28 * t10 - t26 * t11;
t18 = t26 * t10 + t28 * t11;
t13 = cos(t15);
t9 = pkin(2) * t13;
t1 = [(-m(2) - m(3) - m(4) - m(5)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t12 - rSges(3,2) * t13) + g(2) * (rSges(3,1) * t13 - rSges(3,2) * t12)) - m(4) * (g(1) * (t21 - t25) + g(2) * (t22 + t9)) - m(5) * (g(1) * (t18 - t25) + g(2) * (t19 + t9)), -m(4) * (g(1) * t21 + g(2) * t22) - m(5) * (g(1) * t18 + g(2) * t19), -m(5) * (g(3) * t27 + (g(1) * t11 + g(2) * t10) * (-rSges(5,1) * t16 - rSges(5,2) * t17))];
taug = t1(:);
