% Calculate Gravitation load on the joints for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:08
% EndTime: 2019-12-31 17:22:08
% DurationCPUTime: 0.18s
% Computational Cost: add. (179->42), mult. (121->56), div. (0->0), fcn. (90->8), ass. (0->28)
t39 = rSges(5,3) + pkin(7);
t18 = sin(qJ(4));
t20 = cos(qJ(4));
t38 = rSges(5,1) * t20 - rSges(5,2) * t18;
t37 = -pkin(3) - t38;
t19 = sin(qJ(1));
t36 = pkin(1) * t19;
t17 = qJ(1) + qJ(2);
t13 = sin(t17);
t35 = pkin(2) * t13;
t14 = cos(t17);
t32 = t14 * rSges(3,1) - rSges(3,2) * t13;
t15 = qJ(3) + t17;
t11 = sin(t15);
t12 = cos(t15);
t31 = t12 * rSges(4,1) - t11 * rSges(4,2);
t10 = pkin(2) * t14;
t30 = t10 + t31;
t29 = -rSges(3,1) * t13 - rSges(3,2) * t14;
t28 = -rSges(4,1) * t11 - rSges(4,2) * t12;
t26 = t39 * t11 - t37 * t12;
t25 = t28 - t35;
t24 = t10 + t26;
t23 = t37 * t11 + t39 * t12;
t22 = t23 - t35;
t21 = cos(qJ(1));
t16 = t21 * pkin(1);
t1 = [-m(2) * (g(1) * (-t19 * rSges(2,1) - rSges(2,2) * t21) + g(2) * (rSges(2,1) * t21 - t19 * rSges(2,2))) - m(3) * (g(1) * (t29 - t36) + g(2) * (t16 + t32)) - m(4) * (g(1) * (t25 - t36) + g(2) * (t16 + t30)) - m(5) * (g(1) * (t22 - t36) + g(2) * (t16 + t24)), -m(3) * (g(1) * t29 + g(2) * t32) - m(4) * (g(1) * t25 + g(2) * t30) - m(5) * (g(1) * t22 + g(2) * t24), -m(4) * (g(1) * t28 + g(2) * t31) - m(5) * (g(1) * t23 + g(2) * t26), -m(5) * (g(3) * t38 + (g(1) * t12 + g(2) * t11) * (-rSges(5,1) * t18 - rSges(5,2) * t20))];
taug = t1(:);
