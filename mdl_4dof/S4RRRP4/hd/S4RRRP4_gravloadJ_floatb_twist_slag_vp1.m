% Calculate Gravitation load on the joints for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:10
% EndTime: 2019-12-31 17:15:11
% DurationCPUTime: 0.23s
% Computational Cost: add. (127->46), mult. (155->61), div. (0->0), fcn. (122->6), ass. (0->25)
t38 = rSges(5,1) + pkin(3);
t11 = qJ(2) + qJ(3);
t7 = sin(t11);
t8 = cos(t11);
t26 = t8 * rSges(4,1) - rSges(4,2) * t7;
t37 = -rSges(5,2) * t8 - t38 * t7;
t13 = sin(qJ(1));
t15 = cos(qJ(1));
t36 = g(1) * t15 + g(2) * t13;
t25 = -rSges(5,2) * t7 + t38 * t8;
t16 = -pkin(6) - pkin(5);
t14 = cos(qJ(2));
t9 = t14 * pkin(2);
t6 = t9 + pkin(1);
t12 = sin(qJ(2));
t32 = pkin(2) * t12;
t29 = rSges(3,3) + pkin(5);
t28 = rSges(4,3) - t16;
t27 = rSges(5,3) + qJ(4) - t16;
t24 = -rSges(4,1) * t7 - rSges(4,2) * t8;
t22 = rSges(3,1) * t14 - rSges(3,2) * t12;
t21 = t6 + t26;
t20 = t6 + t25;
t19 = pkin(1) + t22;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t13 - rSges(2,2) * t15) + g(2) * (rSges(2,1) * t15 - rSges(2,2) * t13)) - m(3) * ((g(1) * t29 + g(2) * t19) * t15 + (-g(1) * t19 + g(2) * t29) * t13) - m(4) * ((g(1) * t28 + g(2) * t21) * t15 + (-g(1) * t21 + g(2) * t28) * t13) - m(5) * ((g(1) * t27 + g(2) * t20) * t15 + (-g(1) * t20 + g(2) * t27) * t13), (-m(3) * t22 - m(4) * (t26 + t9) - m(5) * (t25 + t9)) * g(3) + t36 * (-m(3) * (-rSges(3,1) * t12 - rSges(3,2) * t14) - m(4) * (t24 - t32) - m(5) * (-t32 + t37)), (-m(4) * t26 - m(5) * t25) * g(3) + t36 * (-m(4) * t24 - m(5) * t37), -m(5) * (g(1) * t13 - g(2) * t15)];
taug = t1(:);
