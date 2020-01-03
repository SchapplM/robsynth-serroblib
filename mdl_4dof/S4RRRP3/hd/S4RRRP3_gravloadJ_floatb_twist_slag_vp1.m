% Calculate Gravitation load on the joints for
% S4RRRP3
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:01
% EndTime: 2019-12-31 17:14:02
% DurationCPUTime: 0.21s
% Computational Cost: add. (154->49), mult. (155->66), div. (0->0), fcn. (125->6), ass. (0->26)
t37 = rSges(5,1) + pkin(3);
t19 = sin(qJ(3));
t21 = cos(qJ(3));
t43 = rSges(4,1) * t21 - rSges(4,2) * t19;
t30 = rSges(5,3) + qJ(4);
t18 = qJ(1) + qJ(2);
t15 = sin(t18);
t16 = cos(t18);
t42 = g(1) * t16 + g(2) * t15;
t41 = t30 * t19 + t37 * t21;
t20 = sin(qJ(1));
t40 = pkin(1) * t20;
t34 = t16 * t19;
t33 = t16 * t21;
t13 = t16 * pkin(6);
t32 = t16 * rSges(5,2) + t13;
t31 = t16 * pkin(2) + t15 * pkin(6);
t29 = t16 * rSges(3,1) - rSges(3,2) * t15;
t28 = -rSges(3,1) * t15 - rSges(3,2) * t16;
t27 = t15 * rSges(5,2) + t30 * t34 + t37 * t33 + t31;
t26 = rSges(4,1) * t33 - rSges(4,2) * t34 + t15 * rSges(4,3) + t31;
t25 = t16 * rSges(4,3) + t13 + (-pkin(2) - t43) * t15;
t24 = g(1) * (-pkin(2) - t41) * t15;
t22 = cos(qJ(1));
t17 = t22 * pkin(1);
t1 = [-m(2) * (g(1) * (-t20 * rSges(2,1) - rSges(2,2) * t22) + g(2) * (rSges(2,1) * t22 - t20 * rSges(2,2))) - m(3) * (g(1) * (t28 - t40) + g(2) * (t17 + t29)) - m(4) * (g(1) * (t25 - t40) + g(2) * (t17 + t26)) - m(5) * (g(1) * (t32 - t40) + g(2) * (t17 + t27) + t24), -m(3) * (g(1) * t28 + g(2) * t29) - m(4) * (g(1) * t25 + g(2) * t26) - m(5) * (g(1) * t32 + g(2) * t27 + t24), (-m(4) * t43 - m(5) * t41) * g(3) + t42 * (-m(4) * (-rSges(4,1) * t19 - rSges(4,2) * t21) - m(5) * (-t37 * t19 + t30 * t21)), -m(5) * (-g(3) * t21 + t42 * t19)];
taug = t1(:);
