% Calculate Gravitation load on the joints for
% S4RRRP2
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:53
% EndTime: 2019-12-31 17:12:54
% DurationCPUTime: 0.23s
% Computational Cost: add. (140->46), mult. (135->64), div. (0->0), fcn. (105->6), ass. (0->25)
t42 = rSges(5,1) + pkin(3);
t23 = cos(qJ(3));
t41 = t42 * t23;
t35 = rSges(4,1) * t23;
t40 = pkin(2) + t35;
t39 = pkin(6) + rSges(4,3);
t38 = pkin(2) + t41;
t37 = qJ(4) + pkin(6) + rSges(5,3);
t22 = sin(qJ(1));
t36 = pkin(1) * t22;
t19 = qJ(1) + qJ(2);
t15 = sin(t19);
t21 = sin(qJ(3));
t34 = t15 * t21;
t16 = cos(t19);
t33 = t16 * t21;
t31 = t16 * rSges(3,1) - rSges(3,2) * t15;
t30 = -rSges(3,1) * t15 - rSges(3,2) * t16;
t29 = -rSges(4,2) * t33 + t39 * t15 + t40 * t16;
t28 = rSges(4,2) * t34 - t40 * t15 + t39 * t16;
t27 = -rSges(5,2) * t33 + t37 * t15 + t38 * t16;
t26 = rSges(5,2) * t34 - t38 * t15 + t37 * t16;
t24 = cos(qJ(1));
t18 = t24 * pkin(1);
t1 = [-m(2) * (g(1) * (-t22 * rSges(2,1) - rSges(2,2) * t24) + g(2) * (rSges(2,1) * t24 - t22 * rSges(2,2))) - m(3) * (g(1) * (t30 - t36) + g(2) * (t18 + t31)) - m(4) * (g(1) * (t28 - t36) + g(2) * (t18 + t29)) - m(5) * (g(1) * (t26 - t36) + g(2) * (t18 + t27)), -m(3) * (g(1) * t30 + g(2) * t31) - m(4) * (g(1) * t28 + g(2) * t29) - m(5) * (g(1) * t26 + g(2) * t27), (-m(4) * (-rSges(4,2) * t21 + t35) - m(5) * (-rSges(5,2) * t21 + t41)) * g(3) + (g(1) * t16 + g(2) * t15) * (-m(4) * (-rSges(4,1) * t21 - rSges(4,2) * t23) - m(5) * (-rSges(5,2) * t23 - t42 * t21)), -m(5) * (g(1) * t15 - g(2) * t16)];
taug = t1(:);
