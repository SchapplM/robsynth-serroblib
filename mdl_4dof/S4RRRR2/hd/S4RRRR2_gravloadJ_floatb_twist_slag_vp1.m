% Calculate Gravitation load on the joints for
% S4RRRR2
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
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:07
% EndTime: 2019-12-31 17:23:08
% DurationCPUTime: 0.19s
% Computational Cost: add. (166->44), mult. (146->60), div. (0->0), fcn. (113->8), ass. (0->28)
t26 = cos(qJ(3));
t22 = qJ(3) + qJ(4);
t16 = sin(t22);
t18 = cos(t22);
t37 = t18 * rSges(5,1) - rSges(5,2) * t16;
t50 = t26 * pkin(3) + t37;
t49 = rSges(4,3) + pkin(6);
t48 = pkin(7) + pkin(6) + rSges(5,3);
t24 = sin(qJ(3));
t47 = rSges(4,1) * t26 - rSges(4,2) * t24;
t46 = -pkin(2) - t50;
t45 = -pkin(2) - t47;
t23 = qJ(1) + qJ(2);
t17 = sin(t23);
t19 = cos(t23);
t44 = g(1) * t19 + g(2) * t17;
t25 = sin(qJ(1));
t41 = t25 * pkin(1);
t36 = t19 * rSges(3,1) - rSges(3,2) * t17;
t35 = -rSges(3,1) * t17 - rSges(3,2) * t19;
t34 = -rSges(5,1) * t16 - rSges(5,2) * t18;
t33 = t49 * t17 - t45 * t19;
t32 = t45 * t17 + t49 * t19;
t31 = t48 * t17 - t46 * t19;
t30 = t46 * t17 + t48 * t19;
t27 = cos(qJ(1));
t21 = t27 * pkin(1);
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t25 - rSges(2,2) * t27) + g(2) * (rSges(2,1) * t27 - rSges(2,2) * t25)) - m(3) * (g(1) * (t35 - t41) + g(2) * (t21 + t36)) - m(4) * (g(1) * (t32 - t41) + g(2) * (t21 + t33)) - m(5) * (g(1) * (t30 - t41) + g(2) * (t21 + t31)), -m(3) * (g(1) * t35 + g(2) * t36) - m(4) * (g(1) * t32 + g(2) * t33) - m(5) * (g(1) * t30 + g(2) * t31), (-m(4) * t47 - m(5) * t50) * g(3) + t44 * (-m(4) * (-rSges(4,1) * t24 - rSges(4,2) * t26) - m(5) * (-pkin(3) * t24 + t34)), -m(5) * (g(3) * t37 + t44 * t34)];
taug = t1(:);
