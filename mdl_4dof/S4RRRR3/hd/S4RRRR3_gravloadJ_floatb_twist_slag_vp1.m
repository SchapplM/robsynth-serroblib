% Calculate Gravitation load on the joints for
% S4RRRR3
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
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:17
% EndTime: 2019-12-31 17:24:17
% DurationCPUTime: 0.24s
% Computational Cost: add. (159->48), mult. (166->63), div. (0->0), fcn. (130->8), ass. (0->29)
t14 = qJ(2) + qJ(3);
t10 = cos(t14);
t11 = qJ(4) + t14;
t6 = sin(t11);
t7 = cos(t11);
t29 = t7 * rSges(5,1) - t6 * rSges(5,2);
t28 = pkin(3) * t10 + t29;
t9 = sin(t14);
t30 = t10 * rSges(4,1) - rSges(4,2) * t9;
t27 = -rSges(5,1) * t6 - rSges(5,2) * t7;
t41 = -pkin(3) * t9 + t27;
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t40 = g(1) * t18 + g(2) * t16;
t19 = -pkin(6) - pkin(5);
t15 = sin(qJ(2));
t37 = pkin(2) * t15;
t33 = rSges(3,3) + pkin(5);
t17 = cos(qJ(2));
t12 = t17 * pkin(2);
t8 = t12 + pkin(1);
t32 = rSges(4,3) - t19;
t31 = rSges(5,3) + pkin(7) - t19;
t26 = -rSges(4,1) * t9 - rSges(4,2) * t10;
t25 = rSges(3,1) * t17 - rSges(3,2) * t15;
t24 = t8 + t28;
t23 = t8 + t30;
t22 = pkin(1) + t25;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t16 - rSges(2,2) * t18) + g(2) * (rSges(2,1) * t18 - rSges(2,2) * t16)) - m(3) * ((g(1) * t33 + g(2) * t22) * t18 + (-g(1) * t22 + g(2) * t33) * t16) - m(4) * ((g(1) * t32 + g(2) * t23) * t18 + (-g(1) * t23 + g(2) * t32) * t16) - m(5) * ((g(1) * t31 + g(2) * t24) * t18 + (-g(1) * t24 + g(2) * t31) * t16), (-m(3) * t25 - m(4) * (t12 + t30) - m(5) * (t12 + t28)) * g(3) + t40 * (-m(3) * (-rSges(3,1) * t15 - rSges(3,2) * t17) - m(4) * (t26 - t37) - m(5) * (-t37 + t41)), (-m(4) * t30 - m(5) * t28) * g(3) + t40 * (-m(4) * t26 - m(5) * t41), -m(5) * (g(3) * t29 + t27 * t40)];
taug = t1(:);
