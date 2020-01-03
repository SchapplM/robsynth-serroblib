% Calculate Gravitation load on the joints for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:07
% EndTime: 2019-12-31 16:48:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (124->37), mult. (89->48), div. (0->0), fcn. (64->8), ass. (0->22)
t33 = rSges(5,3) + pkin(6);
t17 = sin(qJ(4));
t19 = cos(qJ(4));
t32 = rSges(5,1) * t19 - rSges(5,2) * t17;
t31 = -pkin(3) - t32;
t18 = sin(qJ(1));
t30 = pkin(1) * t18;
t16 = qJ(1) + pkin(7);
t13 = cos(t16);
t20 = cos(qJ(1));
t15 = t20 * pkin(1);
t29 = pkin(2) * t13 + t15;
t14 = qJ(3) + t16;
t10 = sin(t14);
t11 = cos(t14);
t26 = t11 * rSges(4,1) - rSges(4,2) * t10;
t12 = sin(t16);
t25 = -pkin(2) * t12 - t30;
t24 = -rSges(4,1) * t10 - rSges(4,2) * t11;
t22 = t33 * t10 - t31 * t11;
t21 = t31 * t10 + t33 * t11;
t1 = [-m(2) * (g(1) * (-t18 * rSges(2,1) - rSges(2,2) * t20) + g(2) * (rSges(2,1) * t20 - t18 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t12 - rSges(3,2) * t13 - t30) + g(2) * (rSges(3,1) * t13 - rSges(3,2) * t12 + t15)) - m(4) * (g(1) * (t24 + t25) + g(2) * (t26 + t29)) - m(5) * (g(1) * (t21 + t25) + g(2) * (t22 + t29)), (-m(3) - m(4) - m(5)) * g(3), -m(4) * (g(1) * t24 + g(2) * t26) - m(5) * (g(1) * t21 + g(2) * t22), -m(5) * (g(3) * t32 + (g(1) * t11 + g(2) * t10) * (-rSges(5,1) * t17 - rSges(5,2) * t19))];
taug = t1(:);
