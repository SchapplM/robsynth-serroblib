% Calculate Gravitation load on the joints for
% S4RRRP5
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
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:38
% EndTime: 2019-12-31 17:16:39
% DurationCPUTime: 0.30s
% Computational Cost: add. (145->50), mult. (176->63), div. (0->0), fcn. (143->6), ass. (0->26)
t50 = rSges(5,3) + qJ(4);
t49 = rSges(5,1) + pkin(3);
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t43 = g(1) * t19 + g(2) * t17;
t15 = qJ(2) + qJ(3);
t12 = sin(t15);
t13 = cos(t15);
t44 = t13 * rSges(4,1) - rSges(4,2) * t12;
t42 = t50 * t12 + t49 * t13;
t41 = t43 * t12;
t18 = cos(qJ(2));
t14 = t18 * pkin(2);
t11 = t14 + pkin(1);
t40 = g(2) * t19 * t11;
t16 = sin(qJ(2));
t39 = pkin(2) * t16;
t35 = rSges(3,3) + pkin(5);
t20 = -pkin(6) - pkin(5);
t32 = rSges(5,2) - t20;
t31 = rSges(4,3) - t20;
t27 = rSges(3,1) * t18 - rSges(3,2) * t16;
t25 = -rSges(4,1) * t12 - rSges(4,2) * t13;
t24 = pkin(1) + t27;
t23 = t43 * t50 * t13;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t17 - rSges(2,2) * t19) + g(2) * (rSges(2,1) * t19 - rSges(2,2) * t17)) - m(3) * ((g(1) * t35 + g(2) * t24) * t19 + (-g(1) * t24 + g(2) * t35) * t17) - m(4) * (t40 + (g(1) * t31 + g(2) * t44) * t19 + (g(1) * (-t11 - t44) + g(2) * t31) * t17) - m(5) * (t40 + (g(1) * t32 + g(2) * t42) * t19 + (g(1) * (-t11 - t42) + g(2) * t32) * t17), -m(5) * t23 + (-m(3) * t27 - m(4) * (t14 + t44) - m(5) * (t14 + t42)) * g(3) + t43 * (-m(3) * (-rSges(3,1) * t16 - rSges(3,2) * t18) - m(4) * (t25 - t39) - m(5) * (-t12 * t49 - t39)), -m(4) * (g(3) * t44 + t43 * t25) - m(5) * (g(3) * t42 - t41 * t49 + t23), -m(5) * (-g(3) * t13 + t41)];
taug = t1(:);
