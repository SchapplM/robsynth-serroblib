% Calculate Gravitation load on the joints for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:33
% EndTime: 2019-12-05 17:39:35
% DurationCPUTime: 0.22s
% Computational Cost: add. (138->60), mult. (153->76), div. (0->0), fcn. (119->8), ass. (0->31)
t15 = pkin(8) + qJ(4);
t10 = qJ(5) + t15;
t6 = sin(t10);
t7 = cos(t10);
t25 = rSges(6,1) * t6 + rSges(6,2) * t7;
t8 = sin(t15);
t41 = -pkin(4) * t8 - t25;
t19 = sin(qJ(1));
t20 = cos(qJ(1));
t40 = g(1) * t19 - g(2) * t20;
t9 = cos(t15);
t38 = pkin(4) * t9;
t37 = rSges(6,1) * t7;
t36 = rSges(6,2) * t6;
t16 = sin(pkin(8));
t35 = pkin(3) * t16;
t18 = -pkin(6) - qJ(3);
t32 = rSges(5,3) - t18;
t31 = rSges(6,3) + pkin(7) - t18;
t30 = t20 * pkin(1) + t19 * qJ(2);
t29 = rSges(4,3) + qJ(3);
t28 = -m(4) - m(5) - m(6);
t26 = -rSges(5,1) * t8 - rSges(5,2) * t9;
t24 = rSges(4,1) * t16 + rSges(4,2) * cos(pkin(8));
t23 = t35 - t41;
t12 = t20 * qJ(2);
t22 = g(1) * t12 + g(2) * t30;
t21 = -t26 + t35;
t3 = t20 * t36;
t2 = t19 * t37;
t1 = [-m(2) * (g(1) * (-t19 * rSges(2,1) - rSges(2,2) * t20) + g(2) * (rSges(2,1) * t20 - t19 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t20 + t12 + (rSges(3,2) - pkin(1)) * t19) + g(2) * (-rSges(3,2) * t20 + t19 * rSges(3,3) + t30)) - m(4) * ((g(1) * t24 + g(2) * t29) * t20 + (g(1) * (-pkin(1) - t29) + g(2) * t24) * t19 + t22) - m(5) * ((g(1) * t21 + g(2) * t32) * t20 + (g(1) * (-pkin(1) - t32) + g(2) * t21) * t19 + t22) - m(6) * ((g(1) * t23 + g(2) * t31) * t20 + (g(1) * (-pkin(1) - t31) + g(2) * t23) * t19 + t22), (-m(3) + t28) * t40, t28 * (g(1) * t20 + g(2) * t19), -m(5) * (g(3) * t26 + t40 * (rSges(5,1) * t9 - rSges(5,2) * t8)) - m(6) * (g(1) * (t2 + (-t36 + t38) * t19) + g(2) * (t3 + (-t37 - t38) * t20) + g(3) * t41), -m(6) * (g(1) * (-t19 * t36 + t2) + g(2) * (-t20 * t37 + t3) - g(3) * t25)];
taug = t1(:);
