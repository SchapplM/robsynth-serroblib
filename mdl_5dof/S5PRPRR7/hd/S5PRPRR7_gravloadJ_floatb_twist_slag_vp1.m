% Calculate Gravitation load on the joints for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:44
% EndTime: 2019-12-05 15:59:46
% DurationCPUTime: 0.32s
% Computational Cost: add. (131->58), mult. (233->88), div. (0->0), fcn. (209->8), ass. (0->27)
t19 = cos(qJ(2));
t41 = g(3) * t19;
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t44 = g(1) * t15 + g(2) * t14;
t40 = rSges(5,3) + pkin(6);
t17 = sin(qJ(2));
t39 = t14 * t17;
t38 = t15 * t17;
t16 = sin(qJ(4));
t37 = t16 * t17;
t18 = cos(qJ(4));
t36 = t17 * t18;
t35 = rSges(6,3) + pkin(7) + pkin(6);
t34 = t19 * pkin(2) + t17 * qJ(3);
t32 = -m(4) - m(5) - m(6);
t29 = t44 * qJ(3) * t19;
t27 = t16 * rSges(5,1) + t18 * rSges(5,2);
t26 = -t14 * t16 + t15 * t36;
t25 = t14 * t36 + t15 * t16;
t13 = qJ(4) + qJ(5);
t10 = cos(t13);
t9 = sin(t13);
t23 = rSges(6,1) * t9 + rSges(6,2) * t10 + pkin(4) * t16;
t22 = g(3) * t34 + t29;
t21 = g(1) * ((t10 * t38 - t14 * t9) * rSges(6,1) + (-t14 * t10 - t9 * t38) * rSges(6,2)) + g(2) * ((t10 * t39 + t15 * t9) * rSges(6,1) + (t15 * t10 - t9 * t39) * rSges(6,2)) + (-rSges(6,1) * t10 + t9 * rSges(6,2)) * t41;
t1 = [(-m(2) - m(3) + t32) * g(3), -m(3) * (g(3) * (t19 * rSges(3,1) - t17 * rSges(3,2)) + t44 * (-rSges(3,1) * t17 - rSges(3,2) * t19)) - m(4) * (g(3) * (-t19 * rSges(4,2) + t17 * rSges(4,3) + t34) + t29 + t44 * (rSges(4,3) * t19 + (rSges(4,2) - pkin(2)) * t17)) - m(5) * ((g(3) * t40 + t44 * t27) * t19 + (g(3) * t27 + t44 * (-pkin(2) - t40)) * t17 + t22) - m(6) * ((g(3) * t35 + t44 * t23) * t19 + (g(3) * t23 + t44 * (-pkin(2) - t35)) * t17 + t22), t32 * (t44 * t17 - t41), -m(5) * (g(1) * (t26 * rSges(5,1) + (-t14 * t18 - t15 * t37) * rSges(5,2)) + g(2) * (t25 * rSges(5,1) + (-t14 * t37 + t15 * t18) * rSges(5,2)) + (-rSges(5,1) * t18 + rSges(5,2) * t16) * t41) - m(6) * ((g(1) * t26 + g(2) * t25 - t18 * t41) * pkin(4) + t21), -m(6) * t21];
taug = t1(:);
