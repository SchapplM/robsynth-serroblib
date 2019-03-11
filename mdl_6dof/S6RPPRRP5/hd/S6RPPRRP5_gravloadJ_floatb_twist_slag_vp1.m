% Calculate Gravitation load on the joints for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:41
% EndTime: 2019-03-09 02:07:43
% DurationCPUTime: 0.42s
% Computational Cost: add. (169->91), mult. (316->124), div. (0->0), fcn. (284->6), ass. (0->33)
t46 = rSges(7,1) + pkin(5);
t45 = rSges(6,3) + pkin(8);
t44 = rSges(7,3) + qJ(6) + pkin(8);
t17 = sin(qJ(1));
t20 = cos(qJ(1));
t28 = g(1) * t20 + g(2) * t17;
t15 = sin(qJ(5));
t18 = cos(qJ(5));
t9 = t18 * pkin(5) + pkin(4);
t43 = m(6) * (rSges(6,1) * t18 - rSges(6,2) * t15 + pkin(4)) + m(7) * (rSges(7,1) * t18 - rSges(7,2) * t15 + t9) + m(5) * rSges(5,1);
t41 = pkin(5) * t15;
t38 = -rSges(5,3) - pkin(7);
t37 = t17 * t15;
t36 = t17 * t18;
t35 = t20 * t15;
t34 = t20 * t18;
t33 = -pkin(1) - qJ(3);
t32 = t20 * pkin(1) + t17 * qJ(2);
t31 = t20 * qJ(3) + t32;
t30 = -m(4) - m(5) - m(6) - m(7);
t29 = -pkin(7) - t41;
t16 = sin(qJ(4));
t19 = cos(qJ(4));
t27 = t16 * rSges(5,1) + t19 * rSges(5,2);
t3 = -t16 * t35 - t36;
t1 = t16 * t37 - t34;
t24 = t16 * pkin(4) - t45 * t19;
t23 = t16 * t9 - t44 * t19;
t22 = m(5) * rSges(5,2) - m(6) * t45 - m(7) * t44;
t12 = t20 * qJ(2);
t4 = t16 * t34 - t37;
t2 = -t16 * t36 - t35;
t5 = [-m(2) * (g(1) * (-t17 * rSges(2,1) - t20 * rSges(2,2)) + g(2) * (t20 * rSges(2,1) - t17 * rSges(2,2))) - m(3) * (g(1) * (t20 * rSges(3,3) + t12 + (rSges(3,2) - pkin(1)) * t17) + g(2) * (-t20 * rSges(3,2) + t17 * rSges(3,3) + t32)) - m(4) * (g(1) * (t20 * rSges(4,2) + t12) + g(2) * (t20 * rSges(4,3) + t31) + (g(1) * (-rSges(4,3) + t33) + g(2) * rSges(4,2)) * t17) - m(5) * (g(1) * t12 + g(2) * t31 + (g(1) * t38 + g(2) * t27) * t20 + (g(1) * (-t27 + t33) + g(2) * t38) * t17) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t12) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t31) + (-g(1) * pkin(7) + g(2) * t24) * t20 + (g(1) * (-t24 + t33) - g(2) * pkin(7)) * t17) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t12) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t31) + (g(1) * t29 + g(2) * t23) * t20 + (g(1) * (-t23 + t33) + g(2) * t29) * t17) (-m(3) + t30) * (g(1) * t17 - g(2) * t20) t30 * t28 (t43 * t16 + t22 * t19) * g(3) + t28 * (t22 * t16 - t43 * t19) -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) + t2 * rSges(6,2))) - m(7) * (g(1) * (-t4 * rSges(7,2) + t46 * t3) + g(2) * (t2 * rSges(7,2) - t46 * t1)) + (-m(6) * (-rSges(6,1) * t15 - rSges(6,2) * t18) - m(7) * (-rSges(7,1) * t15 - rSges(7,2) * t18 - t41)) * g(3) * t19, -m(7) * (g(3) * t16 - t28 * t19)];
taug  = t5(:);
