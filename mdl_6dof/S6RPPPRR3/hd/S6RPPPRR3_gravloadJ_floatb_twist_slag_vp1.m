% Calculate Gravitation load on the joints for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:28
% EndTime: 2019-03-09 01:33:29
% DurationCPUTime: 0.37s
% Computational Cost: add. (243->77), mult. (363->110), div. (0->0), fcn. (399->10), ass. (0->37)
t41 = sin(pkin(9));
t42 = cos(pkin(9));
t47 = sin(qJ(1));
t48 = cos(qJ(1));
t8 = -t47 * t41 - t48 * t42;
t9 = t48 * t41 - t47 * t42;
t55 = g(1) * t9 - g(2) * t8;
t25 = sin(qJ(6));
t26 = cos(qJ(6));
t54 = m(7) * (rSges(7,1) * t26 - rSges(7,2) * t25 + pkin(5)) + m(6) * rSges(6,1);
t21 = pkin(10) + qJ(5);
t16 = cos(t21);
t50 = t16 * pkin(5);
t49 = rSges(7,3) + pkin(8);
t46 = t16 * t25;
t45 = t16 * t26;
t44 = t48 * pkin(1) + t47 * qJ(2);
t43 = rSges(5,3) + qJ(4);
t40 = -m(5) - m(6) - m(7);
t39 = t48 * pkin(2) + t44;
t38 = m(4) - t40;
t23 = cos(pkin(10));
t14 = t23 * pkin(4) + pkin(3);
t24 = -pkin(7) - qJ(4);
t37 = -t8 * t14 - t9 * t24 + t39;
t36 = -t47 * pkin(1) + t48 * qJ(2);
t15 = sin(t21);
t35 = t16 * rSges(6,1) - t15 * rSges(6,2);
t34 = t8 * t25 + t9 * t45;
t33 = -t8 * t26 + t9 * t46;
t32 = rSges(5,1) * t23 - rSges(5,2) * sin(pkin(10)) + pkin(3);
t30 = -m(6) * rSges(6,2) + m(7) * t49;
t29 = -t47 * pkin(2) + t36;
t28 = t9 * t14 - t8 * t24 + t29;
t2 = t9 * t25 - t8 * t45;
t1 = t9 * t26 + t8 * t46;
t3 = [-m(2) * (g(1) * (-t47 * rSges(2,1) - t48 * rSges(2,2)) + g(2) * (t48 * rSges(2,1) - t47 * rSges(2,2))) - m(3) * (g(1) * (-t47 * rSges(3,1) + t48 * rSges(3,3) + t36) + g(2) * (t48 * rSges(3,1) + t47 * rSges(3,3) + t44)) - m(4) * (g(1) * (t9 * rSges(4,1) - t8 * rSges(4,2) + t29) + g(2) * (-t8 * rSges(4,1) - t9 * rSges(4,2) + t39)) - m(5) * (g(1) * t29 + g(2) * t39 + (g(1) * t32 + g(2) * t43) * t9 + (g(1) * t43 - g(2) * t32) * t8) - m(6) * (g(1) * (t8 * rSges(6,3) + t35 * t9 + t28) + g(2) * (t9 * rSges(6,3) - t35 * t8 + t37)) - m(7) * (g(1) * (t34 * rSges(7,1) - t33 * rSges(7,2) + t9 * t50 + t28) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t8 * t50 + t37) + t55 * t15 * t49) (-m(3) - t38) * (g(1) * t47 - g(2) * t48) t38 * g(3), t40 * t55 (t30 * t15 + t16 * t54) * g(3) + (g(1) * t8 + g(2) * t9) * (-t15 * t54 + t30 * t16) -m(7) * (g(1) * (t1 * rSges(7,1) - t2 * rSges(7,2)) + g(2) * (t33 * rSges(7,1) + t34 * rSges(7,2)) + g(3) * (t25 * rSges(7,1) + t26 * rSges(7,2)) * t15)];
taug  = t3(:);
