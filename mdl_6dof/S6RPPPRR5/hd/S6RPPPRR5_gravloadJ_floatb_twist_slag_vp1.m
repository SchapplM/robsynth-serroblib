% Calculate Gravitation load on the joints for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:13
% EndTime: 2019-03-09 01:37:14
% DurationCPUTime: 0.38s
% Computational Cost: add. (179->77), mult. (325->111), div. (0->0), fcn. (345->8), ass. (0->33)
t22 = cos(qJ(1));
t12 = t22 * qJ(2);
t19 = sin(qJ(1));
t38 = -qJ(3) - pkin(1);
t49 = t22 * pkin(3) + t19 * t38 + t12;
t17 = sin(qJ(6));
t20 = cos(qJ(6));
t47 = m(7) * (rSges(7,1) * t20 - rSges(7,2) * t17 + pkin(5)) + m(6) * rSges(6,1);
t41 = rSges(7,3) + pkin(8);
t46 = -m(6) * rSges(6,2) + m(7) * t41;
t21 = cos(qJ(5));
t43 = t21 * pkin(5);
t42 = rSges(6,3) + pkin(7);
t40 = t17 * t21;
t39 = t20 * t21;
t37 = t22 * pkin(1) + t19 * qJ(2);
t35 = sin(pkin(9));
t34 = -m(5) - m(6) - m(7);
t33 = t22 * qJ(3) + t37;
t32 = -m(4) + t34;
t31 = t19 * pkin(3) + t33;
t16 = cos(pkin(9));
t6 = t19 * t16 + t22 * t35;
t30 = t6 * pkin(4) + t31;
t18 = sin(qJ(5));
t29 = t21 * rSges(6,1) - t18 * rSges(6,2);
t5 = t22 * t16 - t19 * t35;
t28 = t6 * t17 + t39 * t5;
t27 = -t6 * t20 + t40 * t5;
t24 = t5 * pkin(4) + t49;
t2 = -t5 * t17 + t39 * t6;
t1 = -t5 * t20 - t40 * t6;
t3 = [-m(2) * (g(1) * (-t19 * rSges(2,1) - t22 * rSges(2,2)) + g(2) * (t22 * rSges(2,1) - t19 * rSges(2,2))) - m(3) * (g(1) * (t22 * rSges(3,3) + t12 + (rSges(3,2) - pkin(1)) * t19) + g(2) * (-t22 * rSges(3,2) + t19 * rSges(3,3) + t37)) - m(4) * (g(1) * (t22 * rSges(4,1) + t12) + g(2) * (t22 * rSges(4,3) + t33) + (g(1) * (-rSges(4,3) + t38) + g(2) * rSges(4,1)) * t19) - m(5) * (g(1) * (t5 * rSges(5,1) - t6 * rSges(5,2) + t49) + g(2) * (t6 * rSges(5,1) + t5 * rSges(5,2) + t31)) - m(6) * (g(1) * t24 + g(2) * t30 + (g(1) * t42 + g(2) * t29) * t6 + (g(1) * t29 - g(2) * t42) * t5) - m(7) * (g(1) * (rSges(7,1) * t28 - rSges(7,2) * t27 + t6 * pkin(7) + t43 * t5 + t24) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t5 * pkin(7) + t43 * t6 + t30) + (g(1) * t5 + g(2) * t6) * t18 * t41) (-m(3) + t32) * (g(1) * t19 - g(2) * t22) t32 * (g(1) * t22 + g(2) * t19) t34 * g(3) (-t18 * t46 - t21 * t47) * g(3) + (g(1) * t6 - g(2) * t5) * (t47 * t18 - t46 * t21) -m(7) * (g(1) * (t1 * rSges(7,1) - t2 * rSges(7,2)) + g(2) * (rSges(7,1) * t27 + rSges(7,2) * t28) + g(3) * (-t17 * rSges(7,1) - t20 * rSges(7,2)) * t18)];
taug  = t3(:);
