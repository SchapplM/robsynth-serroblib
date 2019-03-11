% Calculate Gravitation load on the joints for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:48:31
% EndTime: 2019-03-09 01:48:33
% DurationCPUTime: 0.44s
% Computational Cost: add. (181->82), mult. (286->111), div. (0->0), fcn. (253->8), ass. (0->36)
t19 = sin(qJ(4));
t21 = cos(qJ(4));
t16 = sin(pkin(9));
t17 = cos(pkin(9));
t28 = rSges(6,1) * t17 - rSges(6,2) * t16 + pkin(4);
t35 = rSges(6,3) + qJ(5);
t49 = t28 * t19 - t21 * t35;
t46 = rSges(7,3) + pkin(8) + qJ(5);
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t31 = g(1) * t22 + g(2) * t20;
t15 = pkin(9) + qJ(6);
t10 = cos(t15);
t8 = t17 * pkin(5) + pkin(4);
t9 = sin(t15);
t45 = m(6) * t28 + m(7) * (rSges(7,1) * t10 - rSges(7,2) * t9 + t8) + m(5) * rSges(5,1);
t44 = -m(6) - m(7);
t40 = -rSges(5,3) - pkin(7);
t39 = t20 * t19;
t38 = t22 * t19;
t37 = -pkin(1) - qJ(3);
t36 = t22 * pkin(1) + t20 * qJ(2);
t34 = t22 * qJ(3) + t36;
t33 = -m(4) - m(5) + t44;
t32 = -pkin(5) * t16 - pkin(7);
t30 = t19 * rSges(5,1) + t21 * rSges(5,2);
t27 = -t16 * rSges(6,1) - t17 * rSges(6,2) - pkin(7);
t13 = t22 * qJ(2);
t26 = g(1) * t13 + g(2) * t34;
t25 = t19 * t8 - t46 * t21;
t24 = m(5) * rSges(5,2) - m(6) * t35 - m(7) * t46;
t4 = t10 * t38 - t20 * t9;
t3 = -t20 * t10 - t9 * t38;
t2 = -t10 * t39 - t22 * t9;
t1 = -t22 * t10 + t9 * t39;
t5 = [-m(2) * (g(1) * (-t20 * rSges(2,1) - t22 * rSges(2,2)) + g(2) * (t22 * rSges(2,1) - t20 * rSges(2,2))) - m(3) * (g(1) * (t22 * rSges(3,3) + t13 + (rSges(3,2) - pkin(1)) * t20) + g(2) * (-t22 * rSges(3,2) + t20 * rSges(3,3) + t36)) - m(4) * (g(1) * (t22 * rSges(4,2) + t13) + g(2) * (t22 * rSges(4,3) + t34) + (g(1) * (-rSges(4,3) + t37) + g(2) * rSges(4,2)) * t20) - m(5) * ((g(1) * t40 + g(2) * t30) * t22 + (g(1) * (-t30 + t37) + g(2) * t40) * t20 + t26) - m(6) * ((g(1) * t27 + t49 * g(2)) * t22 + (g(2) * t27 + (t37 - t49) * g(1)) * t20 + t26) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t13) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t34) + (g(1) * t32 + g(2) * t25) * t22 + (g(1) * (-t25 + t37) + g(2) * t32) * t20) (-m(3) + t33) * (g(1) * t20 - g(2) * t22) t33 * t31 (t45 * t19 + t24 * t21) * g(3) + t31 * (t24 * t19 - t45 * t21) t44 * (g(3) * t19 - t31 * t21) -m(7) * (g(1) * (t3 * rSges(7,1) - t4 * rSges(7,2)) + g(2) * (-t1 * rSges(7,1) + t2 * rSges(7,2)) + g(3) * (-rSges(7,1) * t9 - rSges(7,2) * t10) * t21)];
taug  = t5(:);
