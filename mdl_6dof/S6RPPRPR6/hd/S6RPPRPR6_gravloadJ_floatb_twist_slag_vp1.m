% Calculate Gravitation load on the joints for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:43
% EndTime: 2019-03-09 01:50:44
% DurationCPUTime: 0.41s
% Computational Cost: add. (145->92), mult. (267->119), div. (0->0), fcn. (230->6), ass. (0->38)
t24 = cos(qJ(1));
t54 = g(1) * t24;
t53 = rSges(7,3) + pkin(8);
t19 = sin(qJ(6));
t22 = cos(qJ(6));
t52 = rSges(7,1) * t19 + rSges(7,2) * t22;
t21 = sin(qJ(1));
t29 = g(2) * t21 + t54;
t51 = -m(6) - m(7);
t50 = -pkin(5) - pkin(7);
t20 = sin(qJ(4));
t47 = g(3) * t20;
t46 = t20 * pkin(4);
t45 = -rSges(6,1) - pkin(7);
t44 = -rSges(5,3) - pkin(7);
t23 = cos(qJ(4));
t15 = t23 * qJ(5);
t17 = t24 * qJ(2);
t43 = t21 * t15 + t17;
t40 = t21 * t23;
t39 = t23 * rSges(6,3);
t38 = t24 * t19;
t37 = t24 * t22;
t36 = -pkin(1) - qJ(3);
t35 = t24 * pkin(1) + t21 * qJ(2);
t34 = qJ(5) * t20;
t33 = t24 * qJ(3) + t35;
t32 = -m(4) - m(5) + t51;
t31 = t24 * t46 + t33;
t30 = (-pkin(4) - t53) * t20;
t28 = t20 * rSges(5,1) + t23 * rSges(5,2);
t27 = -t20 * rSges(6,2) - t39;
t26 = g(2) * (pkin(4) * t40 + t21 * t34) + (pkin(4) * t23 + t34) * t54;
t5 = t19 * t40 - t37;
t4 = t22 * t40 + t38;
t3 = t21 * t22 + t23 * t38;
t2 = t21 * t19 - t23 * t37;
t1 = [-m(2) * (g(1) * (-t21 * rSges(2,1) - t24 * rSges(2,2)) + g(2) * (t24 * rSges(2,1) - t21 * rSges(2,2))) - m(3) * (g(1) * (t24 * rSges(3,3) + t17 + (rSges(3,2) - pkin(1)) * t21) + g(2) * (-t24 * rSges(3,2) + t21 * rSges(3,3) + t35)) - m(4) * (g(1) * (t24 * rSges(4,2) + t17) + g(2) * (t24 * rSges(4,3) + t33) + (g(1) * (-rSges(4,3) + t36) + g(2) * rSges(4,2)) * t21) - m(5) * (g(1) * t17 + g(2) * t33 + (g(1) * t44 + g(2) * t28) * t24 + (g(1) * (-t28 + t36) + g(2) * t44) * t21) - m(6) * (g(1) * t43 + g(2) * t31 + (g(1) * t45 + g(2) * (t27 - t15)) * t24 + (g(1) * (-t27 + t36 - t46) + g(2) * t45) * t21) - m(7) * ((-t3 * rSges(7,1) + t2 * rSges(7,2) + t31 + (t20 * t53 - t15) * t24 + t50 * t21) * g(2) + (t5 * rSges(7,1) + t4 * rSges(7,2) + t43 + t50 * t24 + (t30 + t36) * t21) * g(1)) (-m(3) + t32) * (g(1) * t21 - g(2) * t24) t32 * t29, m(5) * g(3) * t28 - m(6) * (g(3) * (t39 + t15 + (rSges(6,2) - pkin(4)) * t20) + t26) - m(7) * (g(3) * (t23 * t52 + t15 + t30) + t26) + t29 * (-m(5) * (rSges(5,1) * t23 - rSges(5,2) * t20) - m(6) * (-rSges(6,2) * t23 + rSges(6,3) * t20) - m(7) * (t52 * t20 + t23 * t53)) t51 * (-t29 * t23 + t47) -m(7) * (g(1) * (rSges(7,1) * t2 + rSges(7,2) * t3) + g(2) * (-rSges(7,1) * t4 + rSges(7,2) * t5) + (rSges(7,1) * t22 - rSges(7,2) * t19) * t47)];
taug  = t1(:);
