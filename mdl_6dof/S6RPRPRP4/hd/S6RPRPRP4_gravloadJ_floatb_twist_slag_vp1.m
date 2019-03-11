% Calculate Gravitation load on the joints for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:13
% EndTime: 2019-03-09 03:11:15
% DurationCPUTime: 0.61s
% Computational Cost: add. (352->106), mult. (400->141), div. (0->0), fcn. (375->8), ass. (0->44)
t23 = qJ(1) + pkin(9);
t17 = sin(t23);
t18 = cos(t23);
t63 = g(1) * t18 + g(2) * t17;
t56 = rSges(7,1) + pkin(5);
t25 = sin(qJ(3));
t28 = cos(qJ(3));
t34 = -t28 * rSges(5,2) + t25 * rSges(5,3);
t47 = rSges(7,3) + qJ(6);
t62 = -pkin(3) - pkin(8);
t61 = g(1) * t17;
t58 = g(3) * t28;
t26 = sin(qJ(1));
t57 = t26 * pkin(1);
t21 = t28 * pkin(3);
t55 = rSges(4,2) * t25;
t54 = t18 * t28;
t24 = sin(qJ(5));
t53 = t24 * t25;
t27 = cos(qJ(5));
t51 = t25 * t27;
t19 = t25 * qJ(4);
t49 = t19 + t21;
t46 = -m(5) - m(6) - m(7);
t45 = -rSges(7,2) + t62;
t44 = -rSges(6,3) + t62;
t29 = cos(qJ(1));
t22 = t29 * pkin(1);
t43 = t18 * pkin(2) + t17 * pkin(7) + t22;
t42 = t18 * pkin(7) - t57;
t41 = -pkin(2) - t19;
t40 = t18 * pkin(4) + t42;
t39 = t63 * qJ(4) * t28;
t38 = pkin(3) * t54 + t18 * t19 + t43;
t37 = rSges(4,1) * t28 - t55;
t35 = rSges(6,1) * t24 + rSges(6,2) * t27;
t32 = t17 * pkin(4) + pkin(8) * t54 + t38;
t31 = t56 * t24 - t47 * t27;
t30 = g(3) * (t28 * pkin(8) + t49) + t39;
t5 = -t17 * t53 + t18 * t27;
t4 = t17 * t51 + t18 * t24;
t3 = t17 * t27 + t18 * t53;
t2 = t17 * t24 - t18 * t51;
t1 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - rSges(2,2) * t29) + g(2) * (rSges(2,1) * t29 - t26 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t17 - rSges(3,2) * t18 - t57) + g(2) * (rSges(3,1) * t18 - rSges(3,2) * t17 + t22)) - m(4) * (g(1) * (t18 * rSges(4,3) + t42) + g(2) * (rSges(4,1) * t54 - t18 * t55 + t43) + (g(1) * (-pkin(2) - t37) + g(2) * rSges(4,3)) * t17) - m(5) * (g(1) * (t18 * rSges(5,1) + t42) + g(2) * (t34 * t18 + t38) + (g(1) * (-t34 + t41 - t21) + g(2) * rSges(5,1)) * t17) - m(6) * (g(1) * (t5 * rSges(6,1) - t4 * rSges(6,2) + t40) + g(2) * (rSges(6,1) * t3 - rSges(6,2) * t2 + rSges(6,3) * t54 + t32) + (t44 * t28 + t41) * t61) - m(7) * (g(1) * (t47 * t4 + t56 * t5 + t40) + g(2) * (rSges(7,2) * t54 + t47 * t2 + t56 * t3 + t32) + (t45 * t28 + t41) * t61) (-m(3) - m(4) + t46) * g(3), -m(4) * (g(3) * t37 + t63 * (-rSges(4,1) * t25 - rSges(4,2) * t28)) - m(5) * (g(3) * (t34 + t49) + t39 + t63 * (rSges(5,3) * t28 + (rSges(5,2) - pkin(3)) * t25)) - m(6) * ((g(3) * rSges(6,3) + t63 * t35) * t28 + (g(3) * t35 + t63 * t44) * t25 + t30) - m(7) * ((g(3) * rSges(7,2) + t63 * t31) * t28 + (g(3) * t31 + t63 * t45) * t25 + t30) t46 * (t63 * t25 - t58) -m(6) * (g(1) * (-rSges(6,1) * t2 - rSges(6,2) * t3) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t5)) - m(7) * (g(1) * (-t56 * t2 + t47 * t3) + g(2) * (t56 * t4 - t47 * t5)) + (-m(6) * (-rSges(6,1) * t27 + rSges(6,2) * t24) - m(7) * (-t47 * t24 - t56 * t27)) * t58, -m(7) * (g(1) * t2 - g(2) * t4 + t27 * t58)];
taug  = t1(:);
