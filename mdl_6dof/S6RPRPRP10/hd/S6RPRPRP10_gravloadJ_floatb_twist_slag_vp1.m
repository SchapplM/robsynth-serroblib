% Calculate Gravitation load on the joints for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:30:58
% EndTime: 2019-03-09 03:31:00
% DurationCPUTime: 0.65s
% Computational Cost: add. (205->120), mult. (412->159), div. (0->0), fcn. (387->6), ass. (0->45)
t22 = sin(qJ(5));
t25 = cos(qJ(5));
t42 = rSges(7,3) + qJ(6);
t54 = rSges(7,1) + pkin(5);
t60 = t54 * t22 - t42 * t25;
t26 = cos(qJ(3));
t17 = t26 * qJ(4);
t61 = -t26 * rSges(5,3) - t17;
t27 = cos(qJ(1));
t56 = g(2) * t27;
t24 = sin(qJ(1));
t57 = g(1) * t24;
t34 = -t56 + t57;
t59 = -pkin(1) - pkin(7);
t58 = -pkin(3) - pkin(8);
t23 = sin(qJ(3));
t55 = g(3) * t23;
t53 = rSges(5,2) - pkin(3);
t48 = t24 * t26;
t50 = t23 * t24;
t52 = pkin(3) * t48 + qJ(4) * t50;
t51 = rSges(6,1) * t22;
t49 = t23 * t27;
t47 = t25 * t26;
t45 = t26 * t27;
t18 = t27 * qJ(2);
t44 = pkin(3) * t49 + t18;
t43 = t27 * pkin(1) + t24 * qJ(2);
t41 = -m(5) - m(6) - m(7);
t40 = -rSges(7,2) + t58;
t39 = -rSges(6,3) + t58;
t38 = pkin(8) * t48 + t52;
t37 = t27 * pkin(7) + t43;
t36 = g(1) * (-pkin(4) + t59);
t35 = pkin(3) * t50 + t37;
t32 = rSges(4,1) * t23 + rSges(4,2) * t26;
t31 = rSges(6,2) * t25 + t51;
t30 = t27 * pkin(4) + pkin(8) * t50 + t35;
t29 = pkin(8) * t49 - t27 * t17 + t44;
t28 = -rSges(5,2) * t23 + t61;
t5 = -t22 * t48 + t25 * t27;
t4 = t22 * t27 + t24 * t47;
t3 = t22 * t45 + t24 * t25;
t2 = t22 * t24 - t25 * t45;
t1 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - rSges(2,2) * t27) + g(2) * (rSges(2,1) * t27 - t24 * rSges(2,2))) - m(3) * (g(1) * (t27 * rSges(3,3) + t18 + (rSges(3,2) - pkin(1)) * t24) + g(2) * (-rSges(3,2) * t27 + t24 * rSges(3,3) + t43)) - m(4) * (g(1) * (rSges(4,1) * t49 + rSges(4,2) * t45 + t18) + g(2) * (rSges(4,3) * t27 + t37) + (g(1) * (-rSges(4,3) + t59) + g(2) * t32) * t24) - m(5) * (g(1) * t44 + g(2) * t35 + (g(2) * rSges(5,1) + g(1) * t28) * t27 + (g(1) * (-rSges(5,1) + t59) + g(2) * t28) * t24) - m(6) * (g(1) * (-t3 * rSges(6,1) + t2 * rSges(6,2) + rSges(6,3) * t49 + t29) + g(2) * (rSges(6,1) * t5 - rSges(6,2) * t4 + t30) + (t36 + g(2) * (rSges(6,3) * t23 - t17)) * t24) - m(7) * (g(1) * (rSges(7,2) * t49 - t42 * t2 - t3 * t54 + t29) + g(2) * (t42 * t4 + t54 * t5 + t30) + (t36 + g(2) * (rSges(7,2) * t23 - t17)) * t24) (-m(3) - m(4) + t41) * t34, -m(4) * (-g(3) * t32 + t34 * (rSges(4,1) * t26 - rSges(4,2) * t23)) - m(5) * (g(1) * ((-rSges(5,2) * t26 + rSges(5,3) * t23) * t24 + t52) + g(3) * (t53 * t23 - t61) + (t53 * t26 + (-rSges(5,3) - qJ(4)) * t23) * t56) - m(6) * (g(1) * (rSges(6,3) * t48 + t38) + g(3) * (rSges(6,2) * t47 + t26 * t51 + t17) + (g(3) * t39 + t31 * t57) * t23 + (t39 * t26 + (-qJ(4) - t31) * t23) * t56) - m(7) * (g(1) * t38 + g(3) * t17 + (rSges(7,2) * t57 + g(3) * t60 + t40 * t56) * t26 + (g(3) * t40 + t60 * t57 + (-qJ(4) - t60) * t56) * t23) t41 * (-t34 * t26 + t55) -m(6) * (g(1) * (-rSges(6,1) * t4 - rSges(6,2) * t5) + g(2) * (-rSges(6,1) * t2 - rSges(6,2) * t3)) - m(7) * (g(1) * (-t54 * t4 + t42 * t5) + g(2) * (-t54 * t2 + t42 * t3)) + (-m(6) * (rSges(6,1) * t25 - rSges(6,2) * t22) - m(7) * (t42 * t22 + t25 * t54)) * t55, -m(7) * (g(1) * t4 + g(2) * t2 - t25 * t55)];
taug  = t1(:);
