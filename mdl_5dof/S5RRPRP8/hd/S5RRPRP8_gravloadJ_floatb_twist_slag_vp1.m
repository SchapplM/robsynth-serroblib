% Calculate Gravitation load on the joints for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:21
% EndTime: 2019-12-31 20:03:23
% DurationCPUTime: 0.57s
% Computational Cost: add. (181->107), mult. (379->141), div. (0->0), fcn. (368->6), ass. (0->45)
t63 = -rSges(5,3) - pkin(7);
t25 = sin(qJ(2));
t18 = t25 * qJ(3);
t28 = cos(qJ(2));
t46 = t28 * pkin(2) + t18;
t26 = sin(qJ(1));
t29 = cos(qJ(1));
t62 = g(1) * t29 + g(2) * t26;
t61 = t62 * t25;
t57 = t28 * pkin(3);
t24 = sin(qJ(4));
t54 = t25 * t24;
t27 = cos(qJ(4));
t53 = t25 * t27;
t52 = t25 * t29;
t51 = t28 * rSges(4,1);
t17 = t27 * pkin(4) + pkin(3);
t50 = t28 * t17;
t49 = t28 * t24;
t48 = t28 * t29;
t47 = -rSges(6,3) - qJ(5) - pkin(7);
t45 = t29 * pkin(1) + t26 * pkin(6);
t44 = qJ(3) * t28;
t43 = t26 * t49;
t42 = t24 * t48;
t41 = pkin(2) * t48 + t29 * t18 + t45;
t2 = -t26 * t53 + t43;
t32 = t28 * t27 + t54;
t3 = t32 * t26;
t40 = -t2 * rSges(5,1) - t3 * rSges(5,2);
t4 = -t27 * t52 + t42;
t5 = t32 * t29;
t39 = -t4 * rSges(5,1) - t5 * rSges(5,2);
t7 = -t49 + t53;
t38 = -rSges(5,1) * t32 - t7 * rSges(5,2);
t37 = -t2 * rSges(6,1) - t3 * rSges(6,2);
t36 = -t4 * rSges(6,1) - t5 * rSges(6,2);
t35 = -rSges(6,1) * t32 - t7 * rSges(6,2);
t34 = t28 * rSges(3,1) - t25 * rSges(3,2);
t31 = pkin(4) * t54 + t50;
t30 = -pkin(1) - t46;
t21 = t29 * pkin(6);
t12 = t29 * t44;
t10 = t26 * t44;
t1 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - t29 * rSges(2,2)) + g(2) * (t29 * rSges(2,1) - t26 * rSges(2,2))) - m(3) * (g(1) * (t29 * rSges(3,3) + t21) + g(2) * (rSges(3,1) * t48 - rSges(3,2) * t52 + t45) + (g(1) * (-pkin(1) - t34) + g(2) * rSges(3,3)) * t26) - m(4) * (g(1) * (t29 * rSges(4,2) + t21) + g(2) * (rSges(4,1) * t48 + rSges(4,3) * t52 + t41) + (g(1) * (-t25 * rSges(4,3) + t30 - t51) + g(2) * rSges(4,2)) * t26) - m(5) * (g(1) * (-t3 * rSges(5,1) + t2 * rSges(5,2) + t63 * t29 + t21) + g(2) * (t5 * rSges(5,1) - t4 * rSges(5,2) + pkin(3) * t48 + t41) + (g(1) * (t30 - t57) + g(2) * t63) * t26) - m(6) * (g(1) * (-t3 * rSges(6,1) + t2 * rSges(6,2) + t21) + g(2) * (t5 * rSges(6,1) - t4 * rSges(6,2) + t41) + (g(1) * t47 + g(2) * t31) * t29 + (g(1) * (t30 - t31) + g(2) * t47) * t26), -m(3) * (g(3) * t34 + t62 * (-rSges(3,1) * t25 - rSges(3,2) * t28)) - m(4) * (g(1) * (rSges(4,3) * t48 + t12) + g(2) * (t26 * t28 * rSges(4,3) + t10) + g(3) * (t46 + t51) + (g(3) * rSges(4,3) + t62 * (-rSges(4,1) - pkin(2))) * t25) - m(5) * (g(1) * (t12 - t39) + g(2) * (t10 - t40) + g(3) * (-t38 + t46 + t57) + (-pkin(2) - pkin(3)) * t61) - m(6) * (g(1) * (pkin(4) * t42 + t12 - t36) + g(2) * (pkin(4) * t43 + t10 - t37) + g(3) * (-t35 + t46 + t50) + (g(3) * pkin(4) * t24 + t62 * (-pkin(2) - t17)) * t25), (-m(4) - m(5) - m(6)) * (-g(3) * t28 + t61), -m(5) * (g(1) * t39 + g(2) * t40 + g(3) * t38) + (-g(1) * t36 - g(2) * t37 - g(3) * t35 - (-g(3) * t32 + t62 * t7) * pkin(4)) * m(6), -m(6) * (-g(1) * t26 + g(2) * t29)];
taug = t1(:);
