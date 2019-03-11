% Calculate Gravitation load on the joints for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:24:26
% EndTime: 2019-03-09 03:24:28
% DurationCPUTime: 0.67s
% Computational Cost: add. (290->114), mult. (390->153), div. (0->0), fcn. (365->8), ass. (0->49)
t24 = sin(qJ(5));
t27 = cos(qJ(5));
t42 = rSges(7,3) + qJ(6);
t52 = rSges(7,1) + pkin(5);
t60 = t42 * t24 + t52 * t27;
t65 = -pkin(4) - t60;
t29 = cos(qJ(1));
t55 = g(2) * t29;
t26 = sin(qJ(1));
t58 = g(1) * t26;
t63 = -t55 + t58;
t62 = g(1) * t29 + g(2) * t26;
t22 = qJ(3) + pkin(9);
t18 = cos(t22);
t61 = t62 * t18;
t28 = cos(qJ(3));
t59 = pkin(3) * t28;
t54 = g(3) * t18;
t25 = sin(qJ(3));
t53 = t25 * pkin(3);
t51 = -rSges(7,2) - pkin(8);
t50 = rSges(4,3) + pkin(7);
t49 = -rSges(6,3) - pkin(8);
t17 = sin(t22);
t48 = t17 * t26;
t47 = t24 * t29;
t46 = t26 * t24;
t45 = t26 * t27;
t44 = t29 * t27;
t43 = t29 * pkin(1) + t26 * qJ(2);
t41 = -m(5) - m(6) - m(7);
t12 = t26 * t59;
t40 = t26 * t18 * pkin(4) + pkin(8) * t48 + t12;
t39 = t26 * t53 + t43;
t20 = t29 * qJ(2);
t23 = -qJ(4) - pkin(7);
t38 = t26 * t23 + t29 * t53 + t20;
t37 = t18 * pkin(8) - t53;
t35 = rSges(4,1) * t25 + rSges(4,2) * t28;
t34 = rSges(5,1) * t18 - rSges(5,2) * t17;
t33 = rSges(5,1) * t17 + rSges(5,2) * t18;
t32 = -rSges(6,1) * t27 + rSges(6,2) * t24 - pkin(4);
t31 = t29 * t17 * pkin(4) - t26 * pkin(1) + t38;
t30 = pkin(4) * t48 - t23 * t29 + t39;
t4 = t17 * t44 - t46;
t3 = t17 * t47 + t45;
t2 = t17 * t45 + t47;
t1 = t17 * t46 - t44;
t5 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - rSges(2,2) * t29) + g(2) * (rSges(2,1) * t29 - t26 * rSges(2,2))) - m(3) * (g(1) * (t29 * rSges(3,3) + t20 + (rSges(3,2) - pkin(1)) * t26) + g(2) * (-rSges(3,2) * t29 + t26 * rSges(3,3) + t43)) - m(4) * (g(1) * t20 + g(2) * t43 + (g(1) * t35 + g(2) * t50) * t29 + (g(1) * (-pkin(1) - t50) + g(2) * t35) * t26) - m(5) * (g(1) * t38 + g(2) * t39 + (g(1) * t33 + g(2) * (rSges(5,3) - t23)) * t29 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t33) * t26) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t31) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t30) + t49 * t61) - m(7) * (g(1) * (t42 * t3 + t52 * t4 + t31) + g(2) * (t42 * t1 + t52 * t2 + t30) + t51 * t61) (-m(3) - m(4) + t41) * t63, -m(4) * (-g(3) * t35 + t63 * (rSges(4,1) * t28 - rSges(4,2) * t25)) - m(5) * (g(1) * (t34 * t26 + t12) + g(3) * (-t33 - t53) + (-t34 - t59) * t55) - m(6) * (g(1) * (t40 + (rSges(6,1) * t45 - rSges(6,2) * t46) * t18) + g(3) * (t18 * rSges(6,3) + t37) + (rSges(6,3) * t58 + g(3) * t32) * t17 + (t49 * t17 + t32 * t18 - t59) * t55) - m(7) * (g(1) * t40 + g(3) * t37 + (g(3) * rSges(7,2) + t60 * t58) * t18 + (rSges(7,2) * t58 + g(3) * t65) * t17 + (t51 * t17 + t65 * t18 - t59) * t55) t41 * t62, -m(6) * (g(1) * (-rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (rSges(6,1) * t3 + rSges(6,2) * t4)) - m(7) * (g(1) * (-t52 * t1 + t42 * t2) + g(2) * (t52 * t3 - t42 * t4)) + (-m(6) * (-rSges(6,1) * t24 - rSges(6,2) * t27) - m(7) * (-t52 * t24 + t42 * t27)) * t54, -m(7) * (g(1) * t1 - g(2) * t3 + t24 * t54)];
taug  = t5(:);
