% Calculate Gravitation load on the joints for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:25:51
% EndTime: 2019-03-29 15:25:52
% DurationCPUTime: 0.32s
% Computational Cost: add. (272->78), mult. (269->117), div. (0->0), fcn. (230->10), ass. (0->47)
t39 = cos(qJ(3));
t31 = t39 * pkin(2);
t33 = qJ(3) + qJ(4);
t27 = sin(t33);
t29 = cos(t33);
t50 = rSges(5,1) * t29 - rSges(5,2) * t27;
t72 = t31 + t50;
t36 = sin(qJ(3));
t71 = rSges(4,1) * t39 - rSges(4,2) * t36;
t34 = qJ(1) + qJ(2);
t28 = sin(t34);
t30 = cos(t34);
t70 = g(1) * t30 + g(2) * t28;
t69 = pkin(2) * t36;
t37 = sin(qJ(1));
t66 = t37 * pkin(1);
t35 = sin(qJ(5));
t53 = rSges(6,2) * t27 * t35;
t65 = (rSges(6,3) * t29 + t53) * t28;
t20 = t27 * rSges(6,3);
t63 = t27 * t30;
t62 = t29 * t30;
t61 = t29 * t35;
t38 = cos(qJ(5));
t60 = t29 * t38;
t59 = t30 * t35;
t58 = t30 * t38;
t55 = rSges(6,3) * t62 + t30 * t53;
t54 = rSges(6,1) * t27 * t38;
t19 = t30 * t31;
t7 = t28 * t38 - t29 * t59;
t8 = t28 * t35 + t29 * t58;
t52 = rSges(6,1) * t8 + rSges(6,2) * t7 + rSges(6,3) * t63 + t19;
t51 = rSges(3,1) * t30 - rSges(3,2) * t28;
t49 = -rSges(3,1) * t28 - rSges(3,2) * t30;
t48 = -rSges(5,1) * t27 - rSges(5,2) * t29;
t47 = t30 * rSges(4,3) - t28 * t71;
t46 = t28 * rSges(4,3) + t30 * t71;
t45 = rSges(6,1) * t60 - rSges(6,2) * t61 + t20;
t44 = rSges(5,1) * t62 - rSges(5,2) * t63 + rSges(5,3) * t28 + t19;
t5 = t28 * t61 + t58;
t6 = -t28 * t60 + t59;
t43 = t5 * rSges(6,2) + t6 * rSges(6,1) + (-t20 - t31) * t28;
t42 = t30 * rSges(5,3) - t28 * t72;
t40 = cos(qJ(1));
t32 = t40 * pkin(1);
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t37 - rSges(2,2) * t40) + g(2) * (rSges(2,1) * t40 - rSges(2,2) * t37)) - m(3) * (g(1) * (t49 - t66) + g(2) * (t32 + t51)) - m(4) * (g(1) * (t47 - t66) + g(2) * (t32 + t46)) - m(5) * (g(1) * (t42 - t66) + g(2) * (t32 + t44)) - m(6) * (g(1) * (t43 - t66) + g(2) * (t32 + t52)) -m(3) * (g(1) * t49 + g(2) * t51) - m(4) * (g(1) * t47 + g(2) * t46) - m(5) * (g(1) * t42 + g(2) * t44) - m(6) * (g(1) * t43 + g(2) * t52) -m(6) * (g(1) * t55 + g(2) * t65) + (-m(4) * t71 - m(5) * t72 - m(6) * (t31 + t45)) * g(3) + t70 * (-m(4) * (-rSges(4,1) * t36 - rSges(4,2) * t39) - m(5) * (t48 - t69) - m(6) * (-t54 - t69)) -m(5) * (g(3) * t50 + t48 * t70) - m(6) * (g(1) * (-t30 * t54 + t55) + g(2) * (-t28 * t54 + t65) + g(3) * t45) -m(6) * (g(1) * (rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (-rSges(6,1) * t5 + rSges(6,2) * t6) + g(3) * (-rSges(6,1) * t35 - rSges(6,2) * t38) * t27)];
taug  = t1(:);
