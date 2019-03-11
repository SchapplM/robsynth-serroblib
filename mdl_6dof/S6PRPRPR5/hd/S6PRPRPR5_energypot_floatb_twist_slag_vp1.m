% Calculate potential energy for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRPR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:42:51
% EndTime: 2019-03-08 19:42:52
% DurationCPUTime: 0.43s
% Computational Cost: add. (353->130), mult. (523->159), div. (0->0), fcn. (603->12), ass. (0->47)
t34 = sin(pkin(11));
t65 = pkin(3) * t34;
t64 = pkin(9) + rSges(7,3);
t35 = sin(pkin(10));
t36 = sin(pkin(6));
t63 = t35 * t36;
t38 = cos(pkin(10));
t62 = t36 * t38;
t42 = sin(qJ(2));
t61 = t36 * t42;
t44 = cos(qJ(2));
t60 = t36 * t44;
t39 = cos(pkin(6));
t59 = t39 * t42;
t58 = t39 * t44;
t57 = rSges(6,3) + qJ(5);
t56 = qJ(3) + rSges(4,3);
t55 = t35 * pkin(1) + r_base(2);
t54 = t34 * t63;
t53 = qJ(1) + r_base(3);
t52 = t38 * pkin(1) + pkin(7) * t63 + r_base(1);
t51 = t39 * pkin(7) + t53;
t16 = t35 * t58 + t38 * t42;
t17 = -t35 * t59 + t38 * t44;
t37 = cos(pkin(11));
t27 = t37 * pkin(3) + pkin(2);
t40 = -pkin(8) - qJ(3);
t50 = pkin(3) * t54 - t16 * t40 + t17 * t27 + t52;
t49 = t27 * t61 + t39 * t65 + t40 * t60 + t51;
t33 = pkin(11) + qJ(4);
t28 = sin(t33);
t29 = cos(t33);
t6 = t17 * t29 + t28 * t63;
t48 = t6 * pkin(4) + t50;
t11 = t39 * t28 + t29 * t61;
t47 = t11 * pkin(4) + t49;
t14 = t35 * t42 - t38 * t58;
t15 = t35 * t44 + t38 * t59;
t46 = t15 * t27 + (-pkin(7) - t65) * t62 - t14 * t40 + t55;
t4 = t15 * t29 - t28 * t62;
t45 = t4 * pkin(4) + t46;
t43 = cos(qJ(6));
t41 = sin(qJ(6));
t10 = t28 * t61 - t39 * t29;
t5 = t17 * t28 - t29 * t63;
t3 = t15 * t28 + t29 * t62;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t38 * rSges(2,1) - t35 * rSges(2,2) + r_base(1)) + g(2) * (t35 * rSges(2,1) + t38 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t53)) - m(3) * (g(1) * (t17 * rSges(3,1) - t16 * rSges(3,2) + t52) + g(2) * (t15 * rSges(3,1) - t14 * rSges(3,2) + t55) + g(3) * (t39 * rSges(3,3) + t51) + (g(1) * rSges(3,3) * t35 + g(3) * (rSges(3,1) * t42 + rSges(3,2) * t44) + g(2) * (-rSges(3,3) - pkin(7)) * t38) * t36) - m(4) * (g(1) * (t17 * pkin(2) + (t17 * t37 + t54) * rSges(4,1) + (-t17 * t34 + t37 * t63) * rSges(4,2) + t56 * t16 + t52) + g(2) * (t15 * pkin(2) - pkin(7) * t62 + (t15 * t37 - t34 * t62) * rSges(4,1) + (-t15 * t34 - t37 * t62) * rSges(4,2) + t56 * t14 + t55) + g(3) * ((t34 * rSges(4,1) + t37 * rSges(4,2)) * t39 + (-t56 * t44 + (t37 * rSges(4,1) - t34 * rSges(4,2) + pkin(2)) * t42) * t36 + t51)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t16 * rSges(5,3) + t50) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t14 * rSges(5,3) + t46) + g(3) * (t11 * rSges(5,1) - t10 * rSges(5,2) - rSges(5,3) * t60 + t49)) - m(6) * (g(1) * (t16 * rSges(6,1) - t6 * rSges(6,2) + t57 * t5 + t48) + g(2) * (t14 * rSges(6,1) - t4 * rSges(6,2) + t57 * t3 + t45) + g(3) * (-rSges(6,1) * t60 - t11 * rSges(6,2) + t57 * t10 + t47)) - m(7) * (g(1) * (t16 * pkin(5) + t5 * qJ(5) + (t16 * t43 + t5 * t41) * rSges(7,1) + (-t16 * t41 + t5 * t43) * rSges(7,2) + t64 * t6 + t48) + g(2) * (t14 * pkin(5) + t3 * qJ(5) + (t14 * t43 + t3 * t41) * rSges(7,1) + (-t14 * t41 + t3 * t43) * rSges(7,2) + t64 * t4 + t45) + g(3) * (-pkin(5) * t60 + t10 * qJ(5) + (t10 * t41 - t43 * t60) * rSges(7,1) + (t10 * t43 + t41 * t60) * rSges(7,2) + t64 * t11 + t47));
U  = t1;
