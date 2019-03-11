% Calculate Gravitation load on the joints for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:19
% EndTime: 2019-03-09 01:55:20
% DurationCPUTime: 0.38s
% Computational Cost: add. (225->101), mult. (279->134), div. (0->0), fcn. (242->8), ass. (0->44)
t59 = rSges(7,3) + pkin(8);
t21 = pkin(9) + qJ(4);
t17 = cos(t21);
t14 = t17 * qJ(5);
t58 = -t17 * rSges(6,3) - t14;
t25 = sin(qJ(6));
t27 = cos(qJ(6));
t31 = rSges(7,1) * t25 + rSges(7,2) * t27;
t57 = -m(6) - m(7);
t22 = sin(pkin(9));
t56 = pkin(3) * t22;
t26 = sin(qJ(1));
t55 = g(1) * t26;
t28 = cos(qJ(1));
t54 = g(2) * t28;
t16 = sin(t21);
t53 = g(3) * t16;
t52 = rSges(6,2) - pkin(4);
t49 = t16 * t26;
t47 = t26 * t25;
t46 = t26 * t27;
t45 = t28 * t25;
t44 = t28 * t27;
t43 = t28 * pkin(1) + t26 * qJ(2);
t42 = rSges(4,3) + qJ(3);
t41 = -pkin(4) - t59;
t40 = t26 * t56 + t43;
t19 = t28 * qJ(2);
t24 = -pkin(7) - qJ(3);
t39 = t26 * t24 + t28 * t56 + t19;
t38 = g(1) * (t26 * t17 * pkin(4) + qJ(5) * t49);
t37 = -m(4) - m(5) + t57;
t36 = pkin(4) * t49 + t40;
t35 = t28 * t16 * pkin(4) + t39;
t34 = -t54 + t55;
t33 = rSges(4,1) * t22 + rSges(4,2) * cos(pkin(9));
t32 = t16 * rSges(5,1) + t17 * rSges(5,2);
t30 = t59 * t16 - t14;
t29 = -t16 * rSges(6,2) + t58;
t5 = -t17 * t47 + t44;
t4 = -t17 * t46 - t45;
t3 = -t17 * t45 - t46;
t2 = -t17 * t44 + t47;
t1 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - t28 * rSges(2,2)) + g(2) * (t28 * rSges(2,1) - t26 * rSges(2,2))) - m(3) * (g(1) * (t28 * rSges(3,3) + t19 + (rSges(3,2) - pkin(1)) * t26) + g(2) * (-t28 * rSges(3,2) + t26 * rSges(3,3) + t43)) - m(4) * (g(1) * t19 + g(2) * t43 + (g(1) * t33 + g(2) * t42) * t28 + (g(1) * (-pkin(1) - t42) + g(2) * t33) * t26) - m(5) * (g(1) * t39 + g(2) * t40 + (g(1) * t32 + g(2) * (rSges(5,3) - t24)) * t28 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t32) * t26) - m(6) * (g(1) * t35 + g(2) * t36 + (g(1) * t29 + g(2) * (rSges(6,1) - t24)) * t28 + (g(1) * (-rSges(6,1) - pkin(1)) + g(2) * t29) * t26) - m(7) * (g(1) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t35) + g(2) * (t5 * rSges(7,1) + t4 * rSges(7,2) + t36) + (g(1) * t30 + g(2) * (pkin(5) - t24)) * t28 + (g(1) * (-pkin(1) - pkin(5)) + g(2) * t30) * t26) (-m(3) + t37) * t34, t37 * (g(1) * t28 + g(2) * t26) m(5) * g(3) * t32 - m(6) * (t38 + g(3) * (t16 * t52 - t58)) - m(7) * (t38 + g(3) * (t16 * t41 + t17 * t31 + t14)) + (-m(5) * (rSges(5,1) * t17 - rSges(5,2) * t16) - m(6) * (-rSges(6,2) * t17 + rSges(6,3) * t16) - m(7) * (t31 * t16 + t59 * t17)) * t55 + ((m(5) * rSges(5,1) - m(6) * t52 - m(7) * t41) * t17 + (-m(5) * rSges(5,2) - m(6) * (-rSges(6,3) - qJ(5)) - m(7) * (-qJ(5) - t31)) * t16) * t54, t57 * (-t17 * t34 + t53) -m(7) * (g(1) * (t4 * rSges(7,1) - t5 * rSges(7,2)) + g(2) * (-t2 * rSges(7,1) + t3 * rSges(7,2)) + (rSges(7,1) * t27 - rSges(7,2) * t25) * t53)];
taug  = t1(:);
