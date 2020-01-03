% Calculate Gravitation load on the joints for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:36
% EndTime: 2019-12-31 20:57:37
% DurationCPUTime: 0.48s
% Computational Cost: add. (247->93), mult. (281->115), div. (0->0), fcn. (235->6), ass. (0->49)
t60 = rSges(6,1) + pkin(4);
t23 = qJ(2) + qJ(3);
t20 = sin(t23);
t21 = cos(t23);
t32 = t21 * rSges(5,1) + t20 * rSges(5,3);
t16 = t21 * rSges(4,1);
t33 = -t20 * rSges(4,2) + t16;
t27 = cos(qJ(1));
t25 = sin(qJ(1));
t54 = g(2) * t25;
t59 = g(1) * t27 + t54;
t13 = t20 * rSges(6,2);
t58 = t60 * t21 + t13;
t57 = t59 * t20;
t24 = sin(qJ(2));
t56 = pkin(2) * t24;
t18 = t21 * pkin(3);
t53 = -rSges(5,1) - pkin(3);
t52 = rSges(3,3) + pkin(6);
t50 = t21 * t25;
t49 = t21 * t27;
t28 = -pkin(7) - pkin(6);
t48 = rSges(5,2) - t28;
t47 = rSges(4,3) - t28;
t11 = t20 * qJ(4);
t46 = t11 + t18;
t45 = qJ(4) * t21;
t44 = -pkin(3) - t60;
t43 = -rSges(6,3) - qJ(5) - t28;
t26 = cos(qJ(2));
t22 = t26 * pkin(2);
t19 = t22 + pkin(1);
t42 = -t19 - t11;
t41 = t46 + t32;
t40 = g(1) * t44;
t5 = t27 * t19;
t39 = g(2) * (pkin(3) * t49 + t27 * t11 + t5);
t2 = t25 * t45;
t38 = -t25 * t56 + t2;
t4 = t27 * t45;
t37 = -t27 * t56 + t4;
t36 = t46 + t58;
t35 = t26 * rSges(3,1) - t24 * rSges(3,2);
t31 = pkin(1) + t35;
t9 = rSges(6,2) * t49;
t8 = rSges(5,3) * t49;
t7 = rSges(6,2) * t50;
t6 = rSges(5,3) * t50;
t1 = [-m(2) * (g(1) * (-t25 * rSges(2,1) - t27 * rSges(2,2)) + g(2) * (t27 * rSges(2,1) - t25 * rSges(2,2))) - m(3) * ((g(1) * t52 + g(2) * t31) * t27 + (-g(1) * t31 + g(2) * t52) * t25) - m(4) * (g(2) * t5 + (g(1) * t47 + g(2) * t33) * t27 + (g(1) * (-t19 - t33) + g(2) * t47) * t25) - m(5) * (t39 + (g(1) * t48 + g(2) * t32) * t27 + (g(1) * (-t32 + t42 - t18) + g(2) * t48) * t25) - m(6) * (t39 + (g(1) * t43 + g(2) * t58) * t27 + (g(1) * (t42 - t13) + g(2) * t43 + t21 * t40) * t25), -m(3) * (g(3) * t35 + t59 * (-rSges(3,1) * t24 - rSges(3,2) * t26)) - m(4) * (g(3) * (t22 + t33) + t59 * (-rSges(4,1) * t20 - rSges(4,2) * t21 - t56)) - m(5) * (g(1) * (t37 + t8) + g(2) * (t38 + t6) + g(3) * (t22 + t41) + t53 * t57) - m(6) * (g(1) * (t37 + t9) + g(2) * (t38 + t7) + g(3) * (t22 + t36) + (t27 * t40 + t44 * t54) * t20), -m(4) * (g(3) * t16 + (-g(1) * t49 - g(2) * t50) * rSges(4,2)) - m(5) * (g(1) * (t4 + t8) + g(2) * (t2 + t6) + g(3) * t41) - m(6) * (g(1) * (t4 + t9) + g(2) * (t2 + t7) + g(3) * t36) + (m(4) * g(3) * rSges(4,2) + t59 * (m(4) * rSges(4,1) - m(5) * t53 - m(6) * t44)) * t20, (-m(5) - m(6)) * (-g(3) * t21 + t57), -m(6) * (-g(1) * t25 + g(2) * t27)];
taug = t1(:);
