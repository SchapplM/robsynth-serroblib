% Calculate Gravitation load on the joints for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:22
% EndTime: 2020-01-03 12:07:23
% DurationCPUTime: 0.20s
% Computational Cost: add. (308->60), mult. (162->72), div. (0->0), fcn. (120->10), ass. (0->37)
t54 = rSges(6,3) + pkin(8);
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t53 = rSges(6,1) * t33 - rSges(6,2) * t31;
t52 = pkin(4) + t53;
t30 = qJ(1) + qJ(2);
t27 = qJ(3) + t30;
t23 = sin(t27);
t24 = cos(t27);
t49 = t23 * rSges(4,1) + t24 * rSges(4,2);
t25 = sin(t30);
t26 = cos(t30);
t48 = t25 * rSges(3,1) + t26 * rSges(3,2);
t22 = pkin(9) + t27;
t13 = sin(t22);
t14 = cos(t22);
t15 = pkin(3) * t23;
t47 = t13 * rSges(5,1) + t14 * rSges(5,2) + t15;
t20 = pkin(2) * t25;
t46 = t20 + t49;
t45 = t20 + t47;
t44 = t26 * rSges(3,1) - rSges(3,2) * t25;
t43 = t24 * rSges(4,1) - rSges(4,2) * t23;
t16 = pkin(3) * t24;
t42 = t14 * rSges(5,1) - rSges(5,2) * t13 + t16;
t21 = pkin(2) * t26;
t41 = t21 + t43;
t39 = t21 + t42;
t38 = t54 * t13 + t52 * t14 + t16;
t37 = t21 + t38;
t36 = t52 * t13 - t54 * t14 + t15;
t35 = t20 + t36;
t34 = cos(qJ(1));
t32 = sin(qJ(1));
t29 = t34 * pkin(1);
t28 = t32 * pkin(1);
t1 = [-m(2) * (g(2) * (rSges(2,1) * t34 - t32 * rSges(2,2)) + g(3) * (t32 * rSges(2,1) + rSges(2,2) * t34)) - m(3) * (g(2) * (t29 + t44) + g(3) * (t28 + t48)) - m(4) * (g(2) * (t29 + t41) + g(3) * (t28 + t46)) - m(5) * (g(2) * (t29 + t39) + g(3) * (t28 + t45)) - m(6) * (g(2) * (t29 + t37) + g(3) * (t28 + t35)), -m(3) * (g(2) * t44 + g(3) * t48) - m(4) * (g(2) * t41 + g(3) * t46) - m(5) * (g(2) * t39 + g(3) * t45) - m(6) * (g(2) * t37 + g(3) * t35), -m(4) * (g(2) * t43 + g(3) * t49) - m(5) * (g(2) * t42 + g(3) * t47) - m(6) * (g(2) * t38 + g(3) * t36), (-m(5) - m(6)) * g(1), -m(6) * (g(1) * t53 + (-g(2) * t13 + g(3) * t14) * (rSges(6,1) * t31 + rSges(6,2) * t33))];
taug = t1(:);
