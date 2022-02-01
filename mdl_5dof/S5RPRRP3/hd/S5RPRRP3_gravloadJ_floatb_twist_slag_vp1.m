% Calculate Gravitation load on the joints for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:35
% EndTime: 2022-01-23 09:29:36
% DurationCPUTime: 0.31s
% Computational Cost: add. (204->60), mult. (174->73), div. (0->0), fcn. (134->8), ass. (0->31)
t44 = rSges(6,1) + pkin(4);
t15 = qJ(3) + qJ(4);
t10 = cos(t15);
t9 = sin(t15);
t31 = t10 * rSges(5,1) - t9 * rSges(5,2);
t43 = -rSges(6,2) * t10 - t44 * t9;
t14 = qJ(1) + pkin(8);
t7 = sin(t14);
t8 = cos(t14);
t42 = g(1) * t8 + g(2) * t7;
t30 = -t9 * rSges(6,2) + t44 * t10;
t20 = -pkin(7) - pkin(6);
t16 = sin(qJ(3));
t38 = t16 * pkin(3);
t17 = sin(qJ(1));
t37 = t17 * pkin(1);
t34 = rSges(4,3) + pkin(6);
t18 = cos(qJ(3));
t11 = t18 * pkin(3);
t6 = t11 + pkin(2);
t33 = rSges(5,3) - t20;
t32 = rSges(6,3) + qJ(5) - t20;
t29 = -rSges(5,1) * t9 - rSges(5,2) * t10;
t27 = t18 * rSges(4,1) - t16 * rSges(4,2);
t19 = cos(qJ(1));
t12 = t19 * pkin(1);
t26 = -g(1) * t37 + g(2) * t12;
t25 = t6 + t31;
t24 = t6 + t30;
t23 = pkin(2) + t27;
t1 = [-m(2) * (g(1) * (-t17 * rSges(2,1) - t19 * rSges(2,2)) + g(2) * (t19 * rSges(2,1) - t17 * rSges(2,2))) - m(3) * (g(1) * (-t7 * rSges(3,1) - t8 * rSges(3,2) - t37) + g(2) * (t8 * rSges(3,1) - t7 * rSges(3,2) + t12)) - m(4) * ((g(1) * t34 + g(2) * t23) * t8 + (-g(1) * t23 + g(2) * t34) * t7 + t26) - m(5) * ((g(1) * t33 + g(2) * t25) * t8 + (-g(1) * t25 + g(2) * t33) * t7 + t26) - m(6) * ((g(1) * t32 + g(2) * t24) * t8 + (-g(1) * t24 + g(2) * t32) * t7 + t26), (-m(3) - m(4) - m(5) - m(6)) * g(3), (-m(4) * t27 - m(5) * (t11 + t31) - m(6) * (t11 + t30)) * g(3) + t42 * (-m(4) * (-rSges(4,1) * t16 - rSges(4,2) * t18) - m(5) * (t29 - t38) - m(6) * (-t38 + t43)), (-m(5) * t31 - m(6) * t30) * g(3) + t42 * (-m(5) * t29 - m(6) * t43), -m(6) * (g(1) * t7 - g(2) * t8)];
taug = t1(:);
