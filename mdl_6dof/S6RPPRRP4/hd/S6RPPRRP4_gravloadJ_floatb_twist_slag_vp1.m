% Calculate Gravitation load on the joints for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:11
% EndTime: 2019-03-09 02:05:12
% DurationCPUTime: 0.50s
% Computational Cost: add. (294->97), mult. (578->142), div. (0->0), fcn. (670->8), ass. (0->40)
t64 = m(6) + m(7);
t58 = rSges(7,1) + pkin(5);
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t55 = sin(qJ(1));
t56 = cos(qJ(1));
t18 = -t55 * t45 - t56 * t46;
t19 = t56 * t45 - t55 * t46;
t63 = g(1) * t18 + g(2) * t19;
t47 = rSges(7,3) + qJ(6);
t29 = sin(qJ(4));
t60 = g(3) * t29;
t59 = t29 * pkin(8);
t57 = rSges(5,3) + pkin(7);
t31 = cos(qJ(4));
t54 = t18 * t31;
t53 = t19 * t31;
t28 = sin(qJ(5));
t52 = t28 * t31;
t51 = t29 * rSges(7,2);
t50 = t29 * rSges(6,3);
t30 = cos(qJ(5));
t49 = t30 * t31;
t48 = t56 * pkin(1) + t55 * qJ(2);
t44 = t56 * pkin(2) + t48;
t43 = m(4) + m(5) + t64;
t42 = -t18 * pkin(3) + t44;
t40 = -t55 * pkin(1) + t56 * qJ(2);
t39 = t31 * rSges(5,1) - t29 * rSges(5,2);
t38 = rSges(6,1) * t30 - rSges(6,2) * t28 + pkin(4);
t2 = t18 * t28 + t19 * t49;
t1 = -t18 * t30 + t19 * t52;
t37 = -pkin(4) * t54 + t19 * pkin(7) - t18 * t59 + t42;
t36 = -t55 * pkin(2) + t40;
t35 = t19 * pkin(3) + t36;
t34 = t47 * t28 + t58 * t30 + pkin(4);
t33 = pkin(4) * t53 + t18 * pkin(7) + t19 * t59 + t35;
t6 = -t18 * t49 + t19 * t28;
t5 = -t18 * t52 - t19 * t30;
t3 = [-m(2) * (g(1) * (-t55 * rSges(2,1) - t56 * rSges(2,2)) + g(2) * (t56 * rSges(2,1) - t55 * rSges(2,2))) - m(3) * (g(1) * (-t55 * rSges(3,1) + t56 * rSges(3,3) + t40) + g(2) * (t56 * rSges(3,1) + t55 * rSges(3,3) + t48)) - m(4) * (g(1) * (t19 * rSges(4,1) - t18 * rSges(4,2) + t36) + g(2) * (-t18 * rSges(4,1) - t19 * rSges(4,2) + t44)) - m(5) * (g(1) * t35 + g(2) * t42 + (g(1) * t39 + g(2) * t57) * t19 + (g(1) * t57 - g(2) * t39) * t18) - m(6) * (g(1) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t19 * t50 + t33) + g(2) * (t6 * rSges(6,1) - t5 * rSges(6,2) - t18 * t50 + t37)) - m(7) * (g(1) * (t47 * t1 + t19 * t51 + t58 * t2 + t33) + g(2) * (-t18 * t51 + t47 * t5 + t58 * t6 + t37)) (-m(3) - t43) * (g(1) * t55 - g(2) * t56) t43 * g(3) (-m(5) * (-g(3) * rSges(5,1) + t63 * rSges(5,2)) - m(6) * (-t63 * rSges(6,3) - g(3) * t38) - m(7) * (-t63 * rSges(7,2) - g(3) * t34)) * t31 + ((-m(5) * rSges(5,2) - m(6) * (-rSges(6,3) - pkin(8)) - m(7) * (-rSges(7,2) - pkin(8))) * g(3) + t63 * (-m(5) * rSges(5,1) - m(6) * t38 - m(7) * t34)) * t29 + t64 * (g(1) * t54 + g(2) * t53) * pkin(8), -m(6) * (g(1) * (-t5 * rSges(6,1) - t6 * rSges(6,2)) + g(2) * (rSges(6,1) * t1 + rSges(6,2) * t2)) - m(7) * (g(1) * (t47 * t6 - t58 * t5) + g(2) * (t1 * t58 - t2 * t47)) + (-m(6) * (rSges(6,1) * t28 + rSges(6,2) * t30) - m(7) * (t58 * t28 - t47 * t30)) * t60, -m(7) * (g(1) * t5 - g(2) * t1 - t28 * t60)];
taug  = t3(:);
