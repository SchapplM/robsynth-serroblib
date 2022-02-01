% Calculate Gravitation load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:01:40
% EndTime: 2022-01-20 12:01:41
% DurationCPUTime: 0.30s
% Computational Cost: add. (349->64), mult. (217->81), div. (0->0), fcn. (169->10), ass. (0->41)
t31 = cos(qJ(4));
t27 = qJ(4) + qJ(5);
t20 = sin(t27);
t22 = cos(t27);
t48 = t22 * rSges(6,1) - t20 * rSges(6,2);
t64 = t31 * pkin(4) + t48;
t63 = pkin(8) + rSges(5,3);
t62 = pkin(9) + pkin(8) + rSges(6,3);
t29 = sin(qJ(4));
t61 = t31 * rSges(5,1) - t29 * rSges(5,2);
t60 = -pkin(3) - t61;
t59 = -pkin(3) - t64;
t28 = qJ(1) + qJ(2);
t24 = qJ(3) + t28;
t17 = sin(t24);
t18 = cos(t24);
t58 = g(1) * t18 + g(2) * t17;
t21 = sin(t28);
t57 = pkin(2) * t21;
t30 = sin(qJ(1));
t54 = t30 * pkin(1);
t23 = cos(t28);
t50 = t23 * rSges(3,1) - t21 * rSges(3,2);
t49 = t18 * rSges(4,1) - t17 * rSges(4,2);
t16 = pkin(2) * t23;
t47 = t16 + t49;
t46 = -t21 * rSges(3,1) - t23 * rSges(3,2);
t45 = -t17 * rSges(4,1) - t18 * rSges(4,2);
t44 = -rSges(6,1) * t20 - rSges(6,2) * t22;
t43 = t63 * t17 - t60 * t18;
t42 = t45 - t57;
t41 = t16 + t43;
t40 = t60 * t17 + t63 * t18;
t39 = t62 * t17 - t59 * t18;
t38 = t16 + t39;
t37 = t59 * t17 + t62 * t18;
t36 = t40 - t57;
t35 = t37 - t57;
t32 = cos(qJ(1));
t26 = t32 * pkin(1);
t1 = [-m(2) * (g(1) * (-t30 * rSges(2,1) - t32 * rSges(2,2)) + g(2) * (t32 * rSges(2,1) - t30 * rSges(2,2))) - m(3) * (g(1) * (t46 - t54) + g(2) * (t26 + t50)) - m(4) * (g(1) * (t42 - t54) + g(2) * (t26 + t47)) - m(5) * (g(1) * (t36 - t54) + g(2) * (t26 + t41)) - m(6) * (g(1) * (t35 - t54) + g(2) * (t26 + t38)), -m(3) * (g(1) * t46 + g(2) * t50) - m(4) * (g(1) * t42 + g(2) * t47) - m(5) * (g(1) * t36 + g(2) * t41) - m(6) * (g(1) * t35 + g(2) * t38), -m(4) * (g(1) * t45 + g(2) * t49) - m(5) * (g(1) * t40 + g(2) * t43) - m(6) * (g(1) * t37 + g(2) * t39), (-m(5) * t61 - m(6) * t64) * g(3) + t58 * (-m(5) * (-rSges(5,1) * t29 - rSges(5,2) * t31) - m(6) * (-pkin(4) * t29 + t44)), -m(6) * (g(3) * t48 + t58 * t44)];
taug = t1(:);
