% Calculate Gravitation load on the joints for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:28
% EndTime: 2019-12-31 20:15:28
% DurationCPUTime: 0.25s
% Computational Cost: add. (231->75), mult. (191->93), div. (0->0), fcn. (149->8), ass. (0->40)
t63 = rSges(5,3) + pkin(7);
t62 = rSges(6,3) + pkin(8) + pkin(7);
t31 = qJ(1) + qJ(2);
t28 = cos(t31);
t26 = sin(t31);
t58 = g(1) * t26;
t61 = -g(2) * t28 + t58;
t33 = sin(qJ(1));
t60 = pkin(1) * t33;
t34 = cos(qJ(4));
t59 = pkin(4) * t34;
t30 = qJ(4) + qJ(5);
t27 = cos(t30);
t56 = rSges(6,1) * t27;
t55 = rSges(5,2) * t34;
t25 = sin(t30);
t54 = t25 * t26;
t53 = t25 * t28;
t32 = sin(qJ(4));
t52 = t26 * t32;
t51 = t27 * rSges(6,2);
t50 = t27 * t28;
t49 = t28 * t32;
t48 = t28 * pkin(2) + t26 * qJ(3);
t16 = t28 * qJ(3);
t47 = rSges(5,1) * t49 + t28 * t55 + t16;
t46 = t28 * rSges(3,1) - rSges(3,2) * t26;
t45 = (-pkin(2) - t63) * t58;
t44 = -rSges(3,1) * t26 - rSges(3,2) * t28;
t42 = -rSges(6,1) * t25 - t51;
t41 = t28 * rSges(4,3) + t16 + (rSges(4,2) - pkin(2)) * t26;
t40 = rSges(5,1) * t52 + t26 * t55 + t63 * t28 + t48;
t39 = -rSges(4,2) * t28 + t26 * rSges(4,3) + t48;
t38 = rSges(6,1) * t54 + pkin(4) * t52 + t26 * t51 + t62 * t28 + t48;
t37 = rSges(6,1) * t53 + rSges(6,2) * t50 + pkin(4) * t49 + t16 + (-pkin(2) - t62) * t26;
t35 = cos(qJ(1));
t29 = t35 * pkin(1);
t5 = rSges(6,2) * t53;
t4 = t26 * t56;
t1 = [-m(2) * (g(1) * (-t33 * rSges(2,1) - rSges(2,2) * t35) + g(2) * (rSges(2,1) * t35 - t33 * rSges(2,2))) - m(3) * (g(1) * (t44 - t60) + g(2) * (t29 + t46)) - m(4) * (g(1) * (t41 - t60) + g(2) * (t29 + t39)) - m(5) * (g(1) * (t47 - t60) + g(2) * (t29 + t40) + t45) - m(6) * (g(1) * (t37 - t60) + g(2) * (t29 + t38)), -m(3) * (g(1) * t44 + g(2) * t46) - m(4) * (g(1) * t41 + g(2) * t39) - m(5) * (g(1) * t47 + g(2) * t40 + t45) - m(6) * (g(1) * t37 + g(2) * t38), (-m(4) - m(5) - m(6)) * t61, -m(5) * (g(3) * (-rSges(5,1) * t32 - t55) + t61 * (rSges(5,1) * t34 - rSges(5,2) * t32)) - m(6) * (g(1) * (t4 + (-rSges(6,2) * t25 + t59) * t26) + g(2) * (t5 + (-t56 - t59) * t28) + g(3) * (-t32 * pkin(4) + t42)), -m(6) * (g(1) * (-rSges(6,2) * t54 + t4) + g(2) * (-rSges(6,1) * t50 + t5) + g(3) * t42)];
taug = t1(:);
