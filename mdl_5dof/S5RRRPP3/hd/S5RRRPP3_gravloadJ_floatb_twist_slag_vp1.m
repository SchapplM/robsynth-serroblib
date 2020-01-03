% Calculate Gravitation load on the joints for
% S5RRRPP3
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:08
% EndTime: 2019-12-31 20:53:09
% DurationCPUTime: 0.34s
% Computational Cost: add. (262->81), mult. (260->93), div. (0->0), fcn. (218->6), ass. (0->39)
t68 = rSges(6,1) + pkin(4);
t33 = cos(qJ(3));
t30 = qJ(1) + qJ(2);
t25 = sin(t30);
t26 = cos(t30);
t43 = g(1) * t26 + g(2) * t25;
t67 = t43 * t33;
t31 = sin(qJ(3));
t53 = t33 * rSges(5,2);
t66 = t31 * rSges(5,3) - t53;
t65 = t33 * rSges(4,1) - t31 * rSges(4,2);
t49 = rSges(6,3) + qJ(5);
t64 = t26 * rSges(5,1) + t25 * t53;
t63 = t68 * t26;
t62 = g(1) * t25;
t32 = sin(qJ(1));
t59 = t32 * pkin(1);
t28 = t33 * pkin(3);
t58 = t26 * t33;
t56 = t31 * rSges(6,2);
t52 = t26 * pkin(2) + t25 * pkin(7);
t27 = t31 * qJ(4);
t51 = t27 + t28;
t48 = -pkin(3) - t49;
t22 = t26 * pkin(7);
t47 = t22 - t59;
t46 = t26 * rSges(3,1) - rSges(3,2) * t25;
t45 = pkin(3) * t58 + t26 * t27 + t52;
t44 = qJ(4) * t67;
t42 = -rSges(3,1) * t25 - rSges(3,2) * t26;
t41 = t25 * rSges(4,3) + t26 * t65 + t52;
t40 = t25 * t68 + t26 * t56 + t49 * t58 + t45;
t39 = t26 * rSges(4,3) + t22 + (-pkin(2) - t65) * t25;
t38 = t25 * rSges(5,1) + t26 * t66 + t45;
t37 = (-t28 - pkin(2) + (-rSges(5,3) - qJ(4)) * t31) * t62;
t36 = (-pkin(2) + (-rSges(6,2) - qJ(4)) * t31 + t48 * t33) * t62;
t34 = cos(qJ(1));
t29 = t34 * pkin(1);
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t32 - rSges(2,2) * t34) + g(2) * (rSges(2,1) * t34 - rSges(2,2) * t32)) - m(3) * (g(1) * (t42 - t59) + g(2) * (t29 + t46)) - m(4) * (g(1) * (t39 - t59) + g(2) * (t29 + t41)) - m(5) * (g(1) * (t47 + t64) + g(2) * (t29 + t38) + t37) - m(6) * (g(1) * (t47 + t63) + g(2) * (t29 + t40) + t36), -m(3) * (g(1) * t42 + g(2) * t46) - m(4) * (g(1) * t39 + g(2) * t41) - m(5) * (g(1) * (t22 + t64) + g(2) * t38 + t37) - m(6) * (g(1) * (t22 + t63) + g(2) * t40 + t36), -m(4) * g(3) * t65 - m(5) * (g(3) * (t51 + t66) + t44) - m(6) * (g(3) * (t33 * t49 + t51 + t56) + t44) + t43 * ((m(4) * rSges(4,2) - m(5) * rSges(5,3) - m(6) * rSges(6,2)) * t33 + (m(4) * rSges(4,1) - m(5) * (rSges(5,2) - pkin(3)) - m(6) * t48) * t31), (-m(5) - m(6)) * (-g(3) * t33 + t31 * t43), -m(6) * (g(3) * t31 + t67)];
taug = t1(:);
