% Calculate Gravitation load on the joints for
% S5RRRPP2
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:29
% EndTime: 2019-12-31 20:51:30
% DurationCPUTime: 0.32s
% Computational Cost: add. (261->79), mult. (257->92), div. (0->0), fcn. (215->6), ass. (0->36)
t54 = rSges(6,1) + pkin(4);
t27 = sin(qJ(3));
t29 = cos(qJ(3));
t63 = t29 * rSges(5,1) + t27 * rSges(5,3);
t62 = t29 * rSges(4,1) - t27 * rSges(4,2);
t26 = qJ(1) + qJ(2);
t21 = sin(t26);
t22 = cos(t26);
t61 = g(1) * t22 + g(2) * t21;
t60 = -rSges(6,3) - qJ(5);
t59 = g(1) * t21;
t28 = sin(qJ(1));
t56 = t28 * pkin(1);
t55 = -rSges(5,1) - pkin(3);
t53 = t22 * t29;
t51 = t27 * rSges(6,2);
t19 = t22 * pkin(7);
t47 = t22 * rSges(5,2) + t19;
t46 = t22 * pkin(2) + t21 * pkin(7);
t23 = t27 * qJ(4);
t45 = t29 * pkin(3) + t23;
t43 = -pkin(3) - t54;
t42 = t22 * rSges(3,1) - rSges(3,2) * t21;
t41 = pkin(3) * t53 + t22 * t23 + t46;
t40 = t61 * qJ(4) * t29;
t39 = -rSges(3,1) * t21 - rSges(3,2) * t22;
t38 = t22 * t51 + t53 * t54 + t41;
t37 = t21 * rSges(5,2) + t22 * t63 + t41;
t36 = t22 * t60 + t19;
t35 = t21 * rSges(4,3) + t22 * t62 + t46;
t34 = t22 * rSges(4,3) + t19 + (-pkin(2) - t62) * t21;
t33 = (-pkin(2) + t55 * t29 + (-rSges(5,3) - qJ(4)) * t27) * t59;
t32 = (g(2) * t60 + (t29 * t43 - pkin(2) - t23 - t51) * g(1)) * t21;
t30 = cos(qJ(1));
t25 = t30 * pkin(1);
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t28 - rSges(2,2) * t30) + g(2) * (rSges(2,1) * t30 - rSges(2,2) * t28)) - m(3) * (g(1) * (t39 - t56) + g(2) * (t25 + t42)) - m(4) * (g(1) * (t34 - t56) + g(2) * (t25 + t35)) - m(5) * (g(1) * (t47 - t56) + g(2) * (t25 + t37) + t33) - m(6) * (g(1) * (t36 - t56) + g(2) * (t25 + t38) + t32), -m(3) * (g(1) * t39 + g(2) * t42) - m(4) * (g(1) * t34 + g(2) * t35) - m(5) * (g(1) * t47 + g(2) * t37 + t33) - m(6) * (g(1) * t36 + g(2) * t38 + t32), -m(4) * g(3) * t62 - m(5) * (g(3) * (t45 + t63) + t40) - m(6) * (g(3) * (t29 * t54 + t45 + t51) + t40) + t61 * ((m(4) * rSges(4,2) - m(5) * rSges(5,3) - m(6) * rSges(6,2)) * t29 + (m(4) * rSges(4,1) - m(5) * t55 - m(6) * t43) * t27), (-m(5) - m(6)) * (-g(3) * t29 + t27 * t61), -m(6) * (g(2) * t22 - t59)];
taug = t1(:);
