% Calculate Gravitation load on the joints for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:24
% EndTime: 2019-12-31 20:49:25
% DurationCPUTime: 0.31s
% Computational Cost: add. (284->76), mult. (229->90), div. (0->0), fcn. (187->8), ass. (0->39)
t59 = rSges(6,1) + pkin(4);
t64 = rSges(4,3) + pkin(7);
t28 = qJ(3) + pkin(8);
t23 = cos(t28);
t63 = t59 * t23;
t31 = sin(qJ(3));
t33 = cos(qJ(3));
t62 = t33 * rSges(4,1) - t31 * rSges(4,2);
t22 = sin(t28);
t61 = t23 * rSges(5,1) - t22 * rSges(5,2);
t60 = -pkin(2) - t62;
t29 = qJ(1) + qJ(2);
t24 = sin(t29);
t25 = cos(t29);
t58 = g(1) * t25 + g(2) * t24;
t57 = rSges(6,3) + qJ(5);
t56 = pkin(3) * t31;
t32 = sin(qJ(1));
t53 = t32 * pkin(1);
t51 = t22 * t25;
t49 = t23 * t25;
t30 = -qJ(4) - pkin(7);
t48 = t25 * t30;
t45 = t22 * qJ(5);
t44 = t25 * rSges(3,1) - t24 * rSges(3,2);
t43 = t25 * rSges(6,2) - t48;
t26 = t33 * pkin(3);
t21 = t26 + pkin(2);
t8 = t25 * t21;
t42 = t24 * rSges(6,2) + rSges(6,3) * t51 + t25 * t45 + t59 * t49 + t8;
t41 = -t24 * rSges(3,1) - t25 * rSges(3,2);
t40 = t64 * t24 - t60 * t25;
t39 = t60 * t24 + t64 * t25;
t38 = rSges(5,1) * t49 - rSges(5,2) * t51 + t8 + (rSges(5,3) - t30) * t24;
t37 = t25 * rSges(5,3) - t48 + (-t21 - t61) * t24;
t36 = (g(1) * (-t22 * rSges(6,3) - t21 - t45 - t63) - g(2) * t30) * t24;
t34 = cos(qJ(1));
t27 = t34 * pkin(1);
t1 = [-m(2) * (g(1) * (-t32 * rSges(2,1) - t34 * rSges(2,2)) + g(2) * (t34 * rSges(2,1) - t32 * rSges(2,2))) - m(3) * (g(1) * (t41 - t53) + g(2) * (t27 + t44)) - m(4) * (g(1) * (t39 - t53) + g(2) * (t27 + t40)) - m(5) * (g(1) * (t37 - t53) + g(2) * (t27 + t38)) - m(6) * (g(1) * (t43 - t53) + g(2) * (t27 + t42) + t36), -m(3) * (g(1) * t41 + g(2) * t44) - m(4) * (g(1) * t39 + g(2) * t40) - m(5) * (g(1) * t37 + g(2) * t38) - m(6) * (g(1) * t43 + g(2) * t42 + t36), (-m(4) * t62 - m(5) * (t26 + t61) - m(6) * (t57 * t22 + t26 + t63)) * g(3) + t58 * (-m(4) * (-rSges(4,1) * t31 - rSges(4,2) * t33) - m(5) * (-rSges(5,1) * t22 - rSges(5,2) * t23 - t56) - m(6) * (-t59 * t22 + t57 * t23 - t56)), (-m(5) - m(6)) * (g(1) * t24 - g(2) * t25), -m(6) * (-g(3) * t23 + t58 * t22)];
taug = t1(:);
