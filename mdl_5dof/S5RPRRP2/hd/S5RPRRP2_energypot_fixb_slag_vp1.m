% Calculate potential energy for
% S5RPRRP2
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
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:02
% EndTime: 2020-01-03 11:45:03
% DurationCPUTime: 0.19s
% Computational Cost: add. (130->54), mult. (85->58), div. (0->0), fcn. (61->8), ass. (0->23)
t52 = rSges(6,1) + pkin(4);
t41 = cos(qJ(1));
t51 = t41 * pkin(1);
t50 = -rSges(5,3) - pkin(7);
t49 = -rSges(6,3) - qJ(5) - pkin(7);
t48 = pkin(5) + qJ(2);
t36 = qJ(1) + pkin(8);
t32 = sin(t36);
t39 = sin(qJ(1));
t35 = t39 * pkin(1);
t47 = pkin(2) * t32 + t35;
t46 = pkin(6) + t48;
t33 = cos(t36);
t45 = -pkin(2) * t33 - t51;
t38 = sin(qJ(4));
t40 = cos(qJ(4));
t44 = rSges(5,1) * t40 - rSges(5,2) * t38 + pkin(3);
t43 = -rSges(6,2) * t38 + t52 * t40 + pkin(3);
t42 = g(2) * t47 + g(3) * t45;
t34 = qJ(3) + t36;
t30 = cos(t34);
t29 = sin(t34);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (t39 * rSges(2,1) + t41 * rSges(2,2)) + g(3) * (-t41 * rSges(2,1) + t39 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) + t48) + g(2) * (t32 * rSges(3,1) + t33 * rSges(3,2) + t35) + g(3) * (-t33 * rSges(3,1) + t32 * rSges(3,2) - t51)) - m(4) * (g(1) * (rSges(4,3) + t46) + g(2) * (t29 * rSges(4,1) + t30 * rSges(4,2) + t47) + g(3) * (-t30 * rSges(4,1) + t29 * rSges(4,2) + t45)) - m(5) * (g(1) * (t38 * rSges(5,1) + t40 * rSges(5,2) + t46) + (g(2) * t50 - g(3) * t44) * t30 + (g(2) * t44 + g(3) * t50) * t29 + t42) - m(6) * (g(1) * (t40 * rSges(6,2) + t52 * t38 + t46) + (g(2) * t49 - g(3) * t43) * t30 + (g(2) * t43 + g(3) * t49) * t29 + t42);
U = t1;
