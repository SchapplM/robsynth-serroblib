% Calculate potential energy for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP6_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:06
% EndTime: 2019-12-31 17:18:06
% DurationCPUTime: 0.26s
% Computational Cost: add. (77->60), mult. (125->77), div. (0->0), fcn. (113->6), ass. (0->24)
t61 = rSges(4,3) + pkin(6);
t60 = rSges(5,3) + qJ(4) + pkin(6);
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t59 = g(1) * t46 + g(2) * t43;
t42 = sin(qJ(2));
t55 = rSges(3,2) * t42;
t41 = sin(qJ(3));
t54 = t43 * t41;
t45 = cos(qJ(2));
t53 = t43 * t45;
t52 = t45 * t46;
t51 = t46 * t41;
t44 = cos(qJ(3));
t50 = t46 * t44;
t48 = t46 * pkin(1) + t43 * pkin(5);
t38 = t43 * pkin(1);
t47 = -t46 * pkin(5) + t38;
t36 = t44 * pkin(3) + pkin(2);
t35 = t45 * t50 + t54;
t34 = t43 * t44 - t45 * t51;
t33 = t44 * t53 - t51;
t32 = -t41 * t53 - t50;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t46 * rSges(2,1) - t43 * rSges(2,2)) + g(2) * (t43 * rSges(2,1) + t46 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t43 * rSges(3,3) + t48) + g(2) * (rSges(3,1) * t53 - t43 * t55 + t38) + g(3) * (t42 * rSges(3,1) + t45 * rSges(3,2) + pkin(4)) + (g(1) * (rSges(3,1) * t45 - t55) + g(2) * (-rSges(3,3) - pkin(5))) * t46) - m(4) * (g(1) * (t35 * rSges(4,1) + t34 * rSges(4,2) + pkin(2) * t52 + t48) + g(2) * (t33 * rSges(4,1) + t32 * rSges(4,2) + pkin(2) * t53 + t47) + g(3) * (-t61 * t45 + pkin(4)) + (g(3) * (rSges(4,1) * t44 - rSges(4,2) * t41 + pkin(2)) + t59 * t61) * t42) - m(5) * (g(1) * (t35 * rSges(5,1) + t34 * rSges(5,2) + pkin(3) * t54 + t36 * t52 + t48) + g(2) * (t33 * rSges(5,1) + t32 * rSges(5,2) - pkin(3) * t51 + t36 * t53 + t47) + g(3) * (-t60 * t45 + pkin(4)) + (g(3) * (rSges(5,1) * t44 - rSges(5,2) * t41 + t36) + t59 * t60) * t42);
U = t1;
