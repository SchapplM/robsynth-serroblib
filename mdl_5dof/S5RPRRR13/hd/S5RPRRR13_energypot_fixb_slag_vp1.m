% Calculate potential energy for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR13_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR13_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:28
% EndTime: 2019-12-31 19:14:28
% DurationCPUTime: 0.36s
% Computational Cost: add. (106->77), mult. (143->95), div. (0->0), fcn. (127->8), ass. (0->26)
t71 = rSges(5,3) + pkin(7);
t70 = rSges(6,3) + pkin(8) + pkin(7);
t49 = sin(qJ(1));
t52 = cos(qJ(1));
t69 = -g(1) * t49 + g(2) * t52;
t68 = pkin(2) + pkin(5);
t51 = cos(qJ(3));
t64 = rSges(4,2) * t51;
t47 = sin(qJ(4));
t63 = t47 * t49;
t62 = t47 * t52;
t48 = sin(qJ(3));
t61 = t49 * t48;
t50 = cos(qJ(4));
t60 = t49 * t50;
t59 = t52 * t48;
t57 = t52 * pkin(1) + t49 * qJ(2);
t43 = t49 * pkin(1);
t56 = t49 * pkin(6) + t43;
t55 = t52 * pkin(6) + t57;
t54 = -t52 * qJ(2) + t56;
t46 = qJ(4) + qJ(5);
t40 = cos(t46);
t39 = sin(t46);
t38 = pkin(4) * t50 + pkin(3);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t52 - rSges(2,2) * t49) + g(2) * (rSges(2,1) * t49 + rSges(2,2) * t52) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t52 + rSges(3,3) * t49 + t57) + g(2) * (-rSges(3,2) * t49 + t43 + (-rSges(3,3) - qJ(2)) * t52) + g(3) * (pkin(5) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t61 + t49 * t64 + t55) + g(2) * (rSges(4,3) * t49 + t56) + g(3) * (rSges(4,1) * t51 - rSges(4,2) * t48 + t68) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t48 - qJ(2) - t64)) * t52) - m(5) * (g(1) * (pkin(3) * t61 + (t48 * t60 + t62) * rSges(5,1) + (-t47 * t61 + t50 * t52) * rSges(5,2) + t55) + g(2) * (-pkin(3) * t59 + (-t50 * t59 + t63) * rSges(5,1) + (t47 * t59 + t60) * rSges(5,2) + t54) + g(3) * (t71 * t48 + t68) + (g(3) * (rSges(5,1) * t50 - rSges(5,2) * t47 + pkin(3)) + t69 * t71) * t51) - m(6) * (g(1) * (t38 * t61 + pkin(4) * t62 + (t39 * t52 + t40 * t61) * rSges(6,1) + (-t39 * t61 + t40 * t52) * rSges(6,2) + t55) + g(2) * (-t38 * t59 + pkin(4) * t63 + (t39 * t49 - t40 * t59) * rSges(6,1) + (t39 * t59 + t40 * t49) * rSges(6,2) + t54) + g(3) * (t70 * t48 + t68) + (g(3) * (rSges(6,1) * t40 - rSges(6,2) * t39 + t38) + t69 * t70) * t51);
U = t1;
