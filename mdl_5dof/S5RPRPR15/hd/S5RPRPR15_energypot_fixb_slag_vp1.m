% Calculate potential energy for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR15_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR15_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:33
% EndTime: 2019-12-31 18:36:34
% DurationCPUTime: 0.36s
% Computational Cost: add. (106->77), mult. (143->95), div. (0->0), fcn. (127->8), ass. (0->25)
t68 = rSges(6,3) + pkin(7) + qJ(4);
t49 = sin(qJ(1));
t51 = cos(qJ(1));
t67 = -g(1) * t49 + g(2) * t51;
t66 = rSges(5,3) + qJ(4);
t65 = pkin(2) + pkin(5);
t50 = cos(qJ(3));
t62 = rSges(4,2) * t50;
t45 = sin(pkin(8));
t61 = t45 * t51;
t60 = t49 * t45;
t48 = sin(qJ(3));
t59 = t49 * t48;
t58 = t51 * t48;
t56 = t51 * pkin(1) + t49 * qJ(2);
t41 = t49 * pkin(1);
t55 = t49 * pkin(6) + t41;
t53 = t51 * pkin(6) + t56;
t52 = -t51 * qJ(2) + t55;
t46 = cos(pkin(8));
t44 = pkin(8) + qJ(5);
t38 = cos(t44);
t37 = sin(t44);
t36 = pkin(4) * t46 + pkin(3);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + rSges(1,3) * g(3)) - m(2) * (g(1) * (rSges(2,1) * t51 - rSges(2,2) * t49) + g(2) * (rSges(2,1) * t49 + rSges(2,2) * t51) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t51 + rSges(3,3) * t49 + t56) + g(2) * (-t49 * rSges(3,2) + t41 + (-rSges(3,3) - qJ(2)) * t51) + g(3) * (pkin(5) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t59 + t49 * t62 + t53) + g(2) * (t49 * rSges(4,3) + t55) + g(3) * (rSges(4,1) * t50 - rSges(4,2) * t48 + t65) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t48 - qJ(2) - t62)) * t51) - m(5) * (g(1) * (pkin(3) * t59 + (t46 * t59 + t61) * rSges(5,1) + (-t45 * t59 + t46 * t51) * rSges(5,2) + t53) + g(2) * (-pkin(3) * t58 + (-t46 * t58 + t60) * rSges(5,1) + (t45 * t58 + t46 * t49) * rSges(5,2) + t52) + g(3) * (t66 * t48 + t65) + (g(3) * (rSges(5,1) * t46 - rSges(5,2) * t45 + pkin(3)) + t67 * t66) * t50) - m(6) * (g(1) * (t36 * t59 + pkin(4) * t61 + (t37 * t51 + t38 * t59) * rSges(6,1) + (-t37 * t59 + t38 * t51) * rSges(6,2) + t53) + g(2) * (-t36 * t58 + pkin(4) * t60 + (t37 * t49 - t38 * t58) * rSges(6,1) + (t37 * t58 + t38 * t49) * rSges(6,2) + t52) + g(3) * (t68 * t48 + t65) + (g(3) * (rSges(6,1) * t38 - rSges(6,2) * t37 + t36) + t67 * t68) * t50);
U = t1;
