% Calculate potential energy for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:48:27
% EndTime: 2019-03-09 01:48:27
% DurationCPUTime: 0.38s
% Computational Cost: add. (128->89), mult. (163->105), div. (0->0), fcn. (143->8), ass. (0->28)
t76 = rSges(7,3) + pkin(8) + qJ(5);
t54 = sin(qJ(1));
t56 = cos(qJ(1));
t75 = g(1) * t56 + g(2) * t54;
t74 = rSges(6,3) + qJ(5);
t73 = pkin(2) + pkin(6);
t50 = sin(pkin(9));
t70 = t54 * t50;
t53 = sin(qJ(4));
t69 = t54 * t53;
t68 = t56 * t50;
t67 = t56 * t53;
t46 = t54 * pkin(1);
t65 = t54 * qJ(3) + t46;
t64 = t56 * pkin(1) + t54 * qJ(2);
t62 = pkin(3) + t73;
t61 = t56 * pkin(7) + t65;
t60 = t56 * qJ(3) + t64;
t55 = cos(qJ(4));
t59 = rSges(5,1) * t53 + rSges(5,2) * t55;
t58 = -t54 * pkin(7) + t60;
t57 = -t56 * qJ(2) + t61;
t51 = cos(pkin(9));
t49 = pkin(9) + qJ(6);
t42 = cos(t49);
t41 = sin(t49);
t40 = pkin(5) * t51 + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t56 - t54 * rSges(2,2)) + g(2) * (t54 * rSges(2,1) + rSges(2,2) * t56) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t56 + t54 * rSges(3,3) + t64) + g(2) * (-t54 * rSges(3,2) + t46 + (-rSges(3,3) - qJ(2)) * t56) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (t54 * rSges(4,2) + rSges(4,3) * t56 + t60) + g(2) * (t54 * rSges(4,3) + (-rSges(4,2) - qJ(2)) * t56 + t65) + g(3) * (rSges(4,1) + t73)) - m(5) * (g(1) * t60 + g(2) * t61 + g(3) * (rSges(5,1) * t55 - rSges(5,2) * t53 + t62) + (g(1) * t59 + g(2) * (rSges(5,3) - qJ(2))) * t56 + (g(1) * (-rSges(5,3) - pkin(7)) + g(2) * t59) * t54) - m(6) * (g(1) * (pkin(4) * t67 + (t51 * t67 - t70) * rSges(6,1) + (-t50 * t67 - t54 * t51) * rSges(6,2) + t58) + g(2) * (pkin(4) * t69 + (t51 * t69 + t68) * rSges(6,1) + (-t50 * t69 + t51 * t56) * rSges(6,2) + t57) + g(3) * (t74 * t53 + t62) + (g(3) * (rSges(6,1) * t51 - rSges(6,2) * t50 + pkin(4)) - t75 * t74) * t55) - m(7) * (g(1) * (t40 * t67 - pkin(5) * t70 + (-t54 * t41 + t42 * t67) * rSges(7,1) + (-t41 * t67 - t54 * t42) * rSges(7,2) + t58) + g(2) * (t40 * t69 + pkin(5) * t68 + (t41 * t56 + t42 * t69) * rSges(7,1) + (-t41 * t69 + t42 * t56) * rSges(7,2) + t57) + g(3) * (t76 * t53 + t62) + (g(3) * (rSges(7,1) * t42 - rSges(7,2) * t41 + t40) - t75 * t76) * t55);
U  = t1;
