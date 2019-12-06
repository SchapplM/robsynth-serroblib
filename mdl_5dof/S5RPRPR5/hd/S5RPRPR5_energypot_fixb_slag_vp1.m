% Calculate potential energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:55:46
% EndTime: 2019-12-05 17:55:46
% DurationCPUTime: 0.44s
% Computational Cost: add. (152->90), mult. (180->114), div. (0->0), fcn. (168->10), ass. (0->39)
t88 = rSges(4,3) + pkin(6);
t58 = -qJ(4) - pkin(6);
t87 = rSges(5,3) - t58;
t86 = rSges(6,3) + pkin(7) - t58;
t60 = sin(qJ(1));
t62 = cos(qJ(1));
t85 = -g(2) * t60 + g(3) * t62;
t61 = cos(qJ(3));
t46 = t61 * pkin(3) + pkin(2);
t56 = sin(pkin(8));
t81 = rSges(3,2) * t56;
t57 = cos(pkin(8));
t80 = t57 * t60;
t79 = t57 * t62;
t55 = qJ(3) + pkin(9);
t49 = qJ(5) + t55;
t44 = sin(t49);
t78 = t60 * t44;
t45 = cos(t49);
t77 = t60 * t45;
t47 = sin(t55);
t76 = t60 * t47;
t48 = cos(t55);
t75 = t60 * t48;
t59 = sin(qJ(3));
t74 = t60 * t59;
t73 = t60 * t61;
t72 = t62 * t44;
t71 = t62 * t45;
t70 = t62 * t47;
t69 = t62 * t48;
t68 = t62 * t59;
t67 = t62 * t61;
t64 = t62 * pkin(1) + t60 * qJ(2);
t51 = t62 * qJ(2);
t63 = -t60 * pkin(1) + t51;
t43 = t59 * pkin(3) + pkin(4) * t47;
t42 = pkin(4) * t48 + t46;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (-t60 * rSges(2,1) - t62 * rSges(2,2)) + g(3) * (t62 * rSges(2,1) - t60 * rSges(2,2))) - m(3) * (g(1) * (t56 * rSges(3,1) + t57 * rSges(3,2) + pkin(5)) + g(2) * (t62 * rSges(3,3) + t51) + g(3) * (rSges(3,1) * t79 - t62 * t81 + t64) + (g(2) * (-rSges(3,1) * t57 - pkin(1) + t81) + g(3) * rSges(3,3)) * t60) - m(4) * (g(1) * (-t88 * t57 + pkin(5)) + g(2) * (-pkin(2) * t80 + (-t57 * t73 + t68) * rSges(4,1) + (t57 * t74 + t67) * rSges(4,2) + t63) + g(3) * (pkin(2) * t79 + (t57 * t67 + t74) * rSges(4,1) + (-t57 * t68 + t73) * rSges(4,2) + t64) + (g(1) * (rSges(4,1) * t61 - rSges(4,2) * t59 + pkin(2)) + t85 * t88) * t56) - m(5) * (g(1) * (-t87 * t57 + pkin(5)) + g(2) * (-t46 * t80 + pkin(3) * t68 + (-t57 * t75 + t70) * rSges(5,1) + (t57 * t76 + t69) * rSges(5,2) + t63) + g(3) * (t46 * t79 + pkin(3) * t74 + (t57 * t69 + t76) * rSges(5,1) + (-t57 * t70 + t75) * rSges(5,2) + t64) + (g(1) * (rSges(5,1) * t48 - rSges(5,2) * t47 + t46) + t85 * t87) * t56) - m(6) * (g(1) * (-t86 * t57 + pkin(5)) + g(2) * (-t42 * t80 + t62 * t43 + (-t57 * t77 + t72) * rSges(6,1) + (t57 * t78 + t71) * rSges(6,2) + t63) + g(3) * (t42 * t79 + t60 * t43 + (t57 * t71 + t78) * rSges(6,1) + (-t57 * t72 + t77) * rSges(6,2) + t64) + (g(1) * (rSges(6,1) * t45 - rSges(6,2) * t44 + t42) + t85 * t86) * t56);
U = t1;
