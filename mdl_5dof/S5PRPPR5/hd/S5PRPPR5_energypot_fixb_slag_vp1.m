% Calculate potential energy for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPPR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR5_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:37:58
% EndTime: 2019-12-31 17:37:58
% DurationCPUTime: 0.34s
% Computational Cost: add. (120->83), mult. (215->105), div. (0->0), fcn. (216->8), ass. (0->29)
t86 = pkin(6) + rSges(6,3);
t67 = sin(pkin(7));
t71 = sin(qJ(2));
t85 = t67 * t71;
t73 = cos(qJ(2));
t84 = t67 * t73;
t69 = cos(pkin(7));
t83 = t69 * t73;
t68 = cos(pkin(8));
t82 = t71 * t68;
t81 = t69 * pkin(1) + t67 * pkin(5);
t80 = qJ(3) * t71;
t79 = t71 * pkin(2) + qJ(1);
t62 = t67 * pkin(1);
t78 = pkin(2) * t84 + t67 * t80 + t62;
t77 = pkin(2) * t83 + t69 * t80 + t81;
t76 = pkin(3) * t84 + t69 * qJ(4) + t78;
t75 = pkin(3) * t83 + t77;
t66 = sin(pkin(8));
t50 = t71 * t66 + t68 * t73;
t74 = t71 * pkin(3) - t73 * qJ(3) + t79;
t72 = cos(qJ(5));
t70 = sin(qJ(5));
t51 = -t66 * t73 + t82;
t49 = t50 * t69;
t48 = t66 * t83 - t69 * t82;
t47 = t50 * t67;
t46 = t66 * t84 - t67 * t82;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t69 - rSges(2,2) * t67) + g(2) * (rSges(2,1) * t67 + rSges(2,2) * t69) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t67 * rSges(3,3) + t81) + g(2) * (rSges(3,1) * t84 - rSges(3,2) * t85 + t62) + g(3) * (t71 * rSges(3,1) + rSges(3,2) * t73 + qJ(1)) + (g(1) * (rSges(3,1) * t73 - rSges(3,2) * t71) + g(2) * (-rSges(3,3) - pkin(5))) * t69) - m(4) * (g(1) * (t67 * rSges(4,2) + t77) + g(2) * (rSges(4,1) * t84 + rSges(4,3) * t85 + t78) + g(3) * (t71 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t73 + t79) + (g(1) * (rSges(4,1) * t73 + rSges(4,3) * t71) + g(2) * (-rSges(4,2) - pkin(5))) * t69) - m(5) * (g(1) * (rSges(5,1) * t49 - rSges(5,2) * t48 + (-rSges(5,3) - qJ(4)) * t67 + t75) + g(2) * (rSges(5,1) * t47 - rSges(5,2) * t46 + (rSges(5,3) - pkin(5)) * t69 + t76) + g(3) * (t51 * rSges(5,1) - t50 * rSges(5,2) + t74)) - m(6) * (g(1) * (t49 * pkin(4) - t67 * qJ(4) + (t49 * t72 - t67 * t70) * rSges(6,1) + (-t49 * t70 - t67 * t72) * rSges(6,2) + t86 * t48 + t75) + g(2) * (t47 * pkin(4) - t69 * pkin(5) + (t47 * t72 + t69 * t70) * rSges(6,1) + (-t47 * t70 + t69 * t72) * rSges(6,2) + t86 * t46 + t76) + (t74 + (t72 * rSges(6,1) - t70 * rSges(6,2) + pkin(4)) * t51 + t86 * t50) * g(3));
U = t1;
