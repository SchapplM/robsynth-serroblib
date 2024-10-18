% Calculate potential energy for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR13_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR13_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:17
% EndTime: 2024-09-27 17:30:18
% DurationCPUTime: 0.05s
% Computational Cost: add. (195->74), mult. (128->86), div. (0->0), fcn. (106->16), ass. (0->32)
t85 = pkin(9) + pkin(10) + rSges(6,3);
t84 = rSges(5,3) + pkin(9);
t83 = pkin(6) + pkin(7);
t69 = cos(pkin(5));
t70 = sin(qJ(4));
t81 = t69 * t70;
t72 = cos(qJ(4));
t80 = t69 * t72;
t67 = qJ(1) + qJ(2);
t60 = sin(t67);
t71 = sin(qJ(1));
t64 = t71 * pkin(1);
t79 = pkin(2) * t60 + t64;
t62 = cos(t67);
t73 = cos(qJ(1));
t65 = t73 * pkin(1);
t78 = pkin(2) * t62 + t65;
t66 = qJ(4) + qJ(5);
t77 = pkin(8) + t83;
t76 = cos(t66) * rSges(6,1) - sin(t66) * rSges(6,2) + t72 * pkin(4) + pkin(3);
t57 = pkin(5) + t66;
t48 = sin(t57) / 0.2e1;
t58 = pkin(5) - t66;
t49 = cos(t58) / 0.2e1;
t50 = sin(t58);
t51 = cos(t57);
t68 = sin(pkin(5));
t75 = (t48 - t50 / 0.2e1) * rSges(6,1) + (t49 + t51 / 0.2e1) * rSges(6,2) + pkin(4) * t81 - t85 * t68;
t63 = qJ(3) + t67;
t55 = cos(t63);
t54 = sin(t63);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t73 * rSges(2,1) - t71 * rSges(2,2)) + g(2) * (t71 * rSges(2,1) + t73 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t62 * rSges(3,1) - t60 * rSges(3,2) + t65) + g(2) * (t60 * rSges(3,1) + t62 * rSges(3,2) + t64) + g(3) * (rSges(3,3) + t83)) - m(4) * (g(1) * (t55 * rSges(4,1) - t54 * rSges(4,2) + t78) + g(2) * (t54 * rSges(4,1) + t55 * rSges(4,2) + t79) + g(3) * (rSges(4,3) + t77)) - m(5) * (g(1) * (t55 * pkin(3) + (-t54 * t81 + t55 * t72) * rSges(5,1) + (-t54 * t80 - t55 * t70) * rSges(5,2) + t78) + g(2) * (t54 * pkin(3) + (t54 * t72 + t55 * t81) * rSges(5,1) + (-t54 * t70 + t55 * t80) * rSges(5,2) + t79) + g(3) * (t84 * t69 + t77) + (g(3) * (rSges(5,1) * t70 + rSges(5,2) * t72) + (g(1) * t54 - g(2) * t55) * t84) * t68) - m(6) * (g(1) * (-t75 * t54 + t76 * t55 + t78) + g(2) * (t76 * t54 + t75 * t55 + t79) + g(3) * (t68 * t70 * pkin(4) + (t49 - t51 / 0.2e1) * rSges(6,1) + (t48 + t50 / 0.2e1) * rSges(6,2) + t85 * t69 + t77));
U = t1;
