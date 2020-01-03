% Calculate potential energy for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR12_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:28:56
% EndTime: 2019-12-31 20:28:57
% DurationCPUTime: 0.41s
% Computational Cost: add. (120->82), mult. (215->104), div. (0->0), fcn. (216->8), ass. (0->29)
t86 = pkin(8) + rSges(6,3);
t67 = sin(qJ(4));
t72 = cos(qJ(2));
t68 = sin(qJ(2));
t71 = cos(qJ(4));
t84 = t68 * t71;
t88 = t72 * t67 - t84;
t87 = t68 * pkin(2) + pkin(5);
t69 = sin(qJ(1));
t85 = t68 * t69;
t83 = t69 * t72;
t73 = cos(qJ(1));
t81 = t72 * t73;
t80 = t73 * pkin(1) + t69 * pkin(6);
t79 = qJ(3) * t68;
t63 = t69 * pkin(1);
t78 = pkin(2) * t83 + t69 * t79 + t63;
t77 = pkin(2) * t81 + t73 * t79 + t80;
t76 = pkin(3) * t83 + t73 * pkin(7) + t78;
t75 = pkin(3) * t81 + t77;
t50 = t68 * t67 + t72 * t71;
t74 = t68 * pkin(3) - t72 * qJ(3) + t87;
t70 = cos(qJ(5));
t66 = sin(qJ(5));
t49 = t50 * t73;
t48 = t67 * t81 - t73 * t84;
t47 = t50 * t69;
t46 = t88 * t69;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t73 * rSges(2,1) - t69 * rSges(2,2)) + g(2) * (t69 * rSges(2,1) + t73 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t69 * rSges(3,3) + t80) + g(2) * (rSges(3,1) * t83 - rSges(3,2) * t85 + t63) + g(3) * (t68 * rSges(3,1) + t72 * rSges(3,2) + pkin(5)) + (g(1) * (rSges(3,1) * t72 - rSges(3,2) * t68) + g(2) * (-rSges(3,3) - pkin(6))) * t73) - m(4) * (g(1) * (t69 * rSges(4,2) + t77) + g(2) * (rSges(4,1) * t83 + rSges(4,3) * t85 + t78) + g(3) * (t68 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t72 + t87) + (g(1) * (rSges(4,1) * t72 + rSges(4,3) * t68) + g(2) * (-rSges(4,2) - pkin(6))) * t73) - m(5) * (g(1) * (t49 * rSges(5,1) - t48 * rSges(5,2) + (-rSges(5,3) - pkin(7)) * t69 + t75) + g(2) * (t47 * rSges(5,1) - t46 * rSges(5,2) + (rSges(5,3) - pkin(6)) * t73 + t76) + g(3) * (-rSges(5,1) * t88 - t50 * rSges(5,2) + t74)) - m(6) * (g(1) * (t49 * pkin(4) - t69 * pkin(7) + (t49 * t70 - t69 * t66) * rSges(6,1) + (-t49 * t66 - t69 * t70) * rSges(6,2) + t86 * t48 + t75) + g(2) * (t47 * pkin(4) - t73 * pkin(6) + (t47 * t70 + t73 * t66) * rSges(6,1) + (-t47 * t66 + t73 * t70) * rSges(6,2) + t86 * t46 + t76) + (t74 - (t70 * rSges(6,1) - t66 * rSges(6,2) + pkin(4)) * t88 + t86 * t50) * g(3));
U = t1;
