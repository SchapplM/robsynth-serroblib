% Calculate potential energy for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR10_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:42:57
% EndTime: 2019-12-31 19:42:58
% DurationCPUTime: 0.46s
% Computational Cost: add. (131->87), mult. (242->113), div. (0->0), fcn. (250->8), ass. (0->29)
t82 = -rSges(6,3) - pkin(7);
t63 = sin(qJ(2));
t81 = t63 * pkin(2) + pkin(5);
t64 = sin(qJ(1));
t80 = t64 * t63;
t66 = cos(qJ(2));
t79 = t64 * t66;
t60 = sin(pkin(8));
t67 = cos(qJ(1));
t78 = t67 * t60;
t61 = cos(pkin(8));
t77 = t67 * t61;
t76 = t67 * t63;
t75 = t67 * pkin(1) + t64 * pkin(6);
t74 = qJ(3) * t63;
t73 = rSges(5,3) + qJ(4);
t72 = t81 + (pkin(3) * t61 + qJ(4) * t60) * t63;
t71 = t75 + (pkin(2) * t66 + t74) * t67;
t48 = t64 * t60 + t66 * t77;
t70 = t48 * pkin(3) + t71;
t58 = t64 * pkin(1);
t69 = pkin(2) * t79 - t67 * pkin(6) + t64 * t74 + t58;
t46 = t61 * t79 - t78;
t68 = t46 * pkin(3) + t69;
t65 = cos(qJ(5));
t62 = sin(qJ(5));
t47 = -t64 * t61 + t66 * t78;
t45 = t60 * t79 + t77;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t67 * rSges(2,1) - t64 * rSges(2,2)) + g(2) * (t64 * rSges(2,1) + t67 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t64 * rSges(3,3) + t75) + g(2) * (rSges(3,1) * t79 - rSges(3,2) * t80 + t58) + g(3) * (t63 * rSges(3,1) + t66 * rSges(3,2) + pkin(5)) + (g(1) * (rSges(3,1) * t66 - rSges(3,2) * t63) + g(2) * (-rSges(3,3) - pkin(6))) * t67) - m(4) * (g(1) * (t48 * rSges(4,1) - t47 * rSges(4,2) + rSges(4,3) * t76 + t71) + g(2) * (t46 * rSges(4,1) - t45 * rSges(4,2) + rSges(4,3) * t80 + t69) + g(3) * ((-rSges(4,3) - qJ(3)) * t66 + (rSges(4,1) * t61 - rSges(4,2) * t60) * t63 + t81)) - m(5) * (g(1) * (t48 * rSges(5,1) + rSges(5,2) * t76 + t47 * t73 + t70) + g(2) * (t46 * rSges(5,1) + rSges(5,2) * t80 + t45 * t73 + t68) + g(3) * ((-rSges(5,2) - qJ(3)) * t66 + (rSges(5,1) * t61 + rSges(5,3) * t60) * t63 + t72)) - m(6) * (g(1) * (t48 * pkin(4) + t47 * qJ(4) + (t47 * t62 + t48 * t65) * rSges(6,1) + (t47 * t65 - t48 * t62) * rSges(6,2) + t70) + g(2) * (t46 * pkin(4) + t45 * qJ(4) + (t45 * t62 + t46 * t65) * rSges(6,1) + (t45 * t65 - t46 * t62) * rSges(6,2) + t68) + (g(1) * t67 + g(2) * t64) * t63 * t82 + (t72 + (-qJ(3) - t82) * t66 + (t61 * pkin(4) + (t60 * t62 + t61 * t65) * rSges(6,1) + (t60 * t65 - t61 * t62) * rSges(6,2)) * t63) * g(3));
U = t1;
