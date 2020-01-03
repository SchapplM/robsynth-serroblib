% Calculate potential energy for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:10:55
% EndTime: 2019-12-31 21:10:56
% DurationCPUTime: 0.25s
% Computational Cost: add. (138->68), mult. (131->82), div. (0->0), fcn. (113->8), ass. (0->25)
t73 = pkin(5) + pkin(6);
t72 = -pkin(8) - rSges(6,3);
t54 = qJ(1) + qJ(2);
t49 = sin(t54);
t56 = sin(qJ(3));
t71 = t49 * t56;
t59 = cos(qJ(3));
t70 = t49 * t59;
t57 = sin(qJ(1));
t52 = t57 * pkin(1);
t69 = t49 * pkin(2) + t52;
t68 = qJ(4) * t56;
t67 = t56 * pkin(3) + t73;
t50 = cos(t54);
t60 = cos(qJ(1));
t53 = t60 * pkin(1);
t66 = t50 * pkin(2) + t49 * pkin(7) + t53;
t65 = pkin(3) * t70 + t49 * t68 + t69;
t64 = t66 + (pkin(3) * t59 + t68) * t50;
t55 = sin(qJ(5));
t58 = cos(qJ(5));
t63 = -t59 * t55 + t56 * t58;
t62 = t56 * t55 + t59 * t58;
t61 = t62 * rSges(6,1) + t63 * rSges(6,2) + t59 * pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t60 * rSges(2,1) - t57 * rSges(2,2)) + g(2) * (t57 * rSges(2,1) + t60 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t50 * rSges(3,1) - t49 * rSges(3,2) + t53) + g(2) * (t49 * rSges(3,1) + t50 * rSges(3,2) + t52) + g(3) * (rSges(3,3) + t73)) - m(4) * (g(1) * (t49 * rSges(4,3) + t66) + g(2) * (rSges(4,1) * t70 - rSges(4,2) * t71 + t69) + g(3) * (t56 * rSges(4,1) + t59 * rSges(4,2) + t73) + (g(1) * (rSges(4,1) * t59 - rSges(4,2) * t56) + g(2) * (-rSges(4,3) - pkin(7))) * t50) - m(5) * (g(1) * (t49 * rSges(5,2) + t64) + g(2) * (rSges(5,1) * t70 + rSges(5,3) * t71 + t65) + g(3) * (t56 * rSges(5,1) + (-rSges(5,3) - qJ(4)) * t59 + t67) + (g(1) * (rSges(5,1) * t59 + rSges(5,3) * t56) + g(2) * (-rSges(5,2) - pkin(7))) * t50) - m(6) * (g(1) * t64 + g(2) * t65 + g(3) * (t63 * rSges(6,1) - t62 * rSges(6,2) + t56 * pkin(4) - t59 * qJ(4) + t67) + (g(1) * t72 + g(2) * t61) * t49 + (g(1) * t61 + g(2) * (-pkin(7) - t72)) * t50);
U = t1;
