% Calculate potential energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4RPPP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:22
% EndTime: 2018-11-14 13:45:23
% DurationCPUTime: 0.24s
% Computational Cost: add. (192->66), mult. (208->77), div. (0->0), fcn. (177->10), ass. (0->33)
t90 = rSges(5,2) + qJ(3);
t89 = rSges(5,3) + qJ(4);
t88 = rSges(5,1) + pkin(3);
t69 = cos(pkin(4));
t87 = t69 * qJ(2) + pkin(5);
t67 = sin(pkin(4));
t70 = sin(qJ(1));
t86 = t67 * t70;
t71 = cos(qJ(1));
t85 = t67 * t71;
t68 = cos(pkin(6));
t79 = pkin(4) + pkin(6);
t74 = sin(t79) / 0.2e1;
t80 = pkin(4) - pkin(6);
t76 = sin(t80);
t72 = t74 - t76 / 0.2e1;
t50 = t70 * t68 + t71 * t72;
t64 = t70 * pkin(1);
t84 = t50 * pkin(2) + t64;
t83 = t71 * pkin(1) + qJ(2) * t86;
t82 = rSges(4,3) + qJ(3);
t75 = cos(t80) / 0.2e1;
t77 = cos(t79);
t57 = t75 - t77 / 0.2e1;
t81 = t57 * pkin(2) + t87;
t52 = t71 * t68 - t70 * t72;
t78 = t52 * pkin(2) + t83;
t73 = t75 + t77 / 0.2e1;
t66 = sin(pkin(6));
t56 = t74 + t76 / 0.2e1;
t51 = t71 * t66 + t70 * t73;
t49 = t70 * t66 - t71 * t73;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t71 * rSges(2,1) - t70 * rSges(2,2)) + g(2) * (t70 * rSges(2,1) + t71 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t52 * rSges(3,1) - t51 * rSges(3,2) + rSges(3,3) * t86 + t83) + g(2) * (t50 * rSges(3,1) - t49 * rSges(3,2) + t64 + (-rSges(3,3) - qJ(2)) * t85) + g(3) * (t57 * rSges(3,1) + t56 * rSges(3,2) + t69 * rSges(3,3) + t87)) - m(4) * (g(1) * (rSges(4,1) * t86 - t52 * rSges(4,2) + t82 * t51 + t78) + g(2) * (-t50 * rSges(4,2) + (-rSges(4,1) - qJ(2)) * t85 + t82 * t49 + t84) + g(3) * (t69 * rSges(4,1) - t57 * rSges(4,2) - t82 * t56 + t81)) - m(5) * (g(1) * (t90 * t51 + t89 * t52 + t78) + g(2) * (t90 * t49 + t89 * t50 + t84) + g(3) * (-t90 * t56 + t89 * t57 + t88 * t69 + t81) + (g(1) * t88 * t70 + g(2) * (-qJ(2) - t88) * t71) * t67);
U  = t1;
