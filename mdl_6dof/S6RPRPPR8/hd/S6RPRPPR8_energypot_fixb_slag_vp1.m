% Calculate potential energy for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR8_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:41
% EndTime: 2019-03-09 02:58:41
% DurationCPUTime: 0.34s
% Computational Cost: add. (121->86), mult. (179->97), div. (0->0), fcn. (155->6), ass. (0->29)
t78 = pkin(2) + pkin(6);
t77 = -pkin(3) - pkin(4);
t54 = sin(qJ(1));
t76 = g(1) * t54;
t57 = cos(qJ(1));
t75 = g(2) * t57;
t74 = rSges(7,3) + pkin(8);
t56 = cos(qJ(3));
t73 = rSges(6,1) * t56;
t72 = rSges(4,2) * t56;
t53 = sin(qJ(3));
t71 = t54 * t53;
t70 = t57 * pkin(1) + t54 * qJ(2);
t47 = t54 * pkin(1);
t69 = t54 * pkin(7) + t47;
t68 = t56 * qJ(4);
t67 = -rSges(6,3) - qJ(5);
t66 = t57 * pkin(7) + t70;
t65 = t56 * pkin(3) + t53 * qJ(4) + t78;
t64 = pkin(3) * t71 + t66;
t63 = g(2) * (t57 * t68 + t69);
t62 = t56 * pkin(4) + t65;
t61 = rSges(5,1) * t53 - rSges(5,3) * t56;
t52 = sin(qJ(6));
t55 = cos(qJ(6));
t60 = rSges(7,1) * t55 - rSges(7,2) * t52 + pkin(5);
t59 = -t52 * rSges(7,1) - t55 * rSges(7,2) - qJ(5);
t58 = g(1) * (pkin(4) * t71 + t64) + t63;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t57 - t54 * rSges(2,2)) + g(2) * (t54 * rSges(2,1) + rSges(2,2) * t57) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t57 + t54 * rSges(3,3) + t70) + g(2) * (-t54 * rSges(3,2) + t47 + (-rSges(3,3) - qJ(2)) * t57) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t71 + t54 * t72 + t66) + g(2) * (t54 * rSges(4,3) + t69) + g(3) * (rSges(4,1) * t56 - rSges(4,2) * t53 + t78) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t53 - qJ(2) - t72)) * t57) - m(5) * (g(1) * t64 + t63 + g(3) * (rSges(5,1) * t56 + rSges(5,3) * t53 + t65) + (g(1) * (t61 - t68) + g(2) * rSges(5,2)) * t54 + (g(1) * rSges(5,2) + g(2) * (-t53 * pkin(3) - qJ(2) - t61)) * t57) - m(6) * (g(3) * (rSges(6,1) * t53 - rSges(6,2) * t56 + t62) + (g(1) * (-rSges(6,2) * t53 - t68 - t73) + g(2) * t67) * t54 + (g(1) * t67 + (-qJ(2) + t73 + (rSges(6,2) + t77) * t53) * g(2)) * t57 + t58) - m(7) * (g(3) * t62 + (g(1) * t59 - g(2) * qJ(2)) * t57 + g(2) * t59 * t54 + (g(3) * t74 + t60 * t75 + (-qJ(4) - t60) * t76) * t56 + (g(3) * t60 + t74 * t76 + (-t74 + t77) * t75) * t53 + t58);
U  = t1;
