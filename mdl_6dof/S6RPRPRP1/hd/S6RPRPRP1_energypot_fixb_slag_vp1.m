% Calculate potential energy for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:01:15
% EndTime: 2019-03-09 03:01:15
% DurationCPUTime: 0.36s
% Computational Cost: add. (224->88), mult. (172->103), div. (0->0), fcn. (152->10), ass. (0->39)
t97 = rSges(6,3) + pkin(8);
t96 = rSges(7,3) + qJ(6) + pkin(8);
t68 = qJ(1) + pkin(9);
t61 = sin(t68);
t63 = cos(t68);
t95 = g(1) * t63 + g(2) * t61;
t92 = rSges(4,3) + pkin(7);
t67 = qJ(3) + pkin(10);
t60 = sin(t67);
t90 = rSges(5,2) * t60;
t62 = cos(t67);
t89 = t61 * t62;
t71 = sin(qJ(5));
t88 = t61 * t71;
t74 = cos(qJ(5));
t87 = t61 * t74;
t86 = t62 * t63;
t85 = t63 * t71;
t84 = t63 * t74;
t82 = pkin(6) + qJ(2);
t75 = cos(qJ(3));
t59 = pkin(3) * t75 + pkin(2);
t76 = cos(qJ(1));
t66 = t76 * pkin(1);
t81 = t63 * t59 + t66;
t72 = sin(qJ(3));
t80 = t72 * pkin(3) + t82;
t73 = sin(qJ(1));
t65 = t73 * pkin(1);
t70 = -qJ(4) - pkin(7);
t79 = t61 * t59 + t63 * t70 + t65;
t78 = -t61 * t70 + t81;
t77 = rSges(4,1) * t75 - rSges(4,2) * t72 + pkin(2);
t58 = pkin(5) * t74 + pkin(4);
t54 = t62 * t84 + t88;
t53 = -t62 * t85 + t87;
t52 = t62 * t87 - t85;
t51 = -t62 * t88 - t84;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t76 * rSges(2,1) - t73 * rSges(2,2)) + g(2) * (t73 * rSges(2,1) + t76 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t63 * rSges(3,1) - t61 * rSges(3,2) + t66) + g(2) * (t61 * rSges(3,1) + t63 * rSges(3,2) + t65) + g(3) * (rSges(3,3) + t82)) - m(4) * (g(1) * t66 + g(2) * t65 + g(3) * (t72 * rSges(4,1) + t75 * rSges(4,2) + t82) + (g(1) * t77 - g(2) * t92) * t63 + (g(1) * t92 + g(2) * t77) * t61) - m(5) * (g(1) * (rSges(5,1) * t86 - t63 * t90 + t81) + g(2) * (-t63 * rSges(5,3) + t79) + g(3) * (t60 * rSges(5,1) + t62 * rSges(5,2) + t80) + (g(1) * (rSges(5,3) - t70) + g(2) * (rSges(5,1) * t62 - t90)) * t61) - m(6) * (g(1) * (t54 * rSges(6,1) + t53 * rSges(6,2) + pkin(4) * t86 + t78) + g(2) * (t52 * rSges(6,1) + t51 * rSges(6,2) + pkin(4) * t89 + t79) + g(3) * (-t97 * t62 + t80) + (g(3) * (rSges(6,1) * t74 - rSges(6,2) * t71 + pkin(4)) + t95 * t97) * t60) - m(7) * (g(1) * (t54 * rSges(7,1) + t53 * rSges(7,2) + pkin(5) * t88 + t58 * t86 + t78) + g(2) * (t52 * rSges(7,1) + t51 * rSges(7,2) - pkin(5) * t85 + t58 * t89 + t79) + g(3) * (-t96 * t62 + t80) + (g(3) * (rSges(7,1) * t74 - rSges(7,2) * t71 + t58) + t95 * t96) * t60);
U  = t1;
