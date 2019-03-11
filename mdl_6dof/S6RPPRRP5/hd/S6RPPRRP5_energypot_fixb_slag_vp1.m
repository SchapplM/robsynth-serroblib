% Calculate potential energy for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:37
% EndTime: 2019-03-09 02:07:37
% DurationCPUTime: 0.32s
% Computational Cost: add. (118->84), mult. (163->97), div. (0->0), fcn. (143->6), ass. (0->31)
t82 = rSges(6,3) + pkin(8);
t81 = rSges(7,3) + qJ(6) + pkin(8);
t57 = sin(qJ(1));
t60 = cos(qJ(1));
t80 = g(1) * t60 + g(2) * t57;
t79 = pkin(2) + pkin(6);
t56 = sin(qJ(4));
t75 = t56 * t57;
t74 = t56 * t60;
t55 = sin(qJ(5));
t73 = t57 * t55;
t58 = cos(qJ(5));
t72 = t57 * t58;
t71 = t60 * t55;
t70 = t60 * t58;
t51 = t57 * pkin(1);
t68 = t57 * qJ(3) + t51;
t67 = t60 * pkin(1) + t57 * qJ(2);
t66 = pkin(3) + t79;
t65 = t60 * pkin(7) + t68;
t64 = t60 * qJ(3) + t67;
t59 = cos(qJ(4));
t63 = rSges(5,1) * t56 + rSges(5,2) * t59;
t62 = -t57 * pkin(7) + t64;
t61 = -t60 * qJ(2) + t65;
t47 = t58 * pkin(5) + pkin(4);
t46 = t56 * t70 - t73;
t45 = -t56 * t71 - t72;
t44 = t56 * t72 + t71;
t43 = -t56 * t73 + t70;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t60 * rSges(2,1) - t57 * rSges(2,2)) + g(2) * (t57 * rSges(2,1) + t60 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-t60 * rSges(3,2) + t57 * rSges(3,3) + t67) + g(2) * (-t57 * rSges(3,2) + t51 + (-rSges(3,3) - qJ(2)) * t60) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (t57 * rSges(4,2) + t60 * rSges(4,3) + t64) + g(2) * (t57 * rSges(4,3) + (-rSges(4,2) - qJ(2)) * t60 + t68) + g(3) * (rSges(4,1) + t79)) - m(5) * (g(1) * t64 + g(2) * t65 + g(3) * (t59 * rSges(5,1) - t56 * rSges(5,2) + t66) + (g(1) * t63 + g(2) * (rSges(5,3) - qJ(2))) * t60 + (g(1) * (-rSges(5,3) - pkin(7)) + g(2) * t63) * t57) - m(6) * (g(1) * (t46 * rSges(6,1) + t45 * rSges(6,2) + pkin(4) * t74 + t62) + g(2) * (t44 * rSges(6,1) + t43 * rSges(6,2) + pkin(4) * t75 + t61) + g(3) * (t82 * t56 + t66) + (g(3) * (rSges(6,1) * t58 - rSges(6,2) * t55 + pkin(4)) - t80 * t82) * t59) - m(7) * (g(1) * (t46 * rSges(7,1) + t45 * rSges(7,2) - pkin(5) * t73 + t47 * t74 + t62) + g(2) * (t44 * rSges(7,1) + t43 * rSges(7,2) + pkin(5) * t71 + t47 * t75 + t61) + g(3) * (t81 * t56 + t66) + (g(3) * (rSges(7,1) * t58 - rSges(7,2) * t55 + t47) - t80 * t81) * t59);
U  = t1;
