% Calculate potential energy for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP10_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:19
% EndTime: 2019-12-31 20:09:20
% DurationCPUTime: 0.36s
% Computational Cost: add. (106->79), mult. (172->97), div. (0->0), fcn. (156->6), ass. (0->30)
t84 = rSges(5,3) + pkin(7);
t83 = rSges(6,3) + qJ(5) + pkin(7);
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t82 = g(1) * t63 + g(2) * t60;
t59 = sin(qJ(2));
t78 = t59 * pkin(2) + pkin(5);
t58 = sin(qJ(4));
t77 = t58 * t63;
t76 = t59 * t60;
t75 = t60 * t58;
t61 = cos(qJ(4));
t74 = t60 * t61;
t62 = cos(qJ(2));
t73 = t60 * t62;
t72 = t61 * t63;
t70 = t63 * pkin(1) + t60 * pkin(6);
t69 = qJ(3) * t59;
t68 = t59 * t77;
t67 = t59 * t75;
t55 = t60 * pkin(1);
t66 = pkin(2) * t73 + t60 * t69 + t55;
t65 = t70 + (pkin(2) * t62 + t69) * t63;
t64 = -pkin(6) * t63 + t66;
t52 = pkin(4) * t61 + pkin(3);
t47 = t67 - t72;
t46 = t59 * t74 + t77;
t45 = t68 + t74;
t44 = t59 * t72 - t75;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t63 - t60 * rSges(2,2)) + g(2) * (t60 * rSges(2,1) + rSges(2,2) * t63) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t60 * rSges(3,3) + t70) + g(2) * (rSges(3,1) * t73 - rSges(3,2) * t76 + t55) + g(3) * (rSges(3,1) * t59 + rSges(3,2) * t62 + pkin(5)) + (g(1) * (rSges(3,1) * t62 - rSges(3,2) * t59) + g(2) * (-rSges(3,3) - pkin(6))) * t63) - m(4) * (g(1) * (t60 * rSges(4,1) + t65) + g(2) * (-rSges(4,2) * t73 + rSges(4,3) * t76 + t66) + g(3) * (-rSges(4,2) * t59 + (-rSges(4,3) - qJ(3)) * t62 + t78) + (g(1) * (-rSges(4,2) * t62 + rSges(4,3) * t59) + g(2) * (-rSges(4,1) - pkin(6))) * t63) - m(5) * (g(1) * (t45 * rSges(5,1) + t44 * rSges(5,2) + t60 * pkin(3) + t65) + g(2) * (t47 * rSges(5,1) + t46 * rSges(5,2) - pkin(3) * t63 + t64) + g(3) * (t84 * t59 + t78) + (g(3) * (-rSges(5,1) * t58 - rSges(5,2) * t61 - qJ(3)) + t82 * t84) * t62) - m(6) * (g(1) * (t45 * rSges(6,1) + t44 * rSges(6,2) + pkin(4) * t68 + t60 * t52 + t65) + g(2) * (t47 * rSges(6,1) + t46 * rSges(6,2) + pkin(4) * t67 - t52 * t63 + t64) + g(3) * (t83 * t59 + t78) + (g(3) * (-rSges(6,2) * t61 - qJ(3) + (-rSges(6,1) - pkin(4)) * t58) + t82 * t83) * t62);
U = t1;
