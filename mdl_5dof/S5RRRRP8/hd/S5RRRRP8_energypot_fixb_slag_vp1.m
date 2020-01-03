% Calculate potential energy for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP8_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:34
% EndTime: 2019-12-31 21:59:35
% DurationCPUTime: 0.39s
% Computational Cost: add. (142->85), mult. (180->106), div. (0->0), fcn. (168->8), ass. (0->30)
t78 = rSges(4,3) + pkin(7);
t62 = -pkin(8) - pkin(7);
t77 = rSges(5,3) - t62;
t76 = rSges(6,3) + qJ(5) - t62;
t58 = sin(qJ(1));
t61 = cos(qJ(1));
t75 = g(1) * t61 + g(2) * t58;
t59 = cos(qJ(3));
t47 = t59 * pkin(3) + pkin(2);
t57 = sin(qJ(2));
t71 = rSges(3,2) * t57;
t56 = sin(qJ(3));
t70 = t56 * t58;
t69 = t56 * t61;
t60 = cos(qJ(2));
t68 = t58 * t60;
t67 = t60 * t61;
t64 = t61 * pkin(1) + t58 * pkin(6);
t51 = t58 * pkin(1);
t63 = -t61 * pkin(6) + t51;
t55 = qJ(3) + qJ(4);
t49 = cos(t55);
t48 = sin(t55);
t46 = pkin(3) * t56 + pkin(4) * t48;
t45 = pkin(4) * t49 + t47;
t44 = t48 * t58 + t49 * t67;
t43 = -t48 * t67 + t49 * t58;
t42 = -t48 * t61 + t49 * t68;
t41 = -t48 * t68 - t49 * t61;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t61 - rSges(2,2) * t58) + g(2) * (rSges(2,1) * t58 + rSges(2,2) * t61) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t58 + t64) + g(2) * (rSges(3,1) * t68 - t58 * t71 + t51) + g(3) * (rSges(3,1) * t57 + rSges(3,2) * t60 + pkin(5)) + (g(1) * (rSges(3,1) * t60 - t71) + g(2) * (-rSges(3,3) - pkin(6))) * t61) - m(4) * (g(1) * (pkin(2) * t67 + (t59 * t67 + t70) * rSges(4,1) + (-t56 * t67 + t58 * t59) * rSges(4,2) + t64) + g(2) * (pkin(2) * t68 + (t59 * t68 - t69) * rSges(4,1) + (-t56 * t68 - t59 * t61) * rSges(4,2) + t63) + g(3) * (-t78 * t60 + pkin(5)) + (g(3) * (rSges(4,1) * t59 - rSges(4,2) * t56 + pkin(2)) + t75 * t78) * t57) - m(5) * (g(1) * (t44 * rSges(5,1) + t43 * rSges(5,2) + pkin(3) * t70 + t47 * t67 + t64) + g(2) * (t42 * rSges(5,1) + t41 * rSges(5,2) - pkin(3) * t69 + t47 * t68 + t63) + g(3) * (-t77 * t60 + pkin(5)) + (g(3) * (rSges(5,1) * t49 - rSges(5,2) * t48 + t47) + t75 * t77) * t57) - m(6) * (g(1) * (rSges(6,1) * t44 + rSges(6,2) * t43 + t45 * t67 + t46 * t58 + t64) + g(2) * (rSges(6,1) * t42 + rSges(6,2) * t41 + t45 * t68 - t46 * t61 + t63) + g(3) * (-t76 * t60 + pkin(5)) + (g(3) * (rSges(6,1) * t49 - rSges(6,2) * t48 + t45) + t75 * t76) * t57);
U = t1;
