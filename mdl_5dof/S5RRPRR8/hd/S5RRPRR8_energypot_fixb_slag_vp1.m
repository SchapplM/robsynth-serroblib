% Calculate potential energy for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR8_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:49
% EndTime: 2019-12-31 20:16:49
% DurationCPUTime: 0.33s
% Computational Cost: add. (150->70), mult. (130->86), div. (0->0), fcn. (110->10), ass. (0->32)
t77 = rSges(6,3) + pkin(8);
t76 = rSges(3,3) + pkin(6);
t58 = sin(qJ(2));
t74 = t58 * pkin(2) + pkin(5);
t61 = cos(qJ(2));
t48 = t61 * pkin(2) + pkin(1);
t55 = qJ(2) + pkin(9);
t51 = qJ(4) + t55;
t46 = sin(t51);
t73 = rSges(5,2) * t46;
t57 = sin(qJ(5));
t59 = sin(qJ(1));
t72 = t59 * t57;
t60 = cos(qJ(5));
t71 = t59 * t60;
t47 = cos(t51);
t62 = cos(qJ(1));
t70 = t62 * t47;
t69 = t62 * t57;
t68 = t62 * t60;
t56 = -qJ(3) - pkin(6);
t67 = rSges(4,3) - t56;
t50 = cos(t55);
t43 = pkin(3) * t50 + t48;
t54 = -pkin(7) + t56;
t66 = t59 * t43 + t62 * t54;
t49 = sin(t55);
t65 = pkin(3) * t49 + t74;
t64 = rSges(3,1) * t61 - rSges(3,2) * t58 + pkin(1);
t63 = rSges(4,1) * t50 - rSges(4,2) * t49 + t48;
t42 = t62 * t43;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t62 * rSges(2,1) - t59 * rSges(2,2)) + g(2) * (t59 * rSges(2,1) + t62 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(3) * (t58 * rSges(3,1) + t61 * rSges(3,2) + pkin(5)) + (g(1) * t64 - g(2) * t76) * t62 + (g(1) * t76 + g(2) * t64) * t59) - m(4) * (g(3) * (t49 * rSges(4,1) + t50 * rSges(4,2) + t74) + (g(1) * t63 - g(2) * t67) * t62 + (g(1) * t67 + g(2) * t63) * t59) - m(5) * (g(1) * (rSges(5,1) * t70 - t62 * t73 + t42) + g(2) * (-t62 * rSges(5,3) + t66) + g(3) * (t46 * rSges(5,1) + t47 * rSges(5,2) + t65) + (g(1) * (rSges(5,3) - t54) + g(2) * (rSges(5,1) * t47 - t73)) * t59) - m(6) * (g(1) * (pkin(4) * t70 + t42 - t59 * t54 + (t47 * t68 + t72) * rSges(6,1) + (-t47 * t69 + t71) * rSges(6,2)) + g(2) * (t59 * t47 * pkin(4) + (t47 * t71 - t69) * rSges(6,1) + (-t47 * t72 - t68) * rSges(6,2) + t66) + g(3) * (-t77 * t47 + t65) + (g(3) * (rSges(6,1) * t60 - rSges(6,2) * t57 + pkin(4)) + (g(1) * t62 + g(2) * t59) * t77) * t46);
U = t1;
