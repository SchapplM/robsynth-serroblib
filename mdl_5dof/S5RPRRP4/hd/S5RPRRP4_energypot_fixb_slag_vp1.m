% Calculate potential energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:05:25
% EndTime: 2019-12-05 18:05:25
% DurationCPUTime: 0.38s
% Computational Cost: add. (142->85), mult. (180->106), div. (0->0), fcn. (168->8), ass. (0->32)
t82 = rSges(4,3) + pkin(6);
t64 = -pkin(7) - pkin(6);
t81 = rSges(5,3) - t64;
t80 = rSges(6,3) + qJ(5) - t64;
t61 = sin(qJ(1));
t63 = cos(qJ(1));
t79 = -g(2) * t61 + g(3) * t63;
t62 = cos(qJ(3));
t49 = t62 * pkin(3) + pkin(2);
t58 = sin(pkin(8));
t75 = rSges(3,2) * t58;
t59 = cos(pkin(8));
t74 = t59 * t61;
t73 = t59 * t63;
t60 = sin(qJ(3));
t72 = t60 * t61;
t71 = t60 * t63;
t70 = t61 * t62;
t69 = t62 * t63;
t66 = t63 * pkin(1) + t61 * qJ(2);
t53 = t63 * qJ(2);
t65 = -t61 * pkin(1) + t53;
t57 = qJ(3) + qJ(4);
t51 = cos(t57);
t50 = sin(t57);
t48 = pkin(3) * t60 + pkin(4) * t50;
t47 = pkin(4) * t51 + t49;
t46 = t50 * t61 + t51 * t73;
t45 = -t50 * t73 + t51 * t61;
t44 = t50 * t63 - t51 * t74;
t43 = t50 * t74 + t51 * t63;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (-rSges(2,1) * t61 - rSges(2,2) * t63) + g(3) * (rSges(2,1) * t63 - rSges(2,2) * t61)) - m(3) * (g(1) * (rSges(3,1) * t58 + rSges(3,2) * t59 + pkin(5)) + g(2) * (rSges(3,3) * t63 + t53) + g(3) * (rSges(3,1) * t73 - t63 * t75 + t66) + (g(2) * (-rSges(3,1) * t59 - pkin(1) + t75) + g(3) * rSges(3,3)) * t61) - m(4) * (g(1) * (-t82 * t59 + pkin(5)) + g(2) * (-pkin(2) * t74 + (-t59 * t70 + t71) * rSges(4,1) + (t59 * t72 + t69) * rSges(4,2) + t65) + g(3) * (pkin(2) * t73 + (t59 * t69 + t72) * rSges(4,1) + (-t59 * t71 + t70) * rSges(4,2) + t66) + (g(1) * (rSges(4,1) * t62 - rSges(4,2) * t60 + pkin(2)) + t79 * t82) * t58) - m(5) * (g(1) * (-t81 * t59 + pkin(5)) + g(2) * (t44 * rSges(5,1) + t43 * rSges(5,2) + pkin(3) * t71 - t49 * t74 + t65) + g(3) * (t46 * rSges(5,1) + t45 * rSges(5,2) + pkin(3) * t72 + t49 * t73 + t66) + (g(1) * (rSges(5,1) * t51 - rSges(5,2) * t50 + t49) + t79 * t81) * t58) - m(6) * (g(1) * (-t80 * t59 + pkin(5)) + g(2) * (rSges(6,1) * t44 + rSges(6,2) * t43 - t47 * t74 + t48 * t63 + t65) + g(3) * (rSges(6,1) * t46 + rSges(6,2) * t45 + t47 * t73 + t48 * t61 + t66) + (g(1) * (rSges(6,1) * t51 - rSges(6,2) * t50 + t47) + t79 * t80) * t58);
U = t1;
