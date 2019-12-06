% Calculate potential energy for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP5_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:40
% EndTime: 2019-12-05 16:47:40
% DurationCPUTime: 0.39s
% Computational Cost: add. (142->85), mult. (180->108), div. (0->0), fcn. (168->8), ass. (0->32)
t82 = rSges(4,3) + pkin(6);
t64 = -pkin(7) - pkin(6);
t81 = rSges(5,3) - t64;
t80 = rSges(6,3) + qJ(5) - t64;
t58 = sin(pkin(8));
t59 = cos(pkin(8));
t79 = g(1) * t59 + g(2) * t58;
t62 = cos(qJ(3));
t49 = t62 * pkin(3) + pkin(2);
t61 = sin(qJ(2));
t75 = rSges(3,2) * t61;
t60 = sin(qJ(3));
t74 = t58 * t60;
t63 = cos(qJ(2));
t73 = t58 * t63;
t72 = t59 * t60;
t71 = t59 * t63;
t70 = t60 * t63;
t69 = t62 * t63;
t66 = t59 * pkin(1) + t58 * pkin(5);
t53 = t58 * pkin(1);
t65 = -t59 * pkin(5) + t53;
t57 = qJ(3) + qJ(4);
t51 = cos(t57);
t50 = sin(t57);
t48 = pkin(3) * t60 + pkin(4) * t50;
t47 = pkin(4) * t51 + t49;
t46 = t50 * t58 + t51 * t71;
t45 = -t50 * t71 + t51 * t58;
t44 = -t50 * t59 + t51 * t73;
t43 = -t50 * t73 - t51 * t59;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t59 - rSges(2,2) * t58) + g(2) * (rSges(2,1) * t58 + rSges(2,2) * t59) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t58 + t66) + g(2) * (rSges(3,1) * t73 - t58 * t75 + t53) + g(3) * (rSges(3,1) * t61 + rSges(3,2) * t63 + qJ(1)) + (g(1) * (rSges(3,1) * t63 - t75) + g(2) * (-rSges(3,3) - pkin(5))) * t59) - m(4) * (g(1) * (pkin(2) * t71 + (t59 * t69 + t74) * rSges(4,1) + (t58 * t62 - t59 * t70) * rSges(4,2) + t66) + g(2) * (pkin(2) * t73 + (t58 * t69 - t72) * rSges(4,1) + (-t58 * t70 - t59 * t62) * rSges(4,2) + t65) + g(3) * (-t82 * t63 + qJ(1)) + (g(3) * (rSges(4,1) * t62 - rSges(4,2) * t60 + pkin(2)) + t79 * t82) * t61) - m(5) * (g(1) * (t46 * rSges(5,1) + t45 * rSges(5,2) + pkin(3) * t74 + t49 * t71 + t66) + g(2) * (t44 * rSges(5,1) + t43 * rSges(5,2) - pkin(3) * t72 + t49 * t73 + t65) + g(3) * (-t81 * t63 + qJ(1)) + (g(3) * (rSges(5,1) * t51 - rSges(5,2) * t50 + t49) + t79 * t81) * t61) - m(6) * (g(1) * (rSges(6,1) * t46 + rSges(6,2) * t45 + t47 * t71 + t48 * t58 + t66) + g(2) * (rSges(6,1) * t44 + rSges(6,2) * t43 + t47 * t73 - t48 * t59 + t65) + g(3) * (-t80 * t63 + qJ(1)) + (g(3) * (rSges(6,1) * t51 - rSges(6,2) * t50 + t47) + t79 * t80) * t61);
U = t1;
