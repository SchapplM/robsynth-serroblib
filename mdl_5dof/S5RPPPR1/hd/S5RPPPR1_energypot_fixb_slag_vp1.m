% Calculate potential energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:09
% EndTime: 2022-01-20 09:12:09
% DurationCPUTime: 0.39s
% Computational Cost: add. (157->77), mult. (141->98), div. (0->0), fcn. (125->10), ass. (0->31)
t80 = rSges(6,3) + pkin(6) + qJ(4);
t54 = qJ(1) + pkin(7);
t48 = sin(t54);
t50 = cos(t54);
t79 = g(1) * t50 + g(2) * t48;
t78 = rSges(5,3) + qJ(4);
t56 = sin(pkin(8));
t75 = rSges(4,2) * t56;
t55 = sin(pkin(9));
t74 = t48 * t55;
t58 = cos(pkin(8));
t73 = t48 * t58;
t72 = t50 * t55;
t71 = t50 * t58;
t70 = t55 * t58;
t57 = cos(pkin(9));
t69 = t57 * t58;
t46 = t57 * pkin(4) + pkin(3);
t68 = t58 * t46;
t66 = pkin(5) + qJ(2);
t60 = sin(qJ(1));
t51 = t60 * pkin(1);
t65 = t48 * pkin(2) + t51;
t61 = cos(qJ(1));
t52 = t61 * pkin(1);
t63 = t50 * pkin(2) + t48 * qJ(3) + t52;
t62 = -t50 * qJ(3) + t65;
t53 = pkin(9) + qJ(5);
t49 = cos(t53);
t47 = sin(t53);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t61 * rSges(2,1) - t60 * rSges(2,2)) + g(2) * (t60 * rSges(2,1) + t61 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t50 * rSges(3,1) - t48 * rSges(3,2) + t52) + g(2) * (t48 * rSges(3,1) + t50 * rSges(3,2) + t51) + g(3) * (rSges(3,3) + t66)) - m(4) * (g(1) * (t48 * rSges(4,3) + t63) + g(2) * (rSges(4,1) * t73 - t48 * t75 + t65) + g(3) * (t56 * rSges(4,1) + t58 * rSges(4,2) + t66) + (g(1) * (rSges(4,1) * t58 - t75) + g(2) * (-rSges(4,3) - qJ(3))) * t50) - m(5) * (g(1) * (pkin(3) * t71 + (t50 * t69 + t74) * rSges(5,1) + (t48 * t57 - t50 * t70) * rSges(5,2) + t63) + g(2) * (pkin(3) * t73 + (t48 * t69 - t72) * rSges(5,1) + (-t48 * t70 - t50 * t57) * rSges(5,2) + t62) + g(3) * (-t78 * t58 + t66) + (g(3) * (rSges(5,1) * t57 - rSges(5,2) * t55 + pkin(3)) + t79 * t78) * t56) - m(6) * (g(1) * (t50 * t68 + pkin(4) * t74 + (t48 * t47 + t49 * t71) * rSges(6,1) + (-t47 * t71 + t48 * t49) * rSges(6,2) + t63) + g(2) * (t48 * t68 - pkin(4) * t72 + (-t50 * t47 + t49 * t73) * rSges(6,1) + (-t47 * t73 - t50 * t49) * rSges(6,2) + t62) + g(3) * (-t80 * t58 + t66) + (g(3) * (rSges(6,1) * t49 - rSges(6,2) * t47 + t46) + t79 * t80) * t56);
U = t1;
