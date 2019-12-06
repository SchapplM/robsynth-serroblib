% Calculate potential energy for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:32
% EndTime: 2019-12-05 15:46:33
% DurationCPUTime: 0.42s
% Computational Cost: add. (151->79), mult. (154->99), div. (0->0), fcn. (138->10), ass. (0->36)
t88 = rSges(5,3) + pkin(6);
t87 = rSges(6,3) + pkin(7) + pkin(6);
t58 = sin(pkin(8));
t59 = cos(pkin(8));
t86 = g(1) * t59 + g(2) * t58;
t83 = rSges(3,3) + pkin(5);
t56 = qJ(2) + pkin(9);
t51 = sin(t56);
t81 = rSges(4,2) * t51;
t52 = cos(t56);
t80 = t52 * t58;
t79 = t52 * t59;
t57 = qJ(4) + qJ(5);
t53 = sin(t57);
t78 = t53 * t58;
t77 = t53 * t59;
t54 = cos(t57);
t76 = t54 * t58;
t75 = t54 * t59;
t61 = sin(qJ(4));
t74 = t58 * t61;
t63 = cos(qJ(4));
t73 = t58 * t63;
t72 = t59 * t61;
t71 = t59 * t63;
t64 = cos(qJ(2));
t50 = pkin(2) * t64 + pkin(1);
t60 = -qJ(3) - pkin(5);
t69 = t58 * t50 + t59 * t60;
t62 = sin(qJ(2));
t68 = t62 * pkin(2) + qJ(1);
t47 = t59 * t50;
t67 = -t58 * t60 + t47;
t66 = rSges(3,1) * t64 - rSges(3,2) * t62 + pkin(1);
t49 = pkin(4) * t63 + pkin(3);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t59 - rSges(2,2) * t58) + g(2) * (rSges(2,1) * t58 + rSges(2,2) * t59) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t62 + rSges(3,2) * t64 + qJ(1)) + (g(1) * t66 - g(2) * t83) * t59 + (g(1) * t83 + g(2) * t66) * t58) - m(4) * (g(1) * (rSges(4,1) * t79 - t59 * t81 + t47) + g(2) * (-rSges(4,3) * t59 + t69) + g(3) * (rSges(4,1) * t51 + rSges(4,2) * t52 + t68) + (g(1) * (rSges(4,3) - t60) + g(2) * (rSges(4,1) * t52 - t81)) * t58) - m(5) * (g(1) * (pkin(3) * t79 + (t52 * t71 + t74) * rSges(5,1) + (-t52 * t72 + t73) * rSges(5,2) + t67) + g(2) * (pkin(3) * t80 + (t52 * t73 - t72) * rSges(5,1) + (-t52 * t74 - t71) * rSges(5,2) + t69) + g(3) * (-t88 * t52 + t68) + (g(3) * (rSges(5,1) * t63 - rSges(5,2) * t61 + pkin(3)) + t86 * t88) * t51) - m(6) * (g(1) * (t49 * t79 + pkin(4) * t74 + (t52 * t75 + t78) * rSges(6,1) + (-t52 * t77 + t76) * rSges(6,2) + t67) + g(2) * (t49 * t80 - pkin(4) * t72 + (t52 * t76 - t77) * rSges(6,1) + (-t52 * t78 - t75) * rSges(6,2) + t69) + g(3) * (-t87 * t52 + t68) + (g(3) * (rSges(6,1) * t54 - rSges(6,2) * t53 + t49) + t86 * t87) * t51);
U = t1;
