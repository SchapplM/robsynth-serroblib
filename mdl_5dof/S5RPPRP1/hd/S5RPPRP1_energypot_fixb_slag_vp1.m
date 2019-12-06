% Calculate potential energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:31
% EndTime: 2019-12-05 17:35:32
% DurationCPUTime: 0.30s
% Computational Cost: add. (147->72), mult. (141->89), div. (0->0), fcn. (125->8), ass. (0->31)
t78 = rSges(5,3) + pkin(6);
t77 = rSges(6,3) + qJ(5) + pkin(6);
t52 = qJ(1) + pkin(7);
t49 = sin(t52);
t50 = cos(t52);
t76 = -g(2) * t49 + g(3) * t50;
t57 = sin(qJ(1));
t73 = t57 * pkin(1);
t53 = sin(pkin(8));
t71 = rSges(4,2) * t53;
t54 = cos(pkin(8));
t70 = t49 * t54;
t56 = sin(qJ(4));
t69 = t49 * t56;
t68 = t50 * t54;
t67 = t50 * t56;
t66 = t54 * t56;
t58 = cos(qJ(4));
t65 = t54 * t58;
t63 = pkin(5) + qJ(2);
t59 = cos(qJ(1));
t51 = t59 * pkin(1);
t62 = t50 * pkin(2) + t49 * qJ(3) + t51;
t61 = t50 * qJ(3) - t73;
t60 = -t49 * pkin(2) + t61;
t48 = t58 * pkin(4) + pkin(3);
t44 = t50 * t65 + t69;
t43 = t49 * t58 - t50 * t66;
t42 = -t49 * t65 + t67;
t41 = t49 * t66 + t50 * t58;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (-t57 * rSges(2,1) - t59 * rSges(2,2)) + g(3) * (t59 * rSges(2,1) - t57 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) + t63) + g(2) * (-t49 * rSges(3,1) - t50 * rSges(3,2) - t73) + g(3) * (t50 * rSges(3,1) - t49 * rSges(3,2) + t51)) - m(4) * (g(1) * (t53 * rSges(4,1) + t54 * rSges(4,2) + t63) + g(2) * (t50 * rSges(4,3) + t61) + g(3) * (rSges(4,1) * t68 - t50 * t71 + t62) + (g(2) * (-rSges(4,1) * t54 - pkin(2) + t71) + g(3) * rSges(4,3)) * t49) - m(5) * (g(1) * (-t78 * t54 + t63) + g(2) * (t42 * rSges(5,1) + t41 * rSges(5,2) - pkin(3) * t70 + t60) + g(3) * (t44 * rSges(5,1) + t43 * rSges(5,2) + pkin(3) * t68 + t62) + (g(1) * (rSges(5,1) * t58 - rSges(5,2) * t56 + pkin(3)) + t76 * t78) * t53) - m(6) * (g(1) * (-t77 * t54 + t63) + g(2) * (t42 * rSges(6,1) + t41 * rSges(6,2) + pkin(4) * t67 - t48 * t70 + t60) + g(3) * (t44 * rSges(6,1) + t43 * rSges(6,2) + pkin(4) * t69 + t48 * t68 + t62) + (g(1) * (rSges(6,1) * t58 - rSges(6,2) * t56 + t48) + t76 * t77) * t53);
U = t1;
