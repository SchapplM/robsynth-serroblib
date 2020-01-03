% Calculate potential energy for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR9_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:36
% EndTime: 2019-12-31 19:40:36
% DurationCPUTime: 0.34s
% Computational Cost: add. (99->79), mult. (159->99), div. (0->0), fcn. (139->6), ass. (0->24)
t72 = rSges(6,3) + pkin(7);
t51 = sin(qJ(2));
t70 = t51 * pkin(2) + pkin(5);
t52 = sin(qJ(1));
t69 = t52 * t51;
t53 = cos(qJ(5));
t68 = t52 * t53;
t54 = cos(qJ(2));
t67 = t52 * t54;
t55 = cos(qJ(1));
t66 = t55 * t51;
t65 = t55 * t53;
t64 = t55 * t54;
t63 = t55 * pkin(1) + t52 * pkin(6);
t62 = qJ(3) * t51;
t61 = t51 * pkin(3) + t70;
t48 = t52 * pkin(1);
t60 = pkin(2) * t67 + t52 * t62 + t48;
t59 = pkin(2) * t64 + t55 * t62 + t63;
t58 = pkin(3) * t67 + t55 * qJ(4) + t60;
t57 = pkin(3) * t64 + t59;
t56 = rSges(5,1) * t51 - rSges(5,2) * t54;
t50 = sin(qJ(5));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t55 - rSges(2,2) * t52) + g(2) * (rSges(2,1) * t52 + rSges(2,2) * t55) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t52 * rSges(3,3) + t63) + g(2) * (rSges(3,1) * t67 - rSges(3,2) * t69 + t48) + g(3) * (rSges(3,1) * t51 + rSges(3,2) * t54 + pkin(5)) + (g(1) * (rSges(3,1) * t54 - rSges(3,2) * t51) + g(2) * (-rSges(3,3) - pkin(6))) * t55) - m(4) * (g(1) * (t52 * rSges(4,2) + t59) + g(2) * (rSges(4,1) * t67 + rSges(4,3) * t69 + t60) + g(3) * (t51 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t54 + t70) + (g(1) * (rSges(4,1) * t54 + rSges(4,3) * t51) + g(2) * (-rSges(4,2) - pkin(6))) * t55) - m(5) * (g(1) * t57 + g(2) * t58 + g(3) * (-t51 * rSges(5,2) + (-rSges(5,1) - qJ(3)) * t54 + t61) + (g(1) * t56 + g(2) * (rSges(5,3) - pkin(6))) * t55 + (g(1) * (-rSges(5,3) - qJ(4)) + g(2) * t56) * t52) - m(6) * (g(1) * (pkin(4) * t66 - t52 * qJ(4) + (-t50 * t52 + t51 * t65) * rSges(6,1) + (-t50 * t66 - t68) * rSges(6,2) + t57) + g(2) * (pkin(4) * t69 - t55 * pkin(6) + (t50 * t55 + t51 * t68) * rSges(6,1) + (-t50 * t69 + t65) * rSges(6,2) + t58) + g(3) * (t72 * t51 + t61) + (g(3) * (-rSges(6,1) * t53 + rSges(6,2) * t50 - pkin(4) - qJ(3)) + (g(1) * t55 + g(2) * t52) * t72) * t54);
U = t1;
