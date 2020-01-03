% Calculate potential energy for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP12_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:12
% EndTime: 2019-12-31 18:56:12
% DurationCPUTime: 0.29s
% Computational Cost: add. (96->72), mult. (143->87), div. (0->0), fcn. (127->6), ass. (0->28)
t74 = rSges(5,3) + pkin(7);
t73 = rSges(6,3) + qJ(5) + pkin(7);
t52 = sin(qJ(1));
t55 = cos(qJ(1));
t72 = -g(1) * t52 + g(2) * t55;
t71 = pkin(2) + pkin(5);
t54 = cos(qJ(3));
t67 = rSges(4,2) * t54;
t51 = sin(qJ(3));
t66 = t51 * t52;
t65 = t51 * t55;
t50 = sin(qJ(4));
t64 = t52 * t50;
t53 = cos(qJ(4));
t63 = t52 * t53;
t62 = t55 * t50;
t61 = t55 * t53;
t59 = t55 * pkin(1) + t52 * qJ(2);
t46 = t52 * pkin(1);
t58 = t52 * pkin(6) + t46;
t57 = t55 * pkin(6) + t59;
t56 = -t55 * qJ(2) + t58;
t43 = t53 * pkin(4) + pkin(3);
t42 = -t51 * t61 + t64;
t41 = t51 * t62 + t63;
t40 = t51 * t63 + t62;
t39 = -t51 * t64 + t61;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t55 * rSges(2,1) - t52 * rSges(2,2)) + g(2) * (t52 * rSges(2,1) + t55 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (-t55 * rSges(3,2) + t52 * rSges(3,3) + t59) + g(2) * (-t52 * rSges(3,2) + t46 + (-rSges(3,3) - qJ(2)) * t55) + g(3) * (pkin(5) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t66 + t52 * t67 + t57) + g(2) * (t52 * rSges(4,3) + t58) + g(3) * (t54 * rSges(4,1) - t51 * rSges(4,2) + t71) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t51 - qJ(2) - t67)) * t55) - m(5) * (g(1) * (t40 * rSges(5,1) + t39 * rSges(5,2) + pkin(3) * t66 + t57) + g(2) * (t42 * rSges(5,1) + t41 * rSges(5,2) - pkin(3) * t65 + t56) + g(3) * (t74 * t51 + t71) + (g(3) * (rSges(5,1) * t53 - rSges(5,2) * t50 + pkin(3)) + t72 * t74) * t54) - m(6) * (g(1) * (t40 * rSges(6,1) + t39 * rSges(6,2) + pkin(4) * t62 + t43 * t66 + t57) + g(2) * (t42 * rSges(6,1) + t41 * rSges(6,2) + pkin(4) * t64 - t43 * t65 + t56) + g(3) * (t73 * t51 + t71) + (g(3) * (rSges(6,1) * t53 - rSges(6,2) * t50 + t43) + t72 * t73) * t54);
U = t1;
