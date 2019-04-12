% Calculate potential energy for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR14V3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(1,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_energypot_fixb_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14V3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:17
% EndTime: 2019-04-12 15:03:18
% DurationCPUTime: 0.19s
% Computational Cost: add. (103->76), mult. (227->117), div. (0->0), fcn. (240->10), ass. (0->34)
t50 = sin(qJ(4));
t51 = sin(qJ(2));
t68 = t50 * t51;
t52 = sin(qJ(1));
t67 = t51 * t52;
t55 = cos(qJ(4));
t66 = t51 * t55;
t57 = cos(qJ(1));
t65 = t51 * t57;
t56 = cos(qJ(2));
t64 = t52 * t56;
t63 = t57 * t50;
t62 = t57 * t55;
t61 = qJ(3) * t51;
t60 = t56 * qJ(3);
t59 = rSges(3,1) * t56 - rSges(3,2) * t51;
t58 = rSges(4,1) * t56 + rSges(4,3) * t51;
t54 = cos(qJ(5));
t53 = cos(qJ(6));
t49 = sin(qJ(5));
t48 = sin(qJ(6));
t47 = t57 * t61;
t46 = t52 * t61;
t45 = t52 * t50 + t56 * t62;
t44 = -t52 * t55 + t56 * t63;
t43 = t55 * t64 - t63;
t42 = t50 * t64 + t62;
t41 = -t56 * t49 + t54 * t66;
t40 = t49 * t66 + t56 * t54;
t39 = t45 * t54 + t49 * t65;
t38 = t45 * t49 - t54 * t65;
t37 = t43 * t54 + t49 * t67;
t36 = t43 * t49 - t54 * t67;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t57 * rSges(2,1) - t52 * rSges(2,2)) + g(2) * (t52 * rSges(2,1) + t57 * rSges(2,2)) + g(3) * rSges(2,3)) - m(3) * (g(1) * (t52 * rSges(3,3) + t59 * t57) + g(2) * (-t57 * rSges(3,3) + t59 * t52) + g(3) * (t51 * rSges(3,1) + t56 * rSges(3,2))) - m(4) * (g(1) * (t52 * rSges(4,2) + t58 * t57 + t47) + g(2) * (-t57 * rSges(4,2) + t58 * t52 + t46) + g(3) * (t51 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t56)) - m(5) * (g(1) * (t45 * rSges(5,1) - t44 * rSges(5,2) + rSges(5,3) * t65 + t47) + g(2) * (t43 * rSges(5,1) - t42 * rSges(5,2) + rSges(5,3) * t67 + t46) + g(3) * ((-rSges(5,3) - qJ(3)) * t56 + (rSges(5,1) * t55 - rSges(5,2) * t50) * t51)) - m(6) * (g(1) * (t39 * rSges(6,1) - t38 * rSges(6,2) + t44 * rSges(6,3) + t47) + g(2) * (t37 * rSges(6,1) - t36 * rSges(6,2) + t42 * rSges(6,3) + t46) + g(3) * (t41 * rSges(6,1) - t40 * rSges(6,2) + rSges(6,3) * t68 - t60)) - m(7) * (g(1) * (t47 + (t39 * t53 + t44 * t48) * rSges(7,1) + (-t39 * t48 + t44 * t53) * rSges(7,2) + t38 * rSges(7,3)) + g(2) * (t46 + (t37 * t53 + t42 * t48) * rSges(7,1) + (-t37 * t48 + t42 * t53) * rSges(7,2) + t36 * rSges(7,3)) + g(3) * (-t60 + (t41 * t53 + t48 * t68) * rSges(7,1) + (-t41 * t48 + t53 * t68) * rSges(7,2) + t40 * rSges(7,3)));
U  = t1;
