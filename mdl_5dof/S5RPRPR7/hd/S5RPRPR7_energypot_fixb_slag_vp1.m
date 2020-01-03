% Calculate potential energy for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR7_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:50
% EndTime: 2019-12-31 18:18:50
% DurationCPUTime: 0.28s
% Computational Cost: add. (149->70), mult. (117->84), div. (0->0), fcn. (97->10), ass. (0->30)
t73 = rSges(6,3) + pkin(7);
t72 = rSges(4,3) + pkin(6);
t51 = qJ(3) + pkin(9);
t44 = sin(t51);
t70 = rSges(5,2) * t44;
t52 = qJ(1) + pkin(8);
t45 = sin(t52);
t54 = sin(qJ(5));
t69 = t45 * t54;
t57 = cos(qJ(5));
t68 = t45 * t57;
t46 = cos(t51);
t47 = cos(t52);
t67 = t47 * t46;
t66 = t47 * t54;
t65 = t47 * t57;
t64 = pkin(5) + qJ(2);
t58 = cos(qJ(3));
t43 = t58 * pkin(3) + pkin(2);
t59 = cos(qJ(1));
t50 = t59 * pkin(1);
t63 = t47 * t43 + t50;
t55 = sin(qJ(3));
t62 = t55 * pkin(3) + t64;
t56 = sin(qJ(1));
t49 = t56 * pkin(1);
t53 = -qJ(4) - pkin(6);
t61 = t45 * t43 + t47 * t53 + t49;
t60 = rSges(4,1) * t58 - rSges(4,2) * t55 + pkin(2);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t59 * rSges(2,1) - t56 * rSges(2,2)) + g(2) * (t56 * rSges(2,1) + t59 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t47 * rSges(3,1) - t45 * rSges(3,2) + t50) + g(2) * (t45 * rSges(3,1) + t47 * rSges(3,2) + t49) + g(3) * (rSges(3,3) + t64)) - m(4) * (g(1) * t50 + g(2) * t49 + g(3) * (t55 * rSges(4,1) + t58 * rSges(4,2) + t64) + (g(1) * t60 - g(2) * t72) * t47 + (g(1) * t72 + g(2) * t60) * t45) - m(5) * (g(1) * (rSges(5,1) * t67 - t47 * t70 + t63) + g(2) * (-t47 * rSges(5,3) + t61) + g(3) * (t44 * rSges(5,1) + t46 * rSges(5,2) + t62) + (g(1) * (rSges(5,3) - t53) + g(2) * (rSges(5,1) * t46 - t70)) * t45) - m(6) * (g(1) * (pkin(4) * t67 - t45 * t53 + (t46 * t65 + t69) * rSges(6,1) + (-t46 * t66 + t68) * rSges(6,2) + t63) + g(2) * (t45 * t46 * pkin(4) + (t46 * t68 - t66) * rSges(6,1) + (-t46 * t69 - t65) * rSges(6,2) + t61) + g(3) * (-t73 * t46 + t62) + (g(3) * (rSges(6,1) * t57 - rSges(6,2) * t54 + pkin(4)) + (g(1) * t47 + g(2) * t45) * t73) * t44);
U = t1;
