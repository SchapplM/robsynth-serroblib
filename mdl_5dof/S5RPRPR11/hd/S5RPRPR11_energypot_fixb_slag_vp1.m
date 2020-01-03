% Calculate potential energy for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR11_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:06
% EndTime: 2019-12-31 18:27:06
% DurationCPUTime: 0.26s
% Computational Cost: add. (137->70), mult. (144->86), div. (0->0), fcn. (126->8), ass. (0->27)
t50 = sin(pkin(8));
t70 = t50 * pkin(2) + pkin(5);
t69 = pkin(7) + rSges(6,3);
t49 = pkin(8) + qJ(3);
t46 = sin(t49);
t56 = cos(qJ(1));
t68 = t46 * t56;
t47 = cos(t49);
t67 = t56 * t47;
t51 = cos(pkin(8));
t44 = t51 * pkin(2) + pkin(1);
t52 = -pkin(6) - qJ(2);
t54 = sin(qJ(1));
t66 = t54 * t44 + t56 * t52;
t65 = qJ(4) * t46;
t64 = rSges(3,3) + qJ(2);
t63 = t46 * pkin(3) + t70;
t40 = t56 * t44;
t62 = pkin(3) * t67 + t56 * t65 + t40;
t61 = t66 + (pkin(3) * t47 + t65) * t54;
t53 = sin(qJ(5));
t55 = cos(qJ(5));
t60 = t46 * t55 - t47 * t53;
t59 = t46 * t53 + t47 * t55;
t58 = rSges(3,1) * t51 - rSges(3,2) * t50 + pkin(1);
t57 = rSges(6,1) * t59 + rSges(6,2) * t60 + t47 * pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t56 * rSges(2,1) - t54 * rSges(2,2)) + g(2) * (t54 * rSges(2,1) + t56 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(3) * (t50 * rSges(3,1) + t51 * rSges(3,2) + pkin(5)) + (g(1) * t58 - g(2) * t64) * t56 + (g(1) * t64 + g(2) * t58) * t54) - m(4) * (g(1) * (rSges(4,1) * t67 - rSges(4,2) * t68 + t40) + g(2) * (-t56 * rSges(4,3) + t66) + g(3) * (t46 * rSges(4,1) + t47 * rSges(4,2) + t70) + (g(1) * (rSges(4,3) - t52) + g(2) * (rSges(4,1) * t47 - rSges(4,2) * t46)) * t54) - m(5) * (g(1) * (rSges(5,1) * t67 + rSges(5,3) * t68 + t62) + g(2) * (-t56 * rSges(5,2) + t61) + g(3) * (t46 * rSges(5,1) + (-rSges(5,3) - qJ(4)) * t47 + t63) + (g(1) * (rSges(5,2) - t52) + g(2) * (rSges(5,1) * t47 + rSges(5,3) * t46)) * t54) - m(6) * (g(1) * t62 + g(2) * t61 + g(3) * (rSges(6,1) * t60 - t59 * rSges(6,2) + t46 * pkin(4) - t47 * qJ(4) + t63) + (g(1) * t57 + g(2) * t69) * t56 + (g(1) * (-t52 - t69) + g(2) * t57) * t54);
U = t1;
