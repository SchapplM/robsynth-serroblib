% Calculate potential energy for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:43
% EndTime: 2019-12-31 17:40:43
% DurationCPUTime: 0.24s
% Computational Cost: add. (128->65), mult. (115->75), div. (0->0), fcn. (91->6), ass. (0->22)
t67 = rSges(6,1) + pkin(4);
t51 = pkin(7) + qJ(2);
t46 = sin(t51);
t54 = sin(qJ(3));
t66 = t46 * t54;
t55 = cos(qJ(3));
t65 = t46 * t55;
t64 = pkin(5) + qJ(1);
t52 = sin(pkin(7));
t48 = t52 * pkin(1);
t63 = t46 * pkin(2) + t48;
t62 = qJ(4) * t54;
t61 = -rSges(6,3) - qJ(5);
t60 = t54 * pkin(3) + t64;
t47 = cos(t51);
t53 = cos(pkin(7));
t49 = t53 * pkin(1);
t59 = t47 * pkin(2) + t46 * pkin(6) + t49;
t58 = pkin(3) * t65 + t46 * t62 + t63;
t57 = t59 + (pkin(3) * t55 + t62) * t47;
t56 = rSges(6,2) * t54 + t67 * t55;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t53 - rSges(2,2) * t52) + g(2) * (rSges(2,1) * t52 + rSges(2,2) * t53) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t47 - rSges(3,2) * t46 + t49) + g(2) * (rSges(3,1) * t46 + rSges(3,2) * t47 + t48) + g(3) * (rSges(3,3) + t64)) - m(4) * (g(1) * (t46 * rSges(4,3) + t59) + g(2) * (rSges(4,1) * t65 - rSges(4,2) * t66 + t63) + g(3) * (t54 * rSges(4,1) + rSges(4,2) * t55 + t64) + (g(1) * (rSges(4,1) * t55 - rSges(4,2) * t54) + g(2) * (-rSges(4,3) - pkin(6))) * t47) - m(5) * (g(1) * (t46 * rSges(5,2) + t57) + g(2) * (rSges(5,1) * t65 + rSges(5,3) * t66 + t58) + g(3) * (t54 * rSges(5,1) + (-rSges(5,3) - qJ(4)) * t55 + t60) + (g(1) * (rSges(5,1) * t55 + rSges(5,3) * t54) + g(2) * (-rSges(5,2) - pkin(6))) * t47) - m(6) * (g(1) * t57 + g(2) * t58 + g(3) * ((-rSges(6,2) - qJ(4)) * t55 + t67 * t54 + t60) + (g(1) * t61 + g(2) * t56) * t46 + (g(1) * t56 + g(2) * (-pkin(6) - t61)) * t47);
U = t1;
