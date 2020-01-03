% Calculate potential energy for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR10_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:52
% EndTime: 2019-12-31 18:25:52
% DurationCPUTime: 0.16s
% Computational Cost: add. (120->64), mult. (126->71), div. (0->0), fcn. (120->8), ass. (0->26)
t66 = pkin(5) - pkin(6);
t65 = rSges(6,3) + pkin(7);
t50 = sin(qJ(3));
t51 = sin(qJ(1));
t64 = t51 * t50;
t54 = cos(qJ(1));
t63 = t54 * t50;
t62 = t54 * pkin(1) + t51 * qJ(2);
t61 = qJ(3) + pkin(8);
t60 = -qJ(4) + t66;
t47 = t51 * pkin(1);
t59 = -t54 * qJ(2) + t47;
t53 = cos(qJ(3));
t45 = t53 * pkin(3) + pkin(2);
t58 = pkin(3) * t64 + t54 * t45 + t62;
t57 = cos(t61);
t56 = sin(t61);
t49 = sin(qJ(5));
t52 = cos(qJ(5));
t55 = -rSges(6,1) * t52 + rSges(6,2) * t49 - pkin(4);
t39 = t51 * t45;
t38 = t51 * t53 - t63;
t37 = -t54 * t53 - t64;
t36 = -t51 * t57 + t54 * t56;
t35 = -t51 * t56 - t54 * t57;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t54 * rSges(2,1) - t51 * rSges(2,2)) + g(2) * (t51 * rSges(2,1) + t54 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t54 * rSges(3,1) + t51 * rSges(3,3) + t62) + g(2) * (t51 * rSges(3,1) + t47 + (-rSges(3,3) - qJ(2)) * t54) + g(3) * (pkin(5) + rSges(3,2))) - m(4) * (g(1) * (-t37 * rSges(4,1) + t38 * rSges(4,2) + t54 * pkin(2) + t62) + g(2) * (t38 * rSges(4,1) + t37 * rSges(4,2) + t51 * pkin(2) + t59) + g(3) * (-rSges(4,3) + t66)) - m(5) * (g(1) * (-t35 * rSges(5,1) - t36 * rSges(5,2) + t58) + g(2) * (-t36 * rSges(5,1) + t35 * rSges(5,2) + t39 + t47 + (-pkin(3) * t50 - qJ(2)) * t54) + g(3) * (-rSges(5,3) + t60)) - m(6) * (g(1) * t58 + g(2) * (-pkin(3) * t63 + t39 + t59) + g(3) * (-t49 * rSges(6,1) - t52 * rSges(6,2) + t60) + (g(1) * t65 + g(2) * t55) * t36 + (g(1) * t55 - g(2) * t65) * t35);
U = t1;
