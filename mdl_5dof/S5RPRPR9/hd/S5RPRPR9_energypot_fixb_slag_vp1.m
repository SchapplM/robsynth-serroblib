% Calculate potential energy for
% S5RPRPR9
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR9_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:38
% EndTime: 2019-12-31 18:23:38
% DurationCPUTime: 0.33s
% Computational Cost: add. (137->72), mult. (128->88), div. (0->0), fcn. (108->8), ass. (0->24)
t73 = rSges(6,3) + pkin(7);
t54 = qJ(1) + pkin(8);
t49 = sin(t54);
t56 = sin(qJ(3));
t71 = t49 * t56;
t59 = cos(qJ(3));
t70 = t49 * t59;
t55 = sin(qJ(5));
t69 = t55 * t56;
t58 = cos(qJ(5));
t68 = t56 * t58;
t67 = pkin(5) + qJ(2);
t57 = sin(qJ(1));
t52 = t57 * pkin(1);
t66 = t49 * pkin(2) + t52;
t65 = qJ(4) * t56;
t64 = t56 * pkin(3) + t67;
t50 = cos(t54);
t60 = cos(qJ(1));
t53 = t60 * pkin(1);
t63 = t50 * pkin(2) + t49 * pkin(6) + t53;
t62 = pkin(3) * t70 + t49 * t65 + t66;
t61 = t63 + (pkin(3) * t59 + t65) * t50;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t60 * rSges(2,1) - t57 * rSges(2,2)) + g(2) * (t57 * rSges(2,1) + t60 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t50 * rSges(3,1) - t49 * rSges(3,2) + t53) + g(2) * (t49 * rSges(3,1) + t50 * rSges(3,2) + t52) + g(3) * (rSges(3,3) + t67)) - m(4) * (g(1) * (t49 * rSges(4,3) + t63) + g(2) * (rSges(4,1) * t70 - rSges(4,2) * t71 + t66) + g(3) * (t56 * rSges(4,1) + t59 * rSges(4,2) + t67) + (g(1) * (rSges(4,1) * t59 - rSges(4,2) * t56) + g(2) * (-rSges(4,3) - pkin(6))) * t50) - m(5) * (g(1) * (t49 * rSges(5,1) + t61) + g(2) * (-rSges(5,2) * t70 + rSges(5,3) * t71 + t62) + g(3) * (-t56 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t59 + t64) + (g(1) * (-rSges(5,2) * t59 + rSges(5,3) * t56) + g(2) * (-rSges(5,1) - pkin(6))) * t50) - m(6) * (g(1) * (t49 * pkin(4) + (t49 * t58 + t50 * t69) * rSges(6,1) + (-t49 * t55 + t50 * t68) * rSges(6,2) + t61) + g(2) * (t62 + (rSges(6,1) * t69 + rSges(6,2) * t68) * t49 + (-rSges(6,1) * t58 + rSges(6,2) * t55 - pkin(4) - pkin(6)) * t50) + g(3) * (t73 * t56 + t64) + (g(3) * (-rSges(6,1) * t55 - rSges(6,2) * t58 - qJ(4)) + (g(1) * t50 + g(2) * t49) * t73) * t59);
U = t1;
