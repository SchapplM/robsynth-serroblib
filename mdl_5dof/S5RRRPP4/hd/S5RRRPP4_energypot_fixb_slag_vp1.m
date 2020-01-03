% Calculate potential energy for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:48
% EndTime: 2019-12-31 20:54:48
% DurationCPUTime: 0.23s
% Computational Cost: add. (139->61), mult. (117->70), div. (0->0), fcn. (93->8), ass. (0->27)
t69 = rSges(6,1) + pkin(4);
t51 = qJ(2) + qJ(3);
t45 = pkin(8) + t51;
t41 = sin(t45);
t42 = cos(t45);
t68 = rSges(5,1) * t42 - rSges(5,2) * t41;
t67 = rSges(6,3) + qJ(5);
t56 = -pkin(7) - pkin(6);
t66 = rSges(3,3) + pkin(6);
t52 = sin(qJ(2));
t65 = t52 * pkin(2) + pkin(5);
t54 = cos(qJ(2));
t44 = t54 * pkin(2) + pkin(1);
t62 = rSges(4,3) - t56;
t47 = cos(t51);
t39 = pkin(3) * t47 + t44;
t50 = -qJ(4) + t56;
t53 = sin(qJ(1));
t55 = cos(qJ(1));
t61 = t53 * t39 + t55 * t50;
t46 = sin(t51);
t60 = pkin(3) * t46 + t65;
t59 = rSges(3,1) * t54 - rSges(3,2) * t52 + pkin(1);
t58 = rSges(4,1) * t47 - rSges(4,2) * t46 + t44;
t57 = t67 * t41 + t69 * t42;
t38 = t55 * t39;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t55 * rSges(2,1) - t53 * rSges(2,2)) + g(2) * (t53 * rSges(2,1) + t55 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(3) * (t52 * rSges(3,1) + t54 * rSges(3,2) + pkin(5)) + (g(1) * t59 - g(2) * t66) * t55 + (g(1) * t66 + g(2) * t59) * t53) - m(4) * (g(3) * (t46 * rSges(4,1) + t47 * rSges(4,2) + t65) + (g(1) * t58 - g(2) * t62) * t55 + (g(1) * t62 + g(2) * t58) * t53) - m(5) * (g(1) * (t68 * t55 + t38) + g(2) * (-t55 * rSges(5,3) + t61) + g(3) * (t41 * rSges(5,1) + t42 * rSges(5,2) + t60) + (g(1) * (rSges(5,3) - t50) + g(2) * t68) * t53) - m(6) * (g(1) * t38 + g(2) * t61 + g(3) * (t69 * t41 - t67 * t42 + t60) + (-g(2) * rSges(6,2) + g(1) * t57) * t55 + (g(1) * (rSges(6,2) - t50) + g(2) * t57) * t53);
U = t1;
