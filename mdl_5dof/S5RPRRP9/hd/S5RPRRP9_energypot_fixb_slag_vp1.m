% Calculate potential energy for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP9_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:32
% EndTime: 2019-12-31 18:48:32
% DurationCPUTime: 0.23s
% Computational Cost: add. (139->61), mult. (117->70), div. (0->0), fcn. (93->8), ass. (0->27)
t69 = rSges(6,1) + pkin(4);
t51 = pkin(8) + qJ(3);
t47 = qJ(4) + t51;
t42 = sin(t47);
t43 = cos(t47);
t68 = rSges(5,1) * t43 - rSges(5,2) * t42;
t67 = rSges(6,3) + qJ(5);
t52 = sin(pkin(8));
t66 = t52 * pkin(2) + pkin(5);
t53 = cos(pkin(8));
t44 = t53 * pkin(2) + pkin(1);
t54 = -pkin(6) - qJ(2);
t63 = rSges(4,3) - t54;
t46 = cos(t51);
t39 = pkin(3) * t46 + t44;
t50 = -pkin(7) + t54;
t55 = sin(qJ(1));
t56 = cos(qJ(1));
t62 = t55 * t39 + t56 * t50;
t61 = rSges(3,3) + qJ(2);
t45 = sin(t51);
t60 = pkin(3) * t45 + t66;
t59 = rSges(3,1) * t53 - rSges(3,2) * t52 + pkin(1);
t58 = rSges(4,1) * t46 - rSges(4,2) * t45 + t44;
t57 = t67 * t42 + t69 * t43;
t38 = t56 * t39;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t56 * rSges(2,1) - t55 * rSges(2,2)) + g(2) * (t55 * rSges(2,1) + t56 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(3) * (t52 * rSges(3,1) + t53 * rSges(3,2) + pkin(5)) + (g(1) * t59 - g(2) * t61) * t56 + (g(1) * t61 + g(2) * t59) * t55) - m(4) * (g(3) * (t45 * rSges(4,1) + t46 * rSges(4,2) + t66) + (g(1) * t58 - g(2) * t63) * t56 + (g(1) * t63 + g(2) * t58) * t55) - m(5) * (g(1) * (t68 * t56 + t38) + g(2) * (-t56 * rSges(5,3) + t62) + g(3) * (t42 * rSges(5,1) + t43 * rSges(5,2) + t60) + (g(1) * (rSges(5,3) - t50) + g(2) * t68) * t55) - m(6) * (g(1) * t38 + g(2) * t62 + g(3) * (t69 * t42 - t67 * t43 + t60) + (-g(2) * rSges(6,2) + g(1) * t57) * t56 + (g(1) * (rSges(6,2) - t50) + g(2) * t57) * t55);
U = t1;
