% Calculate potential energy for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:04
% EndTime: 2019-12-31 21:49:04
% DurationCPUTime: 0.23s
% Computational Cost: add. (139->57), mult. (92->62), div. (0->0), fcn. (68->8), ass. (0->23)
t61 = rSges(6,1) + pkin(4);
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t60 = rSges(5,1) * t48 - rSges(5,2) * t46;
t59 = rSges(6,3) + qJ(5);
t58 = pkin(5) + pkin(6);
t45 = qJ(1) + qJ(2);
t40 = sin(t45);
t47 = sin(qJ(1));
t43 = t47 * pkin(1);
t55 = pkin(2) * t40 + t43;
t41 = cos(t45);
t49 = cos(qJ(1));
t44 = t49 * pkin(1);
t54 = pkin(2) * t41 + t44;
t53 = pkin(7) + t58;
t42 = qJ(3) + t45;
t38 = sin(t42);
t52 = t38 * pkin(3) + t55;
t39 = cos(t42);
t51 = t39 * pkin(3) + t38 * pkin(8) + t54;
t50 = t59 * t46 + t61 * t48;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t49 - t47 * rSges(2,2)) + g(2) * (t47 * rSges(2,1) + rSges(2,2) * t49) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t41 - rSges(3,2) * t40 + t44) + g(2) * (rSges(3,1) * t40 + rSges(3,2) * t41 + t43) + g(3) * (rSges(3,3) + t58)) - m(4) * (g(1) * (rSges(4,1) * t39 - rSges(4,2) * t38 + t54) + g(2) * (rSges(4,1) * t38 + rSges(4,2) * t39 + t55) + g(3) * (rSges(4,3) + t53)) - m(5) * (g(1) * (rSges(5,3) * t38 + t51) + g(2) * (t60 * t38 + t52) + g(3) * (rSges(5,1) * t46 + rSges(5,2) * t48 + t53) + (g(1) * t60 + g(2) * (-rSges(5,3) - pkin(8))) * t39) - m(6) * (g(1) * t51 + g(2) * t52 + g(3) * (t61 * t46 - t59 * t48 + t53) + (g(1) * rSges(6,2) + g(2) * t50) * t38 + (g(1) * t50 + g(2) * (-rSges(6,2) - pkin(8))) * t39);
U = t1;
