% Calculate potential energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energypot_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:23:46
% EndTime: 2019-12-05 18:23:47
% DurationCPUTime: 0.26s
% Computational Cost: add. (90->56), mult. (115->79), div. (0->0), fcn. (95->8), ass. (0->26)
t62 = rSges(4,1) + pkin(1);
t61 = rSges(6,3) + pkin(4);
t39 = qJ(2) + qJ(4);
t37 = sin(t39);
t38 = cos(t39);
t60 = rSges(5,1) * t38 - rSges(5,2) * t37;
t41 = sin(qJ(5));
t43 = sin(qJ(1));
t56 = t43 * t41;
t44 = cos(qJ(5));
t55 = t43 * t44;
t45 = cos(qJ(2));
t47 = pkin(2) + pkin(1);
t54 = t45 * t47;
t46 = cos(qJ(1));
t53 = t46 * t41;
t52 = t46 * t44;
t40 = -pkin(3) - qJ(3);
t51 = t46 * t40 + t43 * t54;
t50 = rSges(4,3) + qJ(3);
t42 = sin(qJ(2));
t49 = rSges(3,1) * t45 - rSges(3,2) * t42;
t48 = -rSges(4,2) * t42 + t62 * t45;
t36 = t42 * t47;
t34 = t46 * t54;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t46 * rSges(2,1) - t43 * rSges(2,2)) + g(2) * (t43 * rSges(2,1) + t46 * rSges(2,2)) + g(3) * rSges(2,3)) - m(3) * (g(1) * (t43 * rSges(3,3) + t49 * t46) + g(2) * (-t46 * rSges(3,3) + t49 * t43) + g(3) * (t42 * rSges(3,1) + t45 * rSges(3,2))) - m(4) * (g(3) * (t45 * rSges(4,2) + t62 * t42) + (g(1) * t48 - g(2) * t50) * t46 + (g(1) * t50 + g(2) * t48) * t43) - m(5) * (g(1) * (t60 * t46 + t34) + g(2) * (-t46 * rSges(5,3) + t51) + g(3) * (t37 * rSges(5,1) + t38 * rSges(5,2) + t36) + (g(1) * (rSges(5,3) - t40) + g(2) * t60) * t43) - m(6) * (g(1) * (t34 - t43 * t40 + (t38 * t52 + t56) * rSges(6,1) + (-t38 * t53 + t55) * rSges(6,2)) + g(2) * ((t38 * t55 - t53) * rSges(6,1) + (-t38 * t56 - t52) * rSges(6,2) + t51) + g(3) * (-t61 * t38 + t36) + (g(3) * (rSges(6,1) * t44 - rSges(6,2) * t41) + (g(1) * t46 + g(2) * t43) * t61) * t37);
U = t1;
