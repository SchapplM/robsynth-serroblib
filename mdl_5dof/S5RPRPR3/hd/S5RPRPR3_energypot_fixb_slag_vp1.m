% Calculate potential energy for
% S5RPRPR3
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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:10
% EndTime: 2019-12-05 17:51:11
% DurationCPUTime: 0.32s
% Computational Cost: add. (152->66), mult. (105->78), div. (0->0), fcn. (85->10), ass. (0->26)
t64 = rSges(6,3) + pkin(7);
t49 = sin(qJ(1));
t63 = t49 * pkin(1);
t46 = sin(pkin(9));
t61 = rSges(5,2) * t46;
t45 = qJ(1) + pkin(8);
t43 = qJ(3) + t45;
t40 = cos(t43);
t47 = cos(pkin(9));
t60 = t40 * t47;
t48 = sin(qJ(5));
t59 = t47 * t48;
t50 = cos(qJ(5));
t58 = t47 * t50;
t57 = pkin(5) + qJ(2);
t42 = cos(t45);
t51 = cos(qJ(1));
t44 = t51 * pkin(1);
t56 = pkin(2) * t42 + t44;
t55 = pkin(6) + t57;
t39 = sin(t43);
t54 = t40 * pkin(3) + t39 * qJ(4) + t56;
t41 = sin(t45);
t53 = -pkin(2) * t41 - t63;
t52 = t40 * qJ(4) + t53;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (-t49 * rSges(2,1) - t51 * rSges(2,2)) + g(3) * (t51 * rSges(2,1) - t49 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) + t57) + g(2) * (-t41 * rSges(3,1) - t42 * rSges(3,2) - t63) + g(3) * (t42 * rSges(3,1) - t41 * rSges(3,2) + t44)) - m(4) * (g(1) * (rSges(4,3) + t55) + g(2) * (-t39 * rSges(4,1) - t40 * rSges(4,2) + t53) + g(3) * (t40 * rSges(4,1) - t39 * rSges(4,2) + t56)) - m(5) * (g(1) * (t46 * rSges(5,1) + t47 * rSges(5,2) + t55) + g(2) * (t40 * rSges(5,3) + t52) + g(3) * (rSges(5,1) * t60 - t40 * t61 + t54) + (g(2) * (-rSges(5,1) * t47 - pkin(3) + t61) + g(3) * rSges(5,3)) * t39) - m(6) * (g(1) * (-t64 * t47 + t55) + g(2) * (t52 + (rSges(6,1) * t48 + rSges(6,2) * t50) * t40 + (-rSges(6,1) * t58 + rSges(6,2) * t59 - t47 * pkin(4) - pkin(3)) * t39) + g(3) * (pkin(4) * t60 + (t39 * t48 + t40 * t58) * rSges(6,1) + (t39 * t50 - t40 * t59) * rSges(6,2) + t54) + (g(1) * (rSges(6,1) * t50 - rSges(6,2) * t48 + pkin(4)) + (-g(2) * t39 + g(3) * t40) * t64) * t46);
U = t1;
