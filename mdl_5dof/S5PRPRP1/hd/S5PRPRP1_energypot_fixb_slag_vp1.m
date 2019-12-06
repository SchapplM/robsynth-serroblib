% Calculate potential energy for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:14
% EndTime: 2019-12-05 15:28:15
% DurationCPUTime: 0.24s
% Computational Cost: add. (137->61), mult. (104->68), div. (0->0), fcn. (80->8), ass. (0->25)
t63 = rSges(6,1) + pkin(4);
t45 = pkin(8) + qJ(4);
t38 = sin(t45);
t40 = cos(t45);
t62 = rSges(5,1) * t40 - rSges(5,2) * t38;
t61 = rSges(6,3) + qJ(5);
t58 = pkin(5) + qJ(1);
t49 = cos(pkin(8));
t37 = t49 * pkin(3) + pkin(2);
t46 = pkin(7) + qJ(2);
t41 = cos(t46);
t50 = cos(pkin(7));
t44 = t50 * pkin(1);
t57 = t41 * t37 + t44;
t56 = rSges(4,3) + qJ(3);
t47 = sin(pkin(8));
t55 = t47 * pkin(3) + t58;
t39 = sin(t46);
t48 = sin(pkin(7));
t43 = t48 * pkin(1);
t51 = -pkin(6) - qJ(3);
t54 = t39 * t37 + t41 * t51 + t43;
t53 = rSges(4,1) * t49 - rSges(4,2) * t47 + pkin(2);
t52 = t61 * t38 + t63 * t40;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t50 * rSges(2,1) - t48 * rSges(2,2)) + g(2) * (t48 * rSges(2,1) + t50 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t41 * rSges(3,1) - t39 * rSges(3,2) + t44) + g(2) * (t39 * rSges(3,1) + t41 * rSges(3,2) + t43) + g(3) * (rSges(3,3) + t58)) - m(4) * (g(1) * t44 + g(2) * t43 + g(3) * (t47 * rSges(4,1) + t49 * rSges(4,2) + t58) + (g(1) * t53 - g(2) * t56) * t41 + (g(1) * t56 + g(2) * t53) * t39) - m(5) * (g(1) * (t62 * t41 + t57) + g(2) * (-t41 * rSges(5,3) + t54) + g(3) * (t38 * rSges(5,1) + t40 * rSges(5,2) + t55) + (g(1) * (rSges(5,3) - t51) + g(2) * t62) * t39) - m(6) * (g(1) * t57 + g(2) * t54 + g(3) * (t63 * t38 - t61 * t40 + t55) + (-g(2) * rSges(6,2) + g(1) * t52) * t41 + (g(1) * (rSges(6,2) - t51) + g(2) * t52) * t39);
U = t1;
