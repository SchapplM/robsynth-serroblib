% Calculate potential energy for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR2_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:15
% EndTime: 2019-12-05 16:17:16
% DurationCPUTime: 0.30s
% Computational Cost: add. (152->66), mult. (105->80), div. (0->0), fcn. (85->10), ass. (0->26)
t67 = rSges(6,3) + pkin(7);
t50 = sin(pkin(9));
t65 = rSges(5,2) * t50;
t49 = pkin(8) + qJ(2);
t46 = qJ(3) + t49;
t42 = sin(t46);
t52 = cos(pkin(9));
t64 = t42 * t52;
t54 = sin(qJ(5));
t63 = t52 * t54;
t55 = cos(qJ(5));
t62 = t52 * t55;
t61 = pkin(5) + qJ(1);
t44 = sin(t49);
t51 = sin(pkin(8));
t47 = t51 * pkin(1);
t60 = pkin(2) * t44 + t47;
t45 = cos(t49);
t53 = cos(pkin(8));
t48 = t53 * pkin(1);
t59 = pkin(2) * t45 + t48;
t58 = pkin(6) + t61;
t57 = t42 * pkin(3) + t60;
t43 = cos(t46);
t56 = t43 * pkin(3) + t42 * qJ(4) + t59;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t53 * rSges(2,1) - t51 * rSges(2,2)) + g(2) * (t51 * rSges(2,1) + t53 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t45 * rSges(3,1) - t44 * rSges(3,2) + t48) + g(2) * (t44 * rSges(3,1) + t45 * rSges(3,2) + t47) + g(3) * (rSges(3,3) + t61)) - m(4) * (g(1) * (t43 * rSges(4,1) - t42 * rSges(4,2) + t59) + g(2) * (t42 * rSges(4,1) + t43 * rSges(4,2) + t60) + g(3) * (rSges(4,3) + t58)) - m(5) * (g(1) * (t42 * rSges(5,3) + t56) + g(2) * (rSges(5,1) * t64 - t42 * t65 + t57) + g(3) * (t50 * rSges(5,1) + t52 * rSges(5,2) + t58) + (g(1) * (rSges(5,1) * t52 - t65) + g(2) * (-rSges(5,3) - qJ(4))) * t43) - m(6) * (g(1) * (t43 * t52 * pkin(4) + (t42 * t54 + t43 * t62) * rSges(6,1) + (t42 * t55 - t43 * t63) * rSges(6,2) + t56) + g(2) * (pkin(4) * t64 - t43 * qJ(4) + (t42 * t62 - t43 * t54) * rSges(6,1) + (-t42 * t63 - t43 * t55) * rSges(6,2) + t57) + g(3) * (-t67 * t52 + t58) + (g(3) * (rSges(6,1) * t55 - rSges(6,2) * t54 + pkin(4)) + (g(1) * t43 + g(2) * t42) * t67) * t50);
U = t1;
