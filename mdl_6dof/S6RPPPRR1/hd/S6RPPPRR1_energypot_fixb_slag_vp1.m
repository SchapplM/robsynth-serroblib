% Calculate potential energy for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:47
% EndTime: 2019-03-09 01:29:48
% DurationCPUTime: 0.28s
% Computational Cost: add. (159->78), mult. (127->89), div. (0->0), fcn. (103->8), ass. (0->25)
t65 = rSges(7,3) + pkin(8);
t47 = sin(qJ(5));
t64 = pkin(5) * t47;
t46 = sin(qJ(6));
t62 = t46 * t47;
t49 = cos(qJ(6));
t61 = t47 * t49;
t60 = pkin(6) + qJ(2);
t45 = qJ(1) + pkin(9);
t41 = sin(t45);
t48 = sin(qJ(1));
t43 = t48 * pkin(1);
t59 = t41 * pkin(2) + t43;
t58 = pkin(3) + t60;
t57 = t41 * qJ(4) + t59;
t42 = cos(t45);
t51 = cos(qJ(1));
t44 = t51 * pkin(1);
t56 = t42 * pkin(2) + t41 * qJ(3) + t44;
t55 = pkin(4) + t58;
t54 = t42 * pkin(7) + t57;
t53 = t42 * qJ(4) + t56;
t50 = cos(qJ(5));
t52 = rSges(6,1) * t47 + rSges(6,2) * t50;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t51 - t48 * rSges(2,2)) + g(2) * (t48 * rSges(2,1) + rSges(2,2) * t51) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t42 - rSges(3,2) * t41 + t44) + g(2) * (rSges(3,1) * t41 + rSges(3,2) * t42 + t43) + g(3) * (rSges(3,3) + t60)) - m(4) * (g(1) * (-rSges(4,2) * t42 + rSges(4,3) * t41 + t56) + g(2) * (-rSges(4,2) * t41 + (-rSges(4,3) - qJ(3)) * t42 + t59) + g(3) * (rSges(4,1) + t60)) - m(5) * (g(1) * (rSges(5,2) * t41 + rSges(5,3) * t42 + t53) + g(2) * (rSges(5,3) * t41 + (-rSges(5,2) - qJ(3)) * t42 + t57) + g(3) * (rSges(5,1) + t58)) - m(6) * (g(1) * t53 + g(2) * t54 + g(3) * (rSges(6,1) * t50 - rSges(6,2) * t47 + t55) + (g(1) * t52 + g(2) * (rSges(6,3) - qJ(3))) * t42 + (g(1) * (-rSges(6,3) - pkin(7)) + g(2) * t52) * t41) - m(7) * (g(1) * (t42 * t64 - t41 * pkin(7) + (-t41 * t46 + t42 * t61) * rSges(7,1) + (-t41 * t49 - t42 * t62) * rSges(7,2) + t53) + g(2) * (t41 * t64 - t42 * qJ(3) + (t41 * t61 + t42 * t46) * rSges(7,1) + (-t41 * t62 + t42 * t49) * rSges(7,2) + t54) + g(3) * (t65 * t47 + t55) + (g(3) * (rSges(7,1) * t49 - rSges(7,2) * t46 + pkin(5)) - (g(1) * t42 + g(2) * t41) * t65) * t50);
U  = t1;
