% Calculate potential energy for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR5_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:20
% EndTime: 2019-12-05 15:53:21
% DurationCPUTime: 0.45s
% Computational Cost: add. (152->90), mult. (180->116), div. (0->0), fcn. (168->10), ass. (0->31)
t58 = -pkin(6) - qJ(3);
t78 = rSges(5,3) - t58;
t77 = rSges(6,3) + pkin(7) - t58;
t55 = sin(pkin(8));
t57 = cos(pkin(8));
t76 = g(1) * t57 + g(2) * t55;
t75 = rSges(4,3) + qJ(3);
t56 = cos(pkin(9));
t44 = t56 * pkin(3) + pkin(2);
t59 = sin(qJ(2));
t72 = rSges(3,2) * t59;
t54 = sin(pkin(9));
t71 = t55 * t54;
t60 = cos(qJ(2));
t70 = t55 * t60;
t69 = t57 * t54;
t68 = t57 * t60;
t53 = pkin(9) + qJ(4);
t46 = cos(t53);
t40 = pkin(4) * t46 + t44;
t67 = t60 * t40;
t66 = t60 * t44;
t63 = t57 * pkin(1) + t55 * pkin(5);
t49 = t55 * pkin(1);
t61 = -t57 * pkin(5) + t49;
t47 = qJ(5) + t53;
t45 = sin(t53);
t43 = cos(t47);
t42 = sin(t47);
t41 = t54 * pkin(3) + pkin(4) * t45;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t57 * rSges(2,1) - t55 * rSges(2,2)) + g(2) * (t55 * rSges(2,1) + t57 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t55 * rSges(3,3) + t63) + g(2) * (rSges(3,1) * t70 - t55 * t72 + t49) + g(3) * (t59 * rSges(3,1) + t60 * rSges(3,2) + qJ(1)) + (g(1) * (rSges(3,1) * t60 - t72) + g(2) * (-rSges(3,3) - pkin(5))) * t57) - m(4) * (g(1) * (pkin(2) * t68 + (t56 * t68 + t71) * rSges(4,1) + (-t54 * t68 + t55 * t56) * rSges(4,2) + t63) + g(2) * (pkin(2) * t70 + (t56 * t70 - t69) * rSges(4,1) + (-t54 * t70 - t57 * t56) * rSges(4,2) + t61) + g(3) * (-t75 * t60 + qJ(1)) + (g(3) * (rSges(4,1) * t56 - rSges(4,2) * t54 + pkin(2)) + t76 * t75) * t59) - m(5) * (g(1) * (t57 * t66 + pkin(3) * t71 + (t55 * t45 + t46 * t68) * rSges(5,1) + (-t45 * t68 + t55 * t46) * rSges(5,2) + t63) + g(2) * (t55 * t66 - pkin(3) * t69 + (-t57 * t45 + t46 * t70) * rSges(5,1) + (-t45 * t70 - t57 * t46) * rSges(5,2) + t61) + g(3) * (-t78 * t60 + qJ(1)) + (g(3) * (rSges(5,1) * t46 - rSges(5,2) * t45 + t44) + t76 * t78) * t59) - m(6) * (g(1) * (t57 * t67 + t55 * t41 + (t55 * t42 + t43 * t68) * rSges(6,1) + (-t42 * t68 + t55 * t43) * rSges(6,2) + t63) + g(2) * (t55 * t67 - t57 * t41 + (-t57 * t42 + t43 * t70) * rSges(6,1) + (-t42 * t70 - t57 * t43) * rSges(6,2) + t61) + g(3) * (-t77 * t60 + qJ(1)) + (g(3) * (rSges(6,1) * t43 - rSges(6,2) * t42 + t40) + t76 * t77) * t59);
U = t1;
