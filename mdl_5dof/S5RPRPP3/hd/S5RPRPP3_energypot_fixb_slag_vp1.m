% Calculate potential energy for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:16
% EndTime: 2019-12-31 18:12:16
% DurationCPUTime: 0.26s
% Computational Cost: add. (125->67), mult. (128->79), div. (0->0), fcn. (104->6), ass. (0->24)
t64 = rSges(6,3) + qJ(5);
t63 = -rSges(6,1) - pkin(4);
t47 = sin(pkin(7));
t62 = t47 * pkin(2) + pkin(5);
t46 = pkin(7) + qJ(3);
t43 = sin(t46);
t51 = cos(qJ(1));
t61 = t43 * t51;
t44 = cos(t46);
t60 = t44 * t51;
t48 = cos(pkin(7));
t41 = t48 * pkin(2) + pkin(1);
t49 = -pkin(6) - qJ(2);
t50 = sin(qJ(1));
t59 = t50 * t41 + t51 * t49;
t58 = qJ(4) * t43;
t57 = rSges(3,3) + qJ(2);
t56 = t43 * pkin(3) + t62;
t37 = t51 * t41;
t55 = pkin(3) * t60 + t51 * t58 + t37;
t54 = t59 + (pkin(3) * t44 + t58) * t50;
t53 = rSges(3,1) * t48 - rSges(3,2) * t47 + pkin(1);
t52 = rSges(6,2) * t43 + t44 * t64;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t51 * rSges(2,1) - t50 * rSges(2,2)) + g(2) * (t50 * rSges(2,1) + t51 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(3) * (t47 * rSges(3,1) + t48 * rSges(3,2) + pkin(5)) + (g(1) * t53 - g(2) * t57) * t51 + (g(1) * t57 + g(2) * t53) * t50) - m(4) * (g(1) * (rSges(4,1) * t60 - rSges(4,2) * t61 + t37) + g(2) * (-t51 * rSges(4,3) + t59) + g(3) * (t43 * rSges(4,1) + t44 * rSges(4,2) + t62) + (g(1) * (rSges(4,3) - t49) + g(2) * (rSges(4,1) * t44 - rSges(4,2) * t43)) * t50) - m(5) * (g(1) * (-rSges(5,2) * t60 + rSges(5,3) * t61 + t55) + g(2) * (-t51 * rSges(5,1) + t54) + g(3) * (-t43 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t44 + t56) + (g(1) * (rSges(5,1) - t49) + g(2) * (-rSges(5,2) * t44 + rSges(5,3) * t43)) * t50) - m(6) * (g(1) * t55 + g(2) * t54 + g(3) * ((-rSges(6,2) - qJ(4)) * t44 + t64 * t43 + t56) + (g(1) * t52 + g(2) * t63) * t51 + (g(1) * (-t49 - t63) + g(2) * t52) * t50);
U = t1;
