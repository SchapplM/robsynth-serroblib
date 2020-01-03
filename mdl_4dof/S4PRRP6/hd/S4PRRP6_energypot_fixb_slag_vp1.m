% Calculate potential energy for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP6_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:10
% EndTime: 2019-12-31 16:30:10
% DurationCPUTime: 0.24s
% Computational Cost: add. (79->59), mult. (138->79), div. (0->0), fcn. (130->6), ass. (0->23)
t62 = rSges(5,1) + pkin(3);
t61 = rSges(5,3) + qJ(4);
t46 = sin(pkin(6));
t49 = sin(qJ(2));
t60 = t46 * t49;
t51 = cos(qJ(2));
t59 = t46 * t51;
t47 = cos(pkin(6));
t58 = t47 * t49;
t48 = sin(qJ(3));
t57 = t48 * t51;
t50 = cos(qJ(3));
t56 = t50 * t51;
t55 = t47 * pkin(1) + t46 * pkin(4);
t54 = t49 * pkin(2) + qJ(1);
t53 = t47 * t51 * pkin(2) + pkin(5) * t58 + t55;
t43 = t46 * pkin(1);
t52 = pkin(2) * t59 - t47 * pkin(4) + pkin(5) * t60 + t43;
t36 = t46 * t48 + t47 * t56;
t35 = -t46 * t50 + t47 * t57;
t34 = t46 * t56 - t47 * t48;
t33 = t46 * t57 + t47 * t50;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t47 - rSges(2,2) * t46) + g(2) * (rSges(2,1) * t46 + rSges(2,2) * t47) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t46 * rSges(3,3) + t55) + g(2) * (rSges(3,1) * t59 - rSges(3,2) * t60 + t43) + g(3) * (t49 * rSges(3,1) + rSges(3,2) * t51 + qJ(1)) + (g(1) * (rSges(3,1) * t51 - rSges(3,2) * t49) + g(2) * (-rSges(3,3) - pkin(4))) * t47) - m(4) * (g(1) * (rSges(4,1) * t36 - rSges(4,2) * t35 + rSges(4,3) * t58 + t53) + g(2) * (rSges(4,1) * t34 - rSges(4,2) * t33 + rSges(4,3) * t60 + t52) + g(3) * ((-rSges(4,3) - pkin(5)) * t51 + (rSges(4,1) * t50 - rSges(4,2) * t48) * t49 + t54)) - m(5) * (g(1) * (t61 * t35 + t62 * t36 + t53) + g(2) * (t61 * t33 + t62 * t34 + t52) + g(3) * (t54 + (-rSges(5,2) - pkin(5)) * t51) + (g(3) * (t61 * t48 + t62 * t50) + (g(1) * t47 + g(2) * t46) * rSges(5,2)) * t49);
U = t1;
