% Calculate potential energy for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP7_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:06
% EndTime: 2019-12-31 17:20:06
% DurationCPUTime: 0.23s
% Computational Cost: add. (79->59), mult. (138->77), div. (0->0), fcn. (130->6), ass. (0->22)
t60 = rSges(5,1) + pkin(3);
t59 = rSges(5,3) + qJ(4);
t46 = sin(qJ(2));
t58 = t46 * pkin(2) + pkin(4);
t47 = sin(qJ(1));
t57 = t46 * t47;
t50 = cos(qJ(1));
t56 = t46 * t50;
t49 = cos(qJ(2));
t55 = t47 * t49;
t54 = t49 * t50;
t53 = t50 * pkin(1) + t47 * pkin(5);
t52 = pkin(2) * t54 + pkin(6) * t56 + t53;
t43 = t47 * pkin(1);
t51 = pkin(2) * t55 - pkin(5) * t50 + pkin(6) * t57 + t43;
t48 = cos(qJ(3));
t45 = sin(qJ(3));
t35 = t47 * t45 + t48 * t54;
t34 = t45 * t54 - t47 * t48;
t33 = -t45 * t50 + t48 * t55;
t32 = t45 * t55 + t48 * t50;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t50 - t47 * rSges(2,2)) + g(2) * (t47 * rSges(2,1) + rSges(2,2) * t50) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t47 * rSges(3,3) + t53) + g(2) * (rSges(3,1) * t55 - rSges(3,2) * t57 + t43) + g(3) * (rSges(3,1) * t46 + rSges(3,2) * t49 + pkin(4)) + (g(1) * (rSges(3,1) * t49 - rSges(3,2) * t46) + g(2) * (-rSges(3,3) - pkin(5))) * t50) - m(4) * (g(1) * (t35 * rSges(4,1) - t34 * rSges(4,2) + rSges(4,3) * t56 + t52) + g(2) * (t33 * rSges(4,1) - t32 * rSges(4,2) + rSges(4,3) * t57 + t51) + g(3) * ((-rSges(4,3) - pkin(6)) * t49 + (rSges(4,1) * t48 - rSges(4,2) * t45) * t46 + t58)) - m(5) * (g(1) * (t59 * t34 + t60 * t35 + t52) + g(2) * (t59 * t32 + t60 * t33 + t51) + g(3) * (t58 + (-rSges(5,2) - pkin(6)) * t49) + (g(3) * (t59 * t45 + t60 * t48) + (g(1) * t50 + g(2) * t47) * rSges(5,2)) * t46);
U = t1;
