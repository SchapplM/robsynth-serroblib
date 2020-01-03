% Calculate potential energy for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:23
% EndTime: 2019-12-31 17:34:23
% DurationCPUTime: 0.23s
% Computational Cost: add. (100->56), mult. (141->64), div. (0->0), fcn. (141->6), ass. (0->20)
t58 = -rSges(6,1) - pkin(4);
t57 = rSges(5,3) + pkin(6);
t56 = cos(qJ(3));
t55 = sin(qJ(3));
t54 = rSges(6,3) + qJ(5) + pkin(6);
t53 = qJ(1) - pkin(5);
t42 = cos(pkin(7));
t51 = sin(pkin(7));
t52 = t42 * pkin(1) + t51 * qJ(2);
t50 = t42 * pkin(2) + t52;
t39 = t51 * pkin(1);
t49 = t51 * pkin(2) - t42 * qJ(2) + t39;
t44 = sin(qJ(4));
t45 = cos(qJ(4));
t48 = -rSges(5,1) * t45 + rSges(5,2) * t44 - pkin(3);
t47 = rSges(6,2) * t44 + t58 * t45 - pkin(3);
t46 = g(1) * t50 + g(2) * t49;
t32 = t42 * t55 - t51 * t56;
t31 = -t42 * t56 - t51 * t55;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t42 * rSges(2,1) - t51 * rSges(2,2)) + g(2) * (t51 * rSges(2,1) + t42 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t42 * rSges(3,1) + t51 * rSges(3,3) + t52) + g(2) * (t51 * rSges(3,1) + t39 + (-rSges(3,3) - qJ(2)) * t42) + g(3) * (qJ(1) + rSges(3,2))) - m(4) * (g(1) * (-t31 * rSges(4,1) - t32 * rSges(4,2) + t50) + g(2) * (-t32 * rSges(4,1) + t31 * rSges(4,2) + t49) + g(3) * (-rSges(4,3) + t53)) - m(5) * (g(3) * (-t44 * rSges(5,1) - t45 * rSges(5,2) + t53) + (g(1) * t57 + g(2) * t48) * t32 + (g(1) * t48 - g(2) * t57) * t31 + t46) - m(6) * (g(3) * (-t45 * rSges(6,2) + t58 * t44 + t53) + (g(1) * t54 + g(2) * t47) * t32 + (g(1) * t47 - g(2) * t54) * t31 + t46);
U = t1;
