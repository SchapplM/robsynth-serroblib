% Calculate potential energy for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:10
% EndTime: 2019-12-05 17:59:10
% DurationCPUTime: 0.20s
% Computational Cost: add. (96->59), mult. (103->63), div. (0->0), fcn. (79->6), ass. (0->22)
t54 = rSges(6,1) + pkin(4);
t53 = pkin(2) + pkin(5);
t42 = -pkin(7) - pkin(6);
t38 = sin(qJ(3));
t52 = pkin(3) * t38;
t51 = rSges(4,3) + pkin(6);
t50 = rSges(5,3) - t42;
t49 = rSges(6,3) + qJ(5) - t42;
t39 = sin(qJ(1));
t41 = cos(qJ(1));
t48 = t41 * pkin(1) + t39 * qJ(2);
t40 = cos(qJ(3));
t47 = t40 * pkin(3) + t53;
t46 = rSges(4,1) * t38 + rSges(4,2) * t40;
t37 = qJ(3) + qJ(4);
t30 = sin(t37);
t31 = cos(t37);
t45 = rSges(6,2) * t31 + t54 * t30 + t52;
t33 = t39 * pkin(1);
t44 = g(1) * t48 + g(2) * t33;
t43 = rSges(5,1) * t30 + rSges(5,2) * t31 + t52;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t41 - rSges(2,2) * t39) + g(2) * (rSges(2,1) * t39 + rSges(2,2) * t41) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t41 + rSges(3,3) * t39 + t48) + g(2) * (-rSges(3,2) * t39 + t33 + (-rSges(3,3) - qJ(2)) * t41) + g(3) * (pkin(5) + rSges(3,1))) - m(4) * (g(3) * (rSges(4,1) * t40 - rSges(4,2) * t38 + t53) + (g(1) * t46 + g(2) * t51) * t39 + (g(1) * t51 + g(2) * (-qJ(2) - t46)) * t41 + t44) - m(5) * (g(3) * (rSges(5,1) * t31 - rSges(5,2) * t30 + t47) + (g(1) * t43 + g(2) * t50) * t39 + (g(1) * t50 + g(2) * (-qJ(2) - t43)) * t41 + t44) - m(6) * (g(3) * (-rSges(6,2) * t30 + t54 * t31 + t47) + (g(1) * t45 + g(2) * t49) * t39 + (g(1) * t49 + g(2) * (-qJ(2) - t45)) * t41 + t44);
U = t1;
