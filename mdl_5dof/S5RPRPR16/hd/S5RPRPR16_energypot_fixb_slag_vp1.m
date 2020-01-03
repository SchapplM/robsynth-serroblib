% Calculate potential energy for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR16_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR16_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:40
% EndTime: 2019-12-31 18:38:40
% DurationCPUTime: 0.29s
% Computational Cost: add. (88->72), mult. (130->90), div. (0->0), fcn. (110->6), ass. (0->22)
t58 = rSges(6,3) + pkin(7);
t59 = pkin(2) + pkin(5);
t41 = sin(qJ(3));
t42 = sin(qJ(1));
t57 = t42 * t41;
t44 = cos(qJ(3));
t56 = t42 * t44;
t40 = sin(qJ(5));
t45 = cos(qJ(1));
t55 = t45 * t40;
t43 = cos(qJ(5));
t54 = t45 * t43;
t53 = t45 * pkin(1) + t42 * qJ(2);
t36 = t42 * pkin(1);
t52 = t42 * pkin(6) + t36;
t51 = t44 * qJ(4);
t50 = t45 * t51 + t52;
t49 = t45 * pkin(6) + t53;
t48 = t44 * pkin(3) + t41 * qJ(4) + t59;
t47 = pkin(3) * t57 + t49;
t46 = -rSges(5,2) * t41 - rSges(5,3) * t44;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t45 * rSges(2,1) - t42 * rSges(2,2)) + g(2) * (t42 * rSges(2,1) + t45 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (-t45 * rSges(3,2) + t42 * rSges(3,3) + t53) + g(2) * (-t42 * rSges(3,2) + t36 + (-rSges(3,3) - qJ(2)) * t45) + g(3) * (pkin(5) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t57 + rSges(4,2) * t56 + t49) + g(2) * (t42 * rSges(4,3) + t52) + g(3) * (t44 * rSges(4,1) - t41 * rSges(4,2) + t59) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t41 - rSges(4,2) * t44 - qJ(2))) * t45) - m(5) * (g(1) * t47 + g(2) * t50 + g(3) * (-t44 * rSges(5,2) + t41 * rSges(5,3) + t48) + (g(1) * (t46 - t51) + g(2) * rSges(5,1)) * t42 + (g(1) * rSges(5,1) + g(2) * (-t41 * pkin(3) - qJ(2) - t46)) * t45) - m(6) * (g(1) * (t45 * pkin(4) - t42 * t51 + (-t40 * t56 + t54) * rSges(6,1) + (-t43 * t56 - t55) * rSges(6,2) + t47) + g(2) * (t42 * pkin(4) - t45 * qJ(2) + (t42 * t43 + t44 * t55) * rSges(6,1) + (-t42 * t40 + t44 * t54) * rSges(6,2) + t50) + g(3) * (t58 * t44 + t48) + (g(3) * (rSges(6,1) * t40 + rSges(6,2) * t43) + g(1) * t58 * t42 + g(2) * (-pkin(3) - t58) * t45) * t41);
U = t1;
