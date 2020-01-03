% Calculate potential energy for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:34
% EndTime: 2019-12-31 19:52:34
% DurationCPUTime: 0.23s
% Computational Cost: add. (114->57), mult. (94->62), div. (0->0), fcn. (70->6), ass. (0->20)
t55 = rSges(6,1) + pkin(4);
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t54 = rSges(5,1) * t40 + rSges(5,2) * t42;
t53 = rSges(6,3) + qJ(5);
t52 = pkin(5) + pkin(6);
t39 = qJ(1) + qJ(2);
t35 = sin(t39);
t41 = sin(qJ(1));
t37 = t41 * pkin(1);
t49 = t35 * pkin(2) + t37;
t48 = pkin(3) + t52;
t36 = cos(t39);
t43 = cos(qJ(1));
t38 = t43 * pkin(1);
t47 = t36 * pkin(2) + t35 * qJ(3) + t38;
t46 = t35 * pkin(7) + t49;
t45 = t36 * pkin(7) + t47;
t44 = t55 * t40 - t53 * t42;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t43 - t41 * rSges(2,2)) + g(2) * (t41 * rSges(2,1) + rSges(2,2) * t43) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t36 - rSges(3,2) * t35 + t38) + g(2) * (rSges(3,1) * t35 + rSges(3,2) * t36 + t37) + g(3) * (rSges(3,3) + t52)) - m(4) * (g(1) * (-rSges(4,2) * t36 + rSges(4,3) * t35 + t47) + g(2) * (-rSges(4,2) * t35 + (-rSges(4,3) - qJ(3)) * t36 + t49) + g(3) * (rSges(4,1) + t52)) - m(5) * (g(1) * (t54 * t35 + t45) + g(2) * (rSges(5,3) * t35 + t46) + g(3) * (rSges(5,1) * t42 - rSges(5,2) * t40 + t48) + (g(1) * rSges(5,3) + g(2) * (-qJ(3) - t54)) * t36) - m(6) * (g(1) * t45 + g(2) * t46 + g(3) * (t53 * t40 + t55 * t42 + t48) + (g(2) * rSges(6,2) + g(1) * t44) * t35 + (g(1) * rSges(6,2) + g(2) * (-qJ(3) - t44)) * t36);
U = t1;
