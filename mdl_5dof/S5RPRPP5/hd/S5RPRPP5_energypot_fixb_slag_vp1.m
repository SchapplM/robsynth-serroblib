% Calculate potential energy for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP5_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:04
% EndTime: 2019-12-31 18:16:05
% DurationCPUTime: 0.27s
% Computational Cost: add. (83->65), mult. (117->73), div. (0->0), fcn. (93->4), ass. (0->19)
t53 = rSges(6,1) + pkin(4);
t54 = pkin(2) + pkin(5);
t40 = cos(qJ(3));
t52 = rSges(4,2) * t40;
t51 = rSges(6,2) * t40;
t38 = sin(qJ(3));
t39 = sin(qJ(1));
t50 = t38 * t39;
t41 = cos(qJ(1));
t49 = t41 * pkin(1) + t39 * qJ(2);
t34 = t39 * pkin(1);
t48 = t39 * pkin(6) + t34;
t47 = qJ(4) * t40;
t46 = -rSges(6,3) - qJ(5);
t45 = t41 * pkin(6) + t49;
t44 = t40 * pkin(3) + t38 * qJ(4) + t54;
t43 = rSges(5,1) * t38 - rSges(5,3) * t40;
t42 = g(1) * (pkin(3) * t50 + t45) + g(2) * (t41 * t47 + t48);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t41 * rSges(2,1) - t39 * rSges(2,2)) + g(2) * (t39 * rSges(2,1) + t41 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (-t41 * rSges(3,2) + t39 * rSges(3,3) + t49) + g(2) * (-t39 * rSges(3,2) + t34 + (-rSges(3,3) - qJ(2)) * t41) + g(3) * (pkin(5) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t50 + t39 * t52 + t45) + g(2) * (t39 * rSges(4,3) + t48) + g(3) * (t40 * rSges(4,1) - t38 * rSges(4,2) + t54) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t38 - qJ(2) - t52)) * t41) - m(5) * (g(3) * (t40 * rSges(5,1) + t38 * rSges(5,3) + t44) + (g(1) * (t43 - t47) + g(2) * rSges(5,2)) * t39 + (g(1) * rSges(5,2) + g(2) * (-pkin(3) * t38 - qJ(2) - t43)) * t41 + t42) - m(6) * (g(3) * (t38 * rSges(6,2) + t53 * t40 + t44) + (g(1) * (t53 * t38 - t47 - t51) + g(2) * t46) * t39 + (g(1) * t46 + (-qJ(2) + t51 + (-pkin(3) - t53) * t38) * g(2)) * t41 + t42);
U = t1;
