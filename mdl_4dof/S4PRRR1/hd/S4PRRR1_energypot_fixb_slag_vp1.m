% Calculate potential energy for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR1_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:03
% EndTime: 2019-03-08 18:25:03
% DurationCPUTime: 0.07s
% Computational Cost: add. (76->42), mult. (48->42), div. (0->0), fcn. (28->8), ass. (0->18)
t41 = pkin(4) + qJ(1);
t35 = pkin(7) + qJ(2);
t30 = sin(t35);
t36 = sin(pkin(7));
t33 = t36 * pkin(1);
t40 = pkin(2) * t30 + t33;
t31 = cos(t35);
t37 = cos(pkin(7));
t34 = t37 * pkin(1);
t39 = pkin(2) * t31 + t34;
t38 = pkin(5) + t41;
t32 = qJ(3) + t35;
t29 = qJ(4) + t32;
t28 = cos(t32);
t27 = sin(t32);
t24 = cos(t29);
t23 = sin(t29);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t37 * rSges(2,1) - t36 * rSges(2,2)) + g(2) * (t36 * rSges(2,1) + t37 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t31 * rSges(3,1) - t30 * rSges(3,2) + t34) + g(2) * (t30 * rSges(3,1) + t31 * rSges(3,2) + t33) + g(3) * (rSges(3,3) + t41)) - m(4) * (g(1) * (t28 * rSges(4,1) - t27 * rSges(4,2) + t39) + g(2) * (t27 * rSges(4,1) + t28 * rSges(4,2) + t40) + g(3) * (rSges(4,3) + t38)) - m(5) * (g(1) * (t24 * rSges(5,1) - t23 * rSges(5,2) + pkin(3) * t28 + t39) + g(2) * (t23 * rSges(5,1) + t24 * rSges(5,2) + pkin(3) * t27 + t40) + g(3) * (pkin(6) + rSges(5,3) + t38));
U  = t1;
