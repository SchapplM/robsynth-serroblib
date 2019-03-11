% Calculate potential energy for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP1_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:26
% EndTime: 2019-03-08 18:29:26
% DurationCPUTime: 0.07s
% Computational Cost: add. (77->40), mult. (50->42), div. (0->0), fcn. (30->6), ass. (0->17)
t39 = rSges(5,1) + pkin(3);
t38 = pkin(4) + qJ(2);
t31 = qJ(1) + pkin(6);
t26 = sin(t31);
t32 = sin(qJ(1));
t29 = t32 * pkin(1);
t37 = pkin(2) * t26 + t29;
t27 = cos(t31);
t33 = cos(qJ(1));
t30 = t33 * pkin(1);
t36 = pkin(2) * t27 + t30;
t35 = rSges(5,3) + qJ(4);
t34 = pkin(5) + t38;
t28 = qJ(3) + t31;
t25 = cos(t28);
t24 = sin(t28);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t33 * rSges(2,1) - t32 * rSges(2,2)) + g(2) * (t32 * rSges(2,1) + t33 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t27 * rSges(3,1) - t26 * rSges(3,2) + t30) + g(2) * (t26 * rSges(3,1) + t27 * rSges(3,2) + t29) + g(3) * (rSges(3,3) + t38)) - m(4) * (g(1) * (t25 * rSges(4,1) - t24 * rSges(4,2) + t36) + g(2) * (t24 * rSges(4,1) + t25 * rSges(4,2) + t37) + g(3) * (rSges(4,3) + t34)) - m(5) * (g(1) * t36 + g(2) * t37 + g(3) * (rSges(5,2) + t34) + (g(1) * t39 - g(2) * t35) * t25 + (g(1) * t35 + g(2) * t39) * t24);
U  = t1;
