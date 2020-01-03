% Calculate potential energy for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRP7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_energypot_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP7_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:00
% EndTime: 2019-12-31 16:47:00
% DurationCPUTime: 0.16s
% Computational Cost: add. (56->45), mult. (78->52), div. (0->0), fcn. (58->4), ass. (0->14)
t40 = rSges(5,1) + pkin(3);
t27 = sin(qJ(3));
t29 = cos(qJ(3));
t39 = rSges(4,1) * t27 + rSges(4,2) * t29;
t38 = rSges(5,3) + qJ(4);
t37 = pkin(2) + pkin(4);
t28 = sin(qJ(1));
t30 = cos(qJ(1));
t34 = t30 * pkin(1) + t28 * qJ(2);
t24 = t28 * pkin(1);
t33 = t28 * pkin(5) + t24;
t32 = t30 * pkin(5) + t34;
t31 = t40 * t27 - t38 * t29;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t30 - t28 * rSges(2,2)) + g(2) * (t28 * rSges(2,1) + rSges(2,2) * t30) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t30 + t28 * rSges(3,3) + t34) + g(2) * (-t28 * rSges(3,2) + t24 + (-rSges(3,3) - qJ(2)) * t30) + g(3) * (pkin(4) + rSges(3,1))) - m(4) * (g(1) * (t39 * t28 + t32) + g(2) * (t28 * rSges(4,3) + t33) + g(3) * (rSges(4,1) * t29 - rSges(4,2) * t27 + t37) + (g(1) * rSges(4,3) + g(2) * (-qJ(2) - t39)) * t30) - m(5) * (g(1) * t32 + g(2) * t33 + g(3) * (t38 * t27 + t40 * t29 + t37) + (g(2) * rSges(5,2) + g(1) * t31) * t28 + (g(1) * rSges(5,2) + g(2) * (-qJ(2) - t31)) * t30);
U = t1;
