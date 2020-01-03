% Calculate potential energy for
% S4RPRP6
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_energypot_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP6_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:57
% EndTime: 2019-12-31 16:45:57
% DurationCPUTime: 0.15s
% Computational Cost: add. (55->43), mult. (73->49), div. (0->0), fcn. (53->4), ass. (0->14)
t34 = rSges(5,1) + pkin(3);
t33 = pkin(2) + pkin(4);
t32 = rSges(4,3) + pkin(5);
t31 = rSges(5,3) + qJ(4) + pkin(5);
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t30 = t26 * pkin(1) + t24 * qJ(2);
t23 = sin(qJ(3));
t25 = cos(qJ(3));
t29 = rSges(4,1) * t23 + rSges(4,2) * t25;
t20 = t24 * pkin(1);
t28 = g(1) * t30 + g(2) * t20;
t27 = rSges(5,2) * t25 + t34 * t23;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t26 * rSges(2,1) - t24 * rSges(2,2)) + g(2) * (t24 * rSges(2,1) + t26 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (-t26 * rSges(3,2) + t24 * rSges(3,3) + t30) + g(2) * (-t24 * rSges(3,2) + t20 + (-rSges(3,3) - qJ(2)) * t26) + g(3) * (pkin(4) + rSges(3,1))) - m(4) * (g(3) * (t25 * rSges(4,1) - t23 * rSges(4,2) + t33) + (g(1) * t29 + g(2) * t32) * t24 + (g(1) * t32 + g(2) * (-qJ(2) - t29)) * t26 + t28) - m(5) * (g(3) * (-t23 * rSges(5,2) + t34 * t25 + t33) + (g(1) * t27 + g(2) * t31) * t24 + (g(1) * t31 + g(2) * (-qJ(2) - t27)) * t26 + t28);
U = t1;
