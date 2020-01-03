% Calculate potential energy for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR6_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:27
% EndTime: 2019-12-31 16:52:27
% DurationCPUTime: 0.11s
% Computational Cost: add. (83->45), mult. (81->52), div. (0->0), fcn. (61->8), ass. (0->20)
t37 = sin(pkin(7));
t48 = t37 * pkin(2) + pkin(4);
t38 = cos(pkin(7));
t29 = t38 * pkin(2) + pkin(1);
t39 = -pkin(5) - qJ(2);
t47 = rSges(4,3) - t39;
t46 = rSges(5,3) + pkin(6) - t39;
t45 = rSges(3,3) + qJ(2);
t36 = pkin(7) + qJ(3);
t44 = rSges(3,1) * t38 - rSges(3,2) * t37 + pkin(1);
t30 = sin(t36);
t31 = cos(t36);
t43 = rSges(4,1) * t31 - rSges(4,2) * t30 + t29;
t32 = qJ(4) + t36;
t27 = sin(t32);
t28 = cos(t32);
t42 = rSges(5,1) * t28 - rSges(5,2) * t27 + pkin(3) * t31 + t29;
t41 = cos(qJ(1));
t40 = sin(qJ(1));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t41 * rSges(2,1) - t40 * rSges(2,2)) + g(2) * (t40 * rSges(2,1) + t41 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(3) * (t37 * rSges(3,1) + t38 * rSges(3,2) + pkin(4)) + (g(1) * t44 - g(2) * t45) * t41 + (g(1) * t45 + g(2) * t44) * t40) - m(4) * (g(3) * (t30 * rSges(4,1) + t31 * rSges(4,2) + t48) + (g(1) * t43 - g(2) * t47) * t41 + (g(1) * t47 + g(2) * t43) * t40) - m(5) * (g(3) * (t27 * rSges(5,1) + t28 * rSges(5,2) + pkin(3) * t30 + t48) + (g(1) * t42 - g(2) * t46) * t41 + (g(1) * t46 + g(2) * t42) * t40);
U = t1;
