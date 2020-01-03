% Calculate potential energy for
% S4RPRR3
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR3_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:03
% EndTime: 2019-12-31 16:49:03
% DurationCPUTime: 0.13s
% Computational Cost: add. (84->44), mult. (69->50), div. (0->0), fcn. (49->8), ass. (0->19)
t43 = rSges(4,3) + pkin(5);
t42 = rSges(5,3) + pkin(6) + pkin(5);
t41 = pkin(4) + qJ(2);
t34 = sin(qJ(1));
t29 = t34 * pkin(1);
t36 = cos(qJ(1));
t30 = t36 * pkin(1);
t40 = g(1) * t30 + g(2) * t29;
t33 = sin(qJ(3));
t35 = cos(qJ(3));
t39 = rSges(4,1) * t35 - rSges(4,2) * t33 + pkin(2);
t32 = qJ(3) + qJ(4);
t27 = sin(t32);
t28 = cos(t32);
t38 = rSges(5,1) * t28 - rSges(5,2) * t27 + t35 * pkin(3) + pkin(2);
t31 = qJ(1) + pkin(7);
t26 = cos(t31);
t25 = sin(t31);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t36 * rSges(2,1) - t34 * rSges(2,2)) + g(2) * (t34 * rSges(2,1) + t36 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t26 * rSges(3,1) - t25 * rSges(3,2) + t30) + g(2) * (t25 * rSges(3,1) + t26 * rSges(3,2) + t29) + g(3) * (rSges(3,3) + t41)) - m(4) * (g(3) * (t33 * rSges(4,1) + t35 * rSges(4,2) + t41) + (g(1) * t39 - g(2) * t43) * t26 + (g(1) * t43 + g(2) * t39) * t25 + t40) - m(5) * (g(3) * (t27 * rSges(5,1) + t28 * rSges(5,2) + t33 * pkin(3) + t41) + (g(1) * t38 - g(2) * t42) * t26 + (g(1) * t42 + g(2) * t38) * t25 + t40);
U = t1;
