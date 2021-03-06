% Calculate potential energy for
% S4PRRR3
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR3_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:32
% EndTime: 2019-12-31 16:31:32
% DurationCPUTime: 0.08s
% Computational Cost: add. (84->42), mult. (58->46), div. (0->0), fcn. (38->8), ass. (0->19)
t43 = rSges(5,3) + pkin(6);
t42 = pkin(4) + qJ(1);
t33 = pkin(7) + qJ(2);
t28 = sin(t33);
t34 = sin(pkin(7));
t31 = t34 * pkin(1);
t41 = pkin(2) * t28 + t31;
t29 = cos(t33);
t35 = cos(pkin(7));
t32 = t35 * pkin(1);
t40 = pkin(2) * t29 + t32;
t39 = pkin(5) + t42;
t36 = sin(qJ(4));
t37 = cos(qJ(4));
t38 = rSges(5,1) * t37 - rSges(5,2) * t36 + pkin(3);
t30 = qJ(3) + t33;
t27 = cos(t30);
t26 = sin(t30);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t35 * rSges(2,1) - t34 * rSges(2,2)) + g(2) * (t34 * rSges(2,1) + t35 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t29 * rSges(3,1) - t28 * rSges(3,2) + t32) + g(2) * (t28 * rSges(3,1) + t29 * rSges(3,2) + t31) + g(3) * (rSges(3,3) + t42)) - m(4) * (g(1) * (t27 * rSges(4,1) - t26 * rSges(4,2) + t40) + g(2) * (t26 * rSges(4,1) + t27 * rSges(4,2) + t41) + g(3) * (rSges(4,3) + t39)) - m(5) * (g(1) * t40 + g(2) * t41 + g(3) * (t36 * rSges(5,1) + t37 * rSges(5,2) + t39) + (g(1) * t38 - g(2) * t43) * t27 + (g(1) * t43 + g(2) * t38) * t26);
U = t1;
