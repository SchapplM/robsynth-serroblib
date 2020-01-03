% Calculate potential energy for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR11_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:29
% EndTime: 2019-12-31 18:05:29
% DurationCPUTime: 0.25s
% Computational Cost: add. (80->66), mult. (109->78), div. (0->0), fcn. (89->6), ass. (0->20)
t52 = rSges(6,3) + pkin(7);
t51 = pkin(2) + pkin(5);
t35 = sin(qJ(4));
t36 = sin(qJ(1));
t49 = t36 * t35;
t37 = cos(qJ(5));
t48 = t36 * t37;
t39 = cos(qJ(1));
t47 = t39 * t35;
t46 = t39 * t37;
t31 = t36 * pkin(1);
t45 = t36 * qJ(3) + t31;
t44 = t39 * pkin(1) + t36 * qJ(2);
t43 = pkin(3) + t51;
t42 = t39 * pkin(6) + t45;
t41 = t39 * qJ(3) + t44;
t38 = cos(qJ(4));
t40 = rSges(5,1) * t35 + rSges(5,2) * t38;
t34 = sin(qJ(5));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t39 * rSges(2,1) - t36 * rSges(2,2)) + g(2) * (t36 * rSges(2,1) + t39 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (-t39 * rSges(3,2) + t36 * rSges(3,3) + t44) + g(2) * (-t36 * rSges(3,2) + t31 + (-rSges(3,3) - qJ(2)) * t39) + g(3) * (pkin(5) + rSges(3,1))) - m(4) * (g(1) * (t36 * rSges(4,2) + t39 * rSges(4,3) + t41) + g(2) * (t36 * rSges(4,3) + (-rSges(4,2) - qJ(2)) * t39 + t45) + g(3) * (rSges(4,1) + t51)) - m(5) * (g(1) * t41 + g(2) * t42 + g(3) * (t38 * rSges(5,1) - t35 * rSges(5,2) + t43) + (g(1) * t40 + g(2) * (rSges(5,3) - qJ(2))) * t39 + (g(1) * (-rSges(5,3) - pkin(6)) + g(2) * t40) * t36) - m(6) * (g(1) * (pkin(4) * t47 - t36 * pkin(6) + (-t36 * t34 + t35 * t46) * rSges(6,1) + (-t34 * t47 - t48) * rSges(6,2) + t41) + g(2) * (pkin(4) * t49 - t39 * qJ(2) + (t39 * t34 + t35 * t48) * rSges(6,1) + (-t34 * t49 + t46) * rSges(6,2) + t42) + g(3) * (t52 * t35 + t43) + (g(3) * (rSges(6,1) * t37 - rSges(6,2) * t34 + pkin(4)) - (g(1) * t39 + g(2) * t36) * t52) * t38);
U = t1;
