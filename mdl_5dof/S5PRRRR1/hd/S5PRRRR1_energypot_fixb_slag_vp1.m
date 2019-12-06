% Calculate potential energy for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_energypot_fixb_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:02:48
% EndTime: 2019-12-05 17:02:49
% DurationCPUTime: 0.14s
% Computational Cost: add. (69->51), mult. (87->67), div. (0->0), fcn. (67->8), ass. (0->21)
t34 = cos(qJ(3));
t46 = pkin(2) * t34;
t31 = sin(qJ(3));
t45 = t31 * pkin(2);
t35 = cos(qJ(2));
t44 = t35 * t46 + pkin(1);
t29 = qJ(3) + qJ(4);
t27 = sin(t29);
t43 = t27 * rSges(6,3);
t30 = sin(qJ(5));
t32 = sin(qJ(2));
t42 = t32 * t30;
t33 = cos(qJ(5));
t41 = t32 * t33;
t40 = t35 * t30;
t39 = t35 * t33;
t38 = t32 * t46 + qJ(1);
t37 = rSges(4,1) * t34 - rSges(4,2) * t31;
t28 = cos(t29);
t36 = rSges(5,1) * t28 - rSges(5,2) * t27;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * rSges(2,1) + g(2) * rSges(2,2) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t35 * rSges(3,1) - t32 * rSges(3,2) + pkin(1)) - g(2) * rSges(3,3) + g(3) * (t32 * rSges(3,1) + t35 * rSges(3,2) + qJ(1))) - m(4) * (g(1) * (t32 * rSges(4,3) + t37 * t35 + pkin(1)) + g(2) * (-t31 * rSges(4,1) - t34 * rSges(4,2)) + g(3) * (-t35 * rSges(4,3) + t37 * t32 + qJ(1))) - m(5) * (g(1) * (t32 * rSges(5,3) + t36 * t35 + t44) + g(2) * (-t27 * rSges(5,1) - t28 * rSges(5,2) - t45) + g(3) * (-t35 * rSges(5,3) + t36 * t32 + t38)) - m(6) * (g(1) * ((t28 * t39 + t42) * rSges(6,1) + (-t28 * t40 + t41) * rSges(6,2) + t35 * t43 + t44) + g(2) * (t28 * rSges(6,3) - t45 + (-rSges(6,1) * t33 + rSges(6,2) * t30) * t27) + g(3) * ((t28 * t41 - t40) * rSges(6,1) + (-t28 * t42 - t39) * rSges(6,2) + t32 * t43 + t38));
U = t1;
