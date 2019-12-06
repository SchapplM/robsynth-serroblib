% Calculate potential energy for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(1,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_energypot_fixb_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:33
% EndTime: 2019-12-05 18:08:33
% DurationCPUTime: 0.16s
% Computational Cost: add. (62->55), mult. (123->88), div. (0->0), fcn. (115->8), ass. (0->21)
t32 = sin(qJ(3));
t33 = sin(qJ(1));
t44 = t32 * t33;
t35 = cos(qJ(4));
t43 = t32 * t35;
t37 = cos(qJ(1));
t42 = t32 * t37;
t36 = cos(qJ(3));
t41 = t33 * t36;
t40 = t36 * t37;
t39 = t37 * qJ(2);
t38 = rSges(4,1) * t36 - rSges(4,2) * t32;
t34 = cos(qJ(5));
t31 = sin(qJ(4));
t30 = sin(qJ(5));
t29 = t33 * qJ(2);
t28 = t33 * t31 + t35 * t40;
t27 = t31 * t40 - t33 * t35;
t26 = -t31 * t37 + t35 * t41;
t25 = t31 * t41 + t35 * t37;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t37 - t33 * rSges(2,2)) + g(2) * (t33 * rSges(2,1) + rSges(2,2) * t37) + g(3) * rSges(2,3)) - m(3) * (g(1) * (rSges(3,1) * t37 + t33 * rSges(3,3) + t29) + g(2) * (t33 * rSges(3,1) + (-rSges(3,3) - qJ(2)) * t37) + g(3) * rSges(3,2)) - m(4) * (g(1) * t29 + g(3) * (rSges(4,1) * t32 + rSges(4,2) * t36) + (g(1) * rSges(4,3) + g(2) * t38) * t33 + (g(1) * t38 + g(2) * (-rSges(4,3) - qJ(2))) * t37) - m(5) * (g(1) * (t28 * rSges(5,1) - t27 * rSges(5,2) + rSges(5,3) * t42 + t29) + g(2) * (t26 * rSges(5,1) - t25 * rSges(5,2) + rSges(5,3) * t44 - t39) + g(3) * (-t36 * rSges(5,3) + (rSges(5,1) * t35 - rSges(5,2) * t31) * t32)) - m(6) * (g(1) * (t29 + (t28 * t34 + t30 * t42) * rSges(6,1) + (-t28 * t30 + t34 * t42) * rSges(6,2) + t27 * rSges(6,3)) + g(2) * (-t39 + (t26 * t34 + t30 * t44) * rSges(6,1) + (-t26 * t30 + t34 * t44) * rSges(6,2) + t25 * rSges(6,3)) + g(3) * ((-t30 * t36 + t34 * t43) * rSges(6,1) + (-t30 * t43 - t34 * t36) * rSges(6,2) + t32 * t31 * rSges(6,3)));
U = t1;
