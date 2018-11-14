% Calculate potential energy for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:47
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4RPPR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR1_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR1_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:46:24
% EndTime: 2018-11-14 13:46:24
% DurationCPUTime: 0.08s
% Computational Cost: add. (76->44), mult. (60->48), div. (0->0), fcn. (44->6), ass. (0->15)
t37 = pkin(4) + qJ(2);
t30 = qJ(1) + pkin(6);
t26 = sin(t30);
t32 = sin(qJ(1));
t28 = t32 * pkin(1);
t36 = t26 * pkin(2) + t28;
t27 = cos(t30);
t34 = cos(qJ(1));
t29 = t34 * pkin(1);
t35 = t27 * pkin(2) + t26 * qJ(3) + t29;
t33 = cos(qJ(4));
t31 = sin(qJ(4));
t22 = t26 * t33 - t27 * t31;
t21 = -t26 * t31 - t27 * t33;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t34 * rSges(2,1) - t32 * rSges(2,2)) + g(2) * (t32 * rSges(2,1) + t34 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t27 * rSges(3,1) - t26 * rSges(3,2) + t29) + g(2) * (t26 * rSges(3,1) + t27 * rSges(3,2) + t28) + g(3) * (rSges(3,3) + t37)) - m(4) * (g(1) * (t27 * rSges(4,1) + t26 * rSges(4,3) + t35) + g(2) * (t26 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t27 + t36) + g(3) * (rSges(4,2) + t37)) - m(5) * (g(1) * (-t21 * rSges(5,1) + t22 * rSges(5,2) + t27 * pkin(3) + t35) + g(2) * (t22 * rSges(5,1) + t21 * rSges(5,2) + t26 * pkin(3) - t27 * qJ(3) + t36) + g(3) * (-pkin(5) - rSges(5,3) + t37));
U  = t1;
