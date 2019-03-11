% Calculate potential energy for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PPRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR1_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:27
% EndTime: 2019-03-08 18:15:27
% DurationCPUTime: 0.08s
% Computational Cost: add. (65->50), mult. (74->56), div. (0->0), fcn. (62->6), ass. (0->17)
t32 = sin(pkin(6));
t34 = sin(qJ(3));
t38 = t32 * t34;
t37 = qJ(1) - pkin(4);
t33 = cos(pkin(6));
t36 = t33 * pkin(1) + t32 * qJ(2);
t35 = cos(qJ(3));
t31 = qJ(3) + qJ(4);
t29 = t32 * pkin(1);
t28 = cos(t31);
t27 = sin(t31);
t25 = t35 * pkin(3) + pkin(2);
t24 = t32 * t35 - t33 * t34;
t23 = -t33 * t35 - t38;
t22 = -t33 * t27 + t32 * t28;
t21 = -t32 * t27 - t33 * t28;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t33 * rSges(2,1) - t32 * rSges(2,2)) + g(2) * (t32 * rSges(2,1) + t33 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t33 * rSges(3,1) + t32 * rSges(3,3) + t36) + g(2) * (t32 * rSges(3,1) + t29 + (-rSges(3,3) - qJ(2)) * t33) + g(3) * (qJ(1) + rSges(3,2))) - m(4) * (g(1) * (-t23 * rSges(4,1) + t24 * rSges(4,2) + t33 * pkin(2) + t36) + g(2) * (t24 * rSges(4,1) + t23 * rSges(4,2) + t32 * pkin(2) - t33 * qJ(2) + t29) + g(3) * (-rSges(4,3) + t37)) - m(5) * (g(1) * (-t21 * rSges(5,1) + t22 * rSges(5,2) + pkin(3) * t38 + t33 * t25 + t36) + g(2) * (t22 * rSges(5,1) + t21 * rSges(5,2) + t32 * t25 + t29 + (-pkin(3) * t34 - qJ(2)) * t33) + g(3) * (-pkin(5) - rSges(5,3) + t37));
U  = t1;
