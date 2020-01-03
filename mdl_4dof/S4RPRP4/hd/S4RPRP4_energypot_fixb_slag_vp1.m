% Calculate potential energy for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP4_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:35
% EndTime: 2019-12-31 16:43:35
% DurationCPUTime: 0.15s
% Computational Cost: add. (83->45), mult. (76->52), div. (0->0), fcn. (56->6), ass. (0->17)
t45 = rSges(5,1) + pkin(3);
t33 = sin(qJ(3));
t35 = cos(qJ(3));
t44 = rSges(4,1) * t35 - rSges(4,2) * t33;
t43 = rSges(5,3) + qJ(4);
t40 = pkin(4) + qJ(2);
t32 = qJ(1) + pkin(6);
t28 = sin(t32);
t34 = sin(qJ(1));
t30 = t34 * pkin(1);
t39 = t28 * pkin(2) + t30;
t29 = cos(t32);
t36 = cos(qJ(1));
t31 = t36 * pkin(1);
t38 = t29 * pkin(2) + t28 * pkin(5) + t31;
t37 = t43 * t33 + t45 * t35;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t36 * rSges(2,1) - t34 * rSges(2,2)) + g(2) * (t34 * rSges(2,1) + t36 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t29 * rSges(3,1) - t28 * rSges(3,2) + t31) + g(2) * (t28 * rSges(3,1) + t29 * rSges(3,2) + t30) + g(3) * (rSges(3,3) + t40)) - m(4) * (g(1) * (t28 * rSges(4,3) + t38) + g(2) * (t44 * t28 + t39) + g(3) * (t33 * rSges(4,1) + t35 * rSges(4,2) + t40) + (g(1) * t44 + g(2) * (-rSges(4,3) - pkin(5))) * t29) - m(5) * (g(1) * t38 + g(2) * t39 + g(3) * (t45 * t33 - t43 * t35 + t40) + (g(1) * rSges(5,2) + g(2) * t37) * t28 + (g(1) * t37 + g(2) * (-rSges(5,2) - pkin(5))) * t29);
U = t1;
