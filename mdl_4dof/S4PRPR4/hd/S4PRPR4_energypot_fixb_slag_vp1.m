% Calculate potential energy for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRPR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR4_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:54
% EndTime: 2019-12-31 16:21:54
% DurationCPUTime: 0.08s
% Computational Cost: add. (73->43), mult. (60->47), div. (0->0), fcn. (40->6), ass. (0->15)
t37 = rSges(5,3) + pkin(5);
t36 = pkin(4) + qJ(1);
t28 = pkin(6) + qJ(2);
t24 = sin(t28);
t29 = sin(pkin(6));
t26 = t29 * pkin(1);
t35 = t24 * pkin(2) + t26;
t25 = cos(t28);
t30 = cos(pkin(6));
t27 = t30 * pkin(1);
t34 = t25 * pkin(2) + t24 * qJ(3) + t27;
t31 = sin(qJ(4));
t32 = cos(qJ(4));
t33 = rSges(5,1) * t31 + rSges(5,2) * t32;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t30 * rSges(2,1) - t29 * rSges(2,2)) + g(2) * (t29 * rSges(2,1) + t30 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t25 * rSges(3,1) - t24 * rSges(3,2) + t27) + g(2) * (t24 * rSges(3,1) + t25 * rSges(3,2) + t26) + g(3) * (rSges(3,3) + t36)) - m(4) * (g(1) * (-t25 * rSges(4,2) + t24 * rSges(4,3) + t34) + g(2) * (-t24 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t25 + t35) + g(3) * (rSges(4,1) + t36)) - m(5) * (g(1) * t34 + g(2) * t35 + g(3) * (t32 * rSges(5,1) - t31 * rSges(5,2) + pkin(3) + t36) + (g(1) * t33 + g(2) * t37) * t24 + (g(1) * t37 + g(2) * (-qJ(3) - t33)) * t25);
U = t1;
