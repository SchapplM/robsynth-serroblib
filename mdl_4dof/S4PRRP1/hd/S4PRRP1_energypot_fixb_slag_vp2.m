% Calculate potential energy for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP1_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:22:47
% EndTime: 2019-03-08 18:22:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (77->35), mult. (55->27), div. (0->0), fcn. (30->6), ass. (0->14)
t41 = -m(4) - m(5);
t40 = pkin(4) + qJ(1);
t32 = pkin(6) + qJ(2);
t37 = -m(3) * pkin(1) - mrSges(2,1);
t36 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t35 = m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3);
t34 = cos(pkin(6));
t33 = sin(pkin(6));
t28 = qJ(3) + t32;
t27 = cos(t32);
t26 = sin(t32);
t25 = cos(t28);
t24 = sin(t28);
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t40 - mrSges(3,3) - mrSges(4,3) - mrSges(5,2) + t41 * (pkin(5) + t40)) * g(3) + (-t26 * mrSges(3,1) - t34 * mrSges(2,2) - t27 * mrSges(3,2) + t36 * t24 + t35 * t25 + t37 * t33 - mrSges(1,2) + t41 * (t33 * pkin(1) + pkin(2) * t26)) * g(2) + (-t27 * mrSges(3,1) + t33 * mrSges(2,2) + t26 * mrSges(3,2) - t35 * t24 + t36 * t25 + t37 * t34 - mrSges(1,1) + t41 * (t34 * pkin(1) + pkin(2) * t27)) * g(1);
U  = t1;
