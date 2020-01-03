% Calculate potential energy for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRP6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP6_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:57
% EndTime: 2019-12-31 16:45:57
% DurationCPUTime: 0.17s
% Computational Cost: add. (55->32), mult. (78->21), div. (0->0), fcn. (53->4), ass. (0->11)
t41 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t40 = mrSges(4,2) + mrSges(5,2);
t39 = -m(4) - m(5);
t38 = m(3) - t39;
t36 = -m(4) * pkin(5) + m(5) * (-qJ(4) - pkin(5)) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t27 = sin(qJ(3));
t29 = cos(qJ(3));
t35 = t41 * t27 - t40 * t29 + mrSges(2,2) - mrSges(3,3);
t30 = cos(qJ(1));
t28 = sin(qJ(1));
t1 = (-mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t39 * (pkin(2) + pkin(4)) + t41 * t29 + t40 * t27 + (-m(2) - m(3)) * pkin(4)) * g(3) + (-mrSges(1,2) + (t38 * qJ(2) - t35) * t30 + (-t38 * pkin(1) + t36) * t28) * g(2) + (-mrSges(1,1) - t38 * (t30 * pkin(1) + t28 * qJ(2)) + t36 * t30 + t35 * t28) * g(1);
U = t1;
