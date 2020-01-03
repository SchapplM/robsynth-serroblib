% Calculate potential energy for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:25
% EndTime: 2020-01-03 12:05:26
% DurationCPUTime: 0.36s
% Computational Cost: add. (157->47), mult. (148->38), div. (0->0), fcn. (125->10), ass. (0->20)
t44 = qJ(4) + qJ(5);
t39 = sin(t44);
t41 = cos(t44);
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t69 = -mrSges(4,1) - m(5) * pkin(3) - t50 * mrSges(5,1) + t48 * mrSges(5,2) - m(6) * (pkin(4) * t50 + pkin(3)) - t41 * mrSges(6,1) + t39 * mrSges(6,2);
t68 = mrSges(4,2) - m(5) * pkin(7) + m(6) * (-pkin(8) - pkin(7)) - mrSges(5,3) - mrSges(6,3);
t64 = m(4) + m(5) + m(6);
t46 = sin(pkin(9));
t47 = cos(pkin(9));
t67 = t68 * t46 + t69 * t47 - mrSges(3,1);
t63 = m(3) + t64;
t62 = t39 * mrSges(6,1) + t50 * mrSges(5,2) + t41 * mrSges(6,2) - mrSges(3,2) + mrSges(4,3) + (m(6) * pkin(4) + mrSges(5,1)) * t48 + t64 * qJ(3);
t49 = sin(qJ(1));
t43 = t49 * pkin(1);
t51 = cos(qJ(1));
t45 = qJ(1) + qJ(2);
t42 = cos(t45);
t40 = sin(t45);
t1 = (-mrSges(2,2) * t49 - mrSges(1,3) + (t64 * pkin(2) - t67) * t42 + t62 * t40 + (t63 * pkin(1) + mrSges(2,1)) * t51) * g(3) + (-m(3) * t43 - mrSges(2,1) * t49 - mrSges(2,2) * t51 - mrSges(1,2) - t64 * (t40 * pkin(2) + t43) + t62 * t42 + t67 * t40) * g(2) + (-m(2) * pkin(5) - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) - t63 * (pkin(6) + pkin(5)) - t68 * t47 + t69 * t46) * g(1);
U = t1;
