% Calculate potential energy for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:39
% EndTime: 2019-12-31 19:05:39
% DurationCPUTime: 0.24s
% Computational Cost: add. (106->48), mult. (148->41), div. (0->0), fcn. (141->8), ass. (0->19)
t68 = -mrSges(2,1) - mrSges(3,1);
t67 = mrSges(2,2) - mrSges(3,3);
t66 = -m(4) - m(5) - m(6);
t49 = qJ(4) + qJ(5);
t42 = sin(t49);
t43 = cos(t49);
t50 = sin(qJ(4));
t51 = cos(qJ(4));
t65 = mrSges(4,1) + m(5) * pkin(3) + t51 * mrSges(5,1) - t50 * mrSges(5,2) + m(6) * (t51 * pkin(4) + pkin(3)) + t43 * mrSges(6,1) - t42 * mrSges(6,2);
t64 = m(5) * pkin(7) - m(6) * (-pkin(8) - pkin(7)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t63 = cos(qJ(3));
t62 = sin(qJ(1));
t61 = sin(qJ(3));
t52 = cos(qJ(1));
t60 = t52 * pkin(1) + t62 * qJ(2);
t58 = t62 * pkin(1) - t52 * qJ(2);
t37 = t52 * t61 - t62 * t63;
t36 = -t52 * t63 - t62 * t61;
t1 = (t42 * mrSges(6,1) + t51 * mrSges(5,2) + t43 * mrSges(6,2) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + (m(6) * pkin(4) + mrSges(5,1)) * t50 + t66 * (-pkin(6) + pkin(5)) + (-m(2) - m(3)) * pkin(5)) * g(3) + (-m(3) * t58 - mrSges(1,2) + t68 * t62 - t67 * t52 + t66 * (t62 * pkin(2) + t58) + t65 * t37 + t64 * t36) * g(2) + (-m(3) * t60 - mrSges(1,1) + t67 * t62 + t68 * t52 + t66 * (t52 * pkin(2) + t60) - t64 * t37 + t65 * t36) * g(1);
U = t1;
