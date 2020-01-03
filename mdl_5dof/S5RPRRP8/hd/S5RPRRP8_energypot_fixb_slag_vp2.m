% Calculate potential energy for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:04
% EndTime: 2019-12-31 18:47:04
% DurationCPUTime: 0.24s
% Computational Cost: add. (105->46), mult. (163->39), div. (0->0), fcn. (160->6), ass. (0->20)
t70 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t69 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t68 = -m(5) - m(6);
t67 = -mrSges(2,1) - mrSges(3,1);
t66 = mrSges(2,2) - mrSges(3,3);
t49 = sin(qJ(4));
t50 = cos(qJ(4));
t65 = t69 * t49 + t70 * t50 + mrSges(4,1);
t64 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t63 = cos(qJ(3));
t62 = sin(qJ(1));
t61 = sin(qJ(3));
t51 = cos(qJ(1));
t60 = t51 * pkin(1) + t62 * qJ(2);
t59 = t51 * pkin(2) + t60;
t58 = t62 * pkin(1) - t51 * qJ(2);
t57 = t62 * pkin(2) + t58;
t40 = t51 * t61 - t62 * t63;
t39 = -t51 * t63 - t61 * t62;
t1 = (-mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + (-m(4) + t68) * (-pkin(6) + pkin(5)) - t69 * t50 + t70 * t49 + (-m(2) - m(3)) * pkin(5)) * g(3) + (-m(3) * t58 - m(4) * t57 - mrSges(1,2) + t67 * t62 + t68 * (-t40 * pkin(3) - t39 * pkin(7) + t57) - t66 * t51 + t65 * t40 - t64 * t39) * g(2) + (-m(3) * t60 - m(4) * t59 - mrSges(1,1) + t66 * t62 + t68 * (-t39 * pkin(3) + t40 * pkin(7) + t59) + t67 * t51 + t64 * t40 + t65 * t39) * g(1);
U = t1;
