% Calculate potential energy for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP7_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:30
% EndTime: 2019-12-31 18:44:30
% DurationCPUTime: 0.37s
% Computational Cost: add. (155->46), mult. (161->38), div. (0->0), fcn. (142->8), ass. (0->21)
t73 = sin(qJ(4));
t76 = cos(qJ(4));
t91 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t92 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t102 = t73 * t91 + t76 * t92 - mrSges(4,1);
t99 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t95 = -m(5) - m(6);
t98 = -m(4) + t95;
t97 = t73 * t92 - t76 * t91 + mrSges(3,2) - mrSges(4,3);
t74 = sin(qJ(3));
t77 = cos(qJ(3));
t96 = -mrSges(3,1) + (t95 * pkin(7) + t99) * t74 + (t95 * pkin(3) + t102) * t77;
t75 = sin(qJ(1));
t69 = t75 * pkin(1);
t78 = cos(qJ(1));
t70 = t78 * pkin(1);
t72 = qJ(2) + pkin(5);
t71 = qJ(1) + pkin(8);
t67 = cos(t71);
t66 = sin(t71);
t1 = (-m(2) * pkin(5) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t95 * (t74 * pkin(3) - pkin(7) * t77 + t72) + (-m(3) - m(4)) * t72 - t99 * t77 + t102 * t74) * g(3) + (-m(3) * t69 - t75 * mrSges(2,1) - t78 * mrSges(2,2) - mrSges(1,2) + t98 * (t66 * pkin(2) - pkin(6) * t67 + t69) - t97 * t67 + t96 * t66) * g(2) + (-m(3) * t70 - t78 * mrSges(2,1) + t75 * mrSges(2,2) - mrSges(1,1) + t98 * (t67 * pkin(2) + t66 * pkin(6) + t70) + t97 * t66 + t96 * t67) * g(1);
U = t1;
