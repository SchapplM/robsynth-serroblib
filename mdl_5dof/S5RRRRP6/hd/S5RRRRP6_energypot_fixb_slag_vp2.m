% Calculate potential energy for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP6_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:45
% EndTime: 2019-12-31 21:52:45
% DurationCPUTime: 0.33s
% Computational Cost: add. (141->53), mult. (161->49), div. (0->0), fcn. (138->8), ass. (0->24)
t61 = cos(qJ(4));
t84 = -m(5) * pkin(3) - m(6) * (pkin(4) * t61 + pkin(3)) - mrSges(4,1);
t83 = m(5) * pkin(8) - m(6) * (-qJ(5) - pkin(8)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t82 = m(6) * pkin(4);
t81 = -mrSges(5,1) - mrSges(6,1);
t80 = mrSges(5,2) + mrSges(6,2);
t79 = -m(4) - m(5) - m(6);
t78 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t56 = qJ(2) + qJ(3);
t53 = sin(t56);
t54 = cos(t56);
t59 = sin(qJ(2));
t62 = cos(qJ(2));
t77 = -m(3) * pkin(1) - t62 * mrSges(3,1) + t59 * mrSges(3,2) - t83 * t53 + t84 * t54 - mrSges(2,1);
t58 = sin(qJ(4));
t60 = sin(qJ(1));
t75 = t58 * t60;
t63 = cos(qJ(1));
t74 = t58 * t63;
t73 = t60 * t61;
t72 = t61 * t63;
t64 = -pkin(7) - pkin(6);
t51 = pkin(2) * t62 + pkin(1);
t1 = (-t59 * mrSges(3,1) - t62 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(5) + t79 * (t59 * pkin(2) + pkin(5)) + t83 * t54 + (t80 * t58 + t81 * t61 + t84) * t53) * g(3) + (t74 * t82 - mrSges(1,2) + t79 * (t60 * t51 + t63 * t64) + t81 * (t54 * t73 - t74) - t80 * (-t54 * t75 - t72) - t78 * t63 + t77 * t60) * g(2) + (-t75 * t82 - mrSges(1,1) + t79 * (t63 * t51 - t60 * t64) + t81 * (t54 * t72 + t75) - t80 * (-t54 * t74 + t73) + t78 * t60 + t77 * t63) * g(1);
U = t1;
