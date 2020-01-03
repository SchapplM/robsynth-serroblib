% Calculate potential energy for
% S5RRRRP8
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:34
% EndTime: 2019-12-31 21:59:34
% DurationCPUTime: 0.37s
% Computational Cost: add. (142->55), mult. (187->52), div. (0->0), fcn. (168->8), ass. (0->24)
t64 = cos(qJ(3));
t52 = t64 * pkin(3) + pkin(2);
t60 = qJ(3) + qJ(4);
t54 = cos(t60);
t61 = sin(qJ(3));
t87 = -m(5) * t52 - m(6) * (pkin(4) * t54 + t52) - mrSges(3,1) - m(4) * pkin(2) - t64 * mrSges(4,1) + t61 * mrSges(4,2);
t67 = -pkin(8) - pkin(7);
t86 = m(5) * t67 + m(6) * (-qJ(5) + t67) + mrSges(3,2) - mrSges(5,3) - mrSges(6,3) - m(4) * pkin(7) - mrSges(4,3);
t85 = -mrSges(5,1) - mrSges(6,1);
t84 = mrSges(5,2) + mrSges(6,2);
t83 = -m(3) - m(5) - m(6);
t82 = -m(4) + t83;
t62 = sin(qJ(2));
t65 = cos(qJ(2));
t80 = t86 * t62 + t87 * t65 - mrSges(2,1);
t53 = sin(t60);
t78 = pkin(3) * t61;
t79 = -m(5) * t78 - m(6) * (pkin(4) * t53 + t78) + mrSges(2,2) - mrSges(3,3) - t61 * mrSges(4,1) - t64 * mrSges(4,2);
t63 = sin(qJ(1));
t77 = t63 * t65;
t66 = cos(qJ(1));
t76 = t65 * t66;
t56 = t63 * pkin(1);
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) + t82) * pkin(5) - t86 * t65 + (t84 * t53 + t85 * t54 + t87) * t62) * g(3) + (-m(4) * t56 - mrSges(1,2) + t85 * (-t53 * t66 + t54 * t77) - t84 * (-t53 * t77 - t54 * t66) + t83 * (-t66 * pkin(6) + t56) + (m(4) * pkin(6) - t79) * t66 + t80 * t63) * g(2) + (-mrSges(1,1) + t85 * (t53 * t63 + t54 * t76) - t84 * (-t53 * t76 + t54 * t63) + t82 * (t66 * pkin(1) + t63 * pkin(6)) + t79 * t63 + t80 * t66) * g(1);
U = t1;
