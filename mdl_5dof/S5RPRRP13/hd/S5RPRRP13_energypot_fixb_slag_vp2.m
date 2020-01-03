% Calculate potential energy for
% S5RPRRP13
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP13_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:27
% EndTime: 2019-12-31 18:58:28
% DurationCPUTime: 0.35s
% Computational Cost: add. (98->53), mult. (163->48), div. (0->0), fcn. (144->6), ass. (0->23)
t88 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t87 = -m(5) - m(6);
t86 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t85 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t84 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t83 = -t61 * mrSges(4,1) - t88 * t64 + mrSges(2,2) - mrSges(3,3);
t82 = pkin(2) + pkin(5);
t81 = pkin(3) * t61;
t80 = pkin(7) * t64;
t60 = sin(qJ(4));
t65 = cos(qJ(1));
t79 = t60 * t65;
t62 = sin(qJ(1));
t78 = t62 * t60;
t63 = cos(qJ(4));
t77 = t62 * t63;
t74 = t65 * t63;
t56 = t62 * pkin(1);
t73 = t62 * pkin(6) + t56;
t72 = t65 * pkin(1) + t62 * qJ(2);
t1 = (-m(4) * t82 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t87 * (t64 * pkin(3) + t61 * pkin(7) + t82) + (-m(2) - m(3)) * pkin(5) + (-t84 * t60 + t85 * t63 - mrSges(4,1)) * t64 + t88 * t61) * g(3) + (-m(3) * t56 - m(4) * t73 - mrSges(1,2) + t87 * (t65 * t80 + t73) + t85 * (-t61 * t74 + t78) + t84 * (t61 * t79 + t77) + t86 * t62 + (t87 * (-qJ(2) - t81) + (m(3) + m(4)) * qJ(2) - t83) * t65) * g(2) + (-m(3) * t72 - mrSges(1,1) + t85 * (t61 * t77 + t79) - t84 * (t61 * t78 - t74) + t86 * t65 + (-m(4) + t87) * (t65 * pkin(6) + t72) + (t87 * (-t80 + t81) + t83) * t62) * g(1);
U = t1;
