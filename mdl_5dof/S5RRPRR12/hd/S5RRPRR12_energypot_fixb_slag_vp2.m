% Calculate potential energy for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:28:56
% EndTime: 2019-12-31 20:28:57
% DurationCPUTime: 0.40s
% Computational Cost: add. (120->55), mult. (222->56), div. (0->0), fcn. (216->8), ass. (0->27)
t108 = -mrSges(3,1) - mrSges(4,1);
t107 = mrSges(3,2) - mrSges(4,3);
t76 = sin(qJ(5));
t80 = cos(qJ(5));
t101 = -m(6) * pkin(4) - t80 * mrSges(6,1) + t76 * mrSges(6,2) - mrSges(5,1);
t77 = sin(qJ(4));
t78 = sin(qJ(2));
t81 = cos(qJ(4));
t82 = cos(qJ(2));
t60 = t78 * t77 + t82 * t81;
t106 = t60 * t101 + t107 * t78 + t108 * t82 - mrSges(2,1);
t105 = -m(5) - m(6);
t98 = t78 * t81;
t104 = t82 * t77 - t98;
t102 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t100 = t76 * mrSges(6,1) + t80 * mrSges(6,2) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3);
t79 = sin(qJ(1));
t97 = t79 * t82;
t83 = cos(qJ(1));
t95 = t83 * t82;
t94 = t83 * pkin(1) + t79 * pkin(6);
t93 = qJ(3) * t78;
t92 = t79 * pkin(1) - t83 * pkin(6);
t91 = pkin(2) * t95 + t83 * t93 + t94;
t90 = t78 * pkin(2) - t82 * qJ(3) + pkin(5);
t87 = pkin(2) * t97 + t79 * t93 + t92;
t1 = (-m(4) * t90 - mrSges(1,3) - mrSges(2,3) + t105 * (t78 * pkin(3) + t90) - t107 * t82 + t108 * t78 - t101 * t104 + t102 * t60 + (-m(2) - m(3)) * pkin(5)) * g(3) + (-m(3) * t92 - m(4) * t87 - mrSges(1,2) + t105 * (pkin(3) * t97 + t83 * pkin(7) + t87) - t100 * t83) * g(2) + (-m(3) * t94 - m(4) * t91 - mrSges(1,1) + t105 * (pkin(3) * t95 + t91) + t102 * t77 * t95 + (-t102 * t98 + t106) * t83) * g(1) + ((t102 * t104 + t106) * g(2) + (-t105 * pkin(7) + t100) * g(1)) * t79;
U = t1;
