% Calculate potential energy for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:23
% EndTime: 2019-03-09 04:31:24
% DurationCPUTime: 0.44s
% Computational Cost: add. (232->68), mult. (243->66), div. (0->0), fcn. (228->8), ass. (0->33)
t109 = -mrSges(4,2) + mrSges(5,3) + mrSges(6,2) - mrSges(7,3);
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t108 = pkin(3) * t83 + pkin(8) * t80;
t107 = -m(6) - m(7);
t106 = mrSges(3,2) - mrSges(4,3);
t79 = sin(qJ(4));
t82 = cos(qJ(4));
t105 = (pkin(4) * t82 + qJ(5) * t79) * t80;
t104 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t103 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1) - mrSges(7,1);
t102 = -t83 * mrSges(4,1) - mrSges(3,1) + (m(7) * qJ(6) - t109) * t80;
t81 = sin(qJ(1));
t75 = t81 * pkin(1);
t84 = cos(qJ(1));
t76 = t84 * pkin(1);
t99 = t79 * t83;
t95 = t82 * t83;
t78 = qJ(2) + pkin(6);
t93 = pkin(3) * t80 + t78;
t77 = qJ(1) + pkin(9);
t72 = sin(t77);
t73 = cos(t77);
t92 = pkin(2) * t73 + pkin(7) * t72 + t76;
t91 = pkin(2) * t72 - pkin(7) * t73 + t75;
t90 = t108 * t73 + t92;
t88 = -pkin(8) * t83 + t93;
t87 = t108 * t72 + t91;
t61 = t72 * t79 + t73 * t95;
t60 = -t72 * t82 + t73 * t99;
t59 = t72 * t95 - t73 * t79;
t58 = t72 * t99 + t73 * t82;
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - mrSges(3,3) - m(5) * t88 - m(6) * (t88 + t105) - m(7) * (t93 + t105) + (-m(3) - m(4)) * t78 + (-m(7) * (-pkin(8) + qJ(6)) + t109) * t83 + (t103 * t82 + t104 * t79 - mrSges(4,1)) * t80) * g(3) + (-m(3) * t75 - m(4) * t91 - m(5) * t87 - t81 * mrSges(2,1) - mrSges(2,2) * t84 - mrSges(1,2) + t107 * (pkin(4) * t59 + qJ(5) * t58 + t87) - t106 * t73 + t103 * t59 + t104 * t58 + t102 * t72) * g(2) + (-m(3) * t76 - m(4) * t92 - m(5) * t90 - t84 * mrSges(2,1) + t81 * mrSges(2,2) - mrSges(1,1) + t107 * (pkin(4) * t61 + qJ(5) * t60 + t90) + t106 * t72 + t103 * t61 + t104 * t60 + t102 * t73) * g(1);
U  = t1;
