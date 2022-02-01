% Calculate potential energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:11
% EndTime: 2022-01-23 09:16:12
% DurationCPUTime: 0.52s
% Computational Cost: add. (152->70), mult. (177->63), div. (0->0), fcn. (158->10), ass. (0->30)
t77 = qJ(3) + pkin(6);
t104 = -mrSges(4,3) - mrSges(5,3) - m(4) * qJ(3) - m(5) * t77 + mrSges(3,2) + m(6) * (-pkin(7) - t77) - mrSges(6,3);
t72 = pkin(9) + qJ(4);
t66 = qJ(5) + t72;
t61 = sin(t66);
t62 = cos(t66);
t75 = cos(pkin(9));
t63 = t75 * pkin(3) + pkin(2);
t65 = cos(t72);
t103 = -m(4) * pkin(2) - m(5) * t63 - mrSges(3,1) - m(6) * (pkin(4) * t65 + t63) - mrSges(6,1) * t62 + mrSges(6,2) * t61;
t99 = -m(3) - m(6);
t101 = -m(4) + t99;
t74 = sin(pkin(8));
t76 = cos(pkin(8));
t100 = -mrSges(2,1) + (-m(4) - m(5)) * pkin(1) + t103 * t76 + t104 * t74;
t73 = sin(pkin(9));
t93 = pkin(3) * t73;
t95 = m(5) * (qJ(2) + t93) - mrSges(2,2) + mrSges(3,3) + t61 * mrSges(6,1) + t62 * mrSges(6,2);
t79 = cos(qJ(1));
t91 = t76 * t79;
t64 = sin(t72);
t78 = sin(qJ(1));
t90 = t78 * t64;
t89 = t78 * t65;
t88 = t78 * t73;
t86 = t78 * t75;
t84 = qJ(2) * t79;
t69 = t78 * pkin(1);
t58 = pkin(4) * t64 + t93;
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) - m(5) + t101) * pkin(5) - t104 * t76 + (-mrSges(4,1) * t75 - mrSges(5,1) * t65 + mrSges(4,2) * t73 + mrSges(5,2) * t64 + t103) * t74) * g(3) + (-mrSges(1,2) - m(3) * (t69 - t84) + m(4) * t84 - m(6) * t69 + (-t86 * mrSges(4,1) - t89 * mrSges(5,1) + t88 * mrSges(4,2) + t90 * mrSges(5,2)) * t76 + (t73 * mrSges(4,1) + t75 * mrSges(4,2) + t64 * mrSges(5,1) + t65 * mrSges(5,2) - m(6) * (-qJ(2) - t58) + t95) * t79 + t100 * t78) * g(2) + (-mrSges(1,1) - (t75 * t91 + t88) * mrSges(4,1) - (-t73 * t91 + t86) * mrSges(4,2) - (t65 * t91 + t90) * mrSges(5,1) - (-t64 * t91 + t89) * mrSges(5,2) + (t99 * pkin(1) + t100) * t79 + (-m(6) * t58 + t101 * qJ(2) - t95) * t78) * g(1);
U = t1;
