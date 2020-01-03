% Calculate potential energy for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:02
% EndTime: 2019-12-31 19:38:02
% DurationCPUTime: 0.41s
% Computational Cost: add. (120->60), mult. (186->55), div. (0->0), fcn. (167->8), ass. (0->27)
t93 = mrSges(3,2) - mrSges(4,3);
t64 = sin(qJ(2));
t66 = cos(qJ(2));
t92 = pkin(2) * t66 + qJ(3) * t64;
t62 = cos(pkin(8));
t91 = -m(6) * (pkin(4) * t62 + pkin(3)) - mrSges(3,1) - mrSges(4,1) - m(5) * pkin(3);
t88 = -m(5) - m(6);
t90 = -m(4) + t88;
t65 = sin(qJ(1));
t87 = t92 * t65;
t61 = sin(pkin(8));
t82 = t61 * t64;
t70 = t62 * t66 + t82;
t71 = -t61 * t66 + t62 * t64;
t60 = pkin(8) + qJ(5);
t54 = sin(t60);
t55 = cos(t60);
t72 = t54 * t64 + t55 * t66;
t73 = -t54 * t66 + t55 * t64;
t86 = -m(6) * pkin(4) * t82 - t70 * mrSges(5,1) - t72 * mrSges(6,1) - t71 * mrSges(5,2) - t73 * mrSges(6,2) + t93 * t64 + t91 * t66 - mrSges(2,1);
t85 = mrSges(4,2) + mrSges(3,3) - mrSges(2,2) - mrSges(5,3) - mrSges(6,3);
t83 = t64 * pkin(2) + pkin(5);
t58 = t65 * pkin(1);
t67 = cos(qJ(1));
t78 = -pkin(6) * t67 + t58;
t63 = -pkin(7) - qJ(4);
t1 = (-m(6) * t83 - t71 * mrSges(5,1) - t73 * mrSges(6,1) + t70 * mrSges(5,2) + t72 * mrSges(6,2) - mrSges(1,3) - mrSges(2,3) + (-m(4) - m(5)) * (-qJ(3) * t66 + t83) + (-m(6) * (-pkin(4) * t61 - qJ(3)) - t93) * t66 + t91 * t64 + (-m(2) - m(3)) * pkin(5)) * g(3) + (-mrSges(1,2) - m(3) * t78 - m(4) * (t78 + t87) + t88 * (t58 + t87) + (-m(5) * (-pkin(6) + qJ(4)) - m(6) * (-pkin(6) - t63) + t85) * t67 + t86 * t65) * g(2) + (-mrSges(1,1) + (m(5) * qJ(4) - m(6) * t63 - t85) * t65 + (-m(3) + t90) * (t67 * pkin(1) + t65 * pkin(6)) + (t90 * t92 + t86) * t67) * g(1);
U = t1;
