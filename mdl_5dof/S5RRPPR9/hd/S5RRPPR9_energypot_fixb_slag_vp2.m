% Calculate potential energy for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:36
% EndTime: 2019-12-31 19:40:36
% DurationCPUTime: 0.32s
% Computational Cost: add. (99->57), mult. (166->48), div. (0->0), fcn. (139->6), ass. (0->25)
t87 = -m(6) * pkin(7) - mrSges(3,1) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t86 = t61 * mrSges(6,1) - t58 * mrSges(6,2) + mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t85 = m(5) + m(6);
t60 = sin(qJ(1));
t59 = sin(qJ(2));
t75 = qJ(3) * t59;
t62 = cos(qJ(2));
t78 = t60 * t62;
t84 = pkin(2) * t78 + t60 * t75;
t63 = cos(qJ(1));
t83 = pkin(3) * t78 + t63 * qJ(4);
t81 = -mrSges(2,1) + t87 * t62 + (-m(6) * pkin(4) - t86) * t59;
t80 = t58 * mrSges(6,1) + t61 * mrSges(6,2) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3);
t79 = t59 * pkin(2) + pkin(5);
t77 = t63 * t62;
t76 = t63 * pkin(1) + t60 * pkin(6);
t56 = t60 * pkin(1);
t74 = -pkin(6) * t63 + t56;
t73 = pkin(2) * t77 + t63 * t75 + t76;
t72 = -qJ(3) * t62 + t79;
t65 = t74 + t84;
t53 = t59 * pkin(3);
t1 = (-mrSges(1,3) - mrSges(2,3) - m(4) * t72 - m(5) * (t53 + t72) - m(6) * (t53 + t79) + (-m(2) - m(3)) * pkin(5) + (-m(6) * (-pkin(4) - qJ(3)) + t86) * t62 + t87 * t59) * g(3) + (-mrSges(1,2) - m(3) * t74 - m(4) * t65 - m(5) * (t65 + t83) - m(6) * (t56 + t83 + t84) + (m(6) * pkin(6) - t80) * t63 + t81 * t60) * g(2) + (-m(3) * t76 - m(4) * t73 - mrSges(1,1) - t85 * (pkin(3) * t77 + t73) + t81 * t63 + (t85 * qJ(4) + t80) * t60) * g(1);
U = t1;
