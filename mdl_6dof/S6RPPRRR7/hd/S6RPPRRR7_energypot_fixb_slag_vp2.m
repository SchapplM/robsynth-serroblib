% Calculate potential energy for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:32:53
% EndTime: 2019-03-09 02:32:53
% DurationCPUTime: 0.40s
% Computational Cost: add. (172->64), mult. (165->49), div. (0->0), fcn. (132->10), ass. (0->26)
t83 = m(6) + m(7);
t61 = sin(qJ(6));
t63 = cos(qJ(6));
t82 = -m(7) * pkin(5) - t63 * mrSges(7,1) + t61 * mrSges(7,2) - mrSges(6,1);
t80 = m(3) + m(4) + m(5);
t79 = -m(7) * pkin(9) - mrSges(7,3);
t57 = pkin(10) + qJ(4);
t51 = qJ(5) + t57;
t47 = sin(t51);
t48 = cos(t51);
t49 = sin(t57);
t50 = cos(t57);
t58 = sin(pkin(10));
t59 = cos(pkin(10));
t75 = pkin(3) * t58;
t78 = -m(5) * t75 - mrSges(4,1) * t58 - t49 * mrSges(5,1) - mrSges(4,2) * t59 - t50 * mrSges(5,2) - t48 * mrSges(6,2) + t82 * t47 + mrSges(2,2) - mrSges(3,3);
t60 = -pkin(7) - qJ(3);
t77 = -m(4) * qJ(3) + m(5) * t60 - t61 * mrSges(7,1) - t63 * mrSges(7,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) + t83 * (-pkin(8) + t60);
t76 = pkin(2) + pkin(6);
t62 = sin(qJ(1));
t64 = cos(qJ(1));
t74 = t64 * pkin(1) + t62 * qJ(2);
t45 = pkin(4) * t49 + t75;
t73 = -qJ(2) - t45;
t72 = t59 * pkin(3) + t76;
t1 = (-m(4) * t76 - m(5) * t72 - t59 * mrSges(4,1) - t50 * mrSges(5,1) + t58 * mrSges(4,2) + t49 * mrSges(5,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) - t83 * (pkin(4) * t50 + t72) + t82 * t48 + (mrSges(6,2) + t79) * t47 + (-m(2) - m(3)) * pkin(6)) * g(3) + (-mrSges(1,2) + (-m(6) * t73 - m(7) * (pkin(9) * t48 + t73) - t48 * mrSges(7,3) + t80 * qJ(2) - t78) * t64 + ((-t80 - t83) * pkin(1) + t77) * t62) * g(2) + (-mrSges(1,1) - t83 * (t62 * t45 + t74) - t80 * t74 + t77 * t64 + (-t79 * t48 + t78) * t62) * g(1);
U  = t1;
