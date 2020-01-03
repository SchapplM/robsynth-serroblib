% Calculate potential energy for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR13_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:03
% EndTime: 2019-12-31 20:32:03
% DurationCPUTime: 0.34s
% Computational Cost: add. (152->49), mult. (187->41), div. (0->0), fcn. (168->10), ass. (0->21)
t58 = pkin(9) + qJ(4);
t52 = qJ(5) + t58;
t47 = sin(t52);
t48 = cos(t52);
t60 = cos(pkin(9));
t49 = t60 * pkin(3) + pkin(2);
t50 = sin(t58);
t51 = cos(t58);
t59 = sin(pkin(9));
t85 = -mrSges(3,1) - m(4) * pkin(2) - t60 * mrSges(4,1) + t59 * mrSges(4,2) - m(5) * t49 - t51 * mrSges(5,1) + t50 * mrSges(5,2) - m(6) * (pkin(4) * t51 + t49) - t48 * mrSges(6,1) + t47 * mrSges(6,2);
t61 = -pkin(7) - qJ(3);
t84 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) + m(6) * (-pkin(8) + t61) - mrSges(6,3) + m(5) * t61 - mrSges(5,3);
t83 = m(3) + m(4) + m(5) + m(6);
t62 = sin(qJ(2));
t64 = cos(qJ(2));
t79 = t84 * t62 + t85 * t64 - mrSges(2,1);
t77 = pkin(3) * t59;
t78 = m(5) * t77 + m(6) * (pkin(4) * t50 + t77) - mrSges(2,2) + mrSges(3,3) + t47 * mrSges(6,1) + t48 * mrSges(6,2) + t50 * mrSges(5,1) + t51 * mrSges(5,2) + t59 * mrSges(4,1) + t60 * mrSges(4,2);
t65 = cos(qJ(1));
t63 = sin(qJ(1));
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) - t83) * pkin(5) - t84 * t64 + t85 * t62) * g(3) + (-mrSges(1,2) + (t83 * pkin(6) + t78) * t65 + (-t83 * pkin(1) + t79) * t63) * g(2) + (-mrSges(1,1) - t83 * (t65 * pkin(1) + t63 * pkin(6)) + t79 * t65 - t78 * t63) * g(1);
U = t1;
