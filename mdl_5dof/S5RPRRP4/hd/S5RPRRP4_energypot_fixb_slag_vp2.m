% Calculate potential energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:05:25
% EndTime: 2019-12-05 18:05:25
% DurationCPUTime: 0.44s
% Computational Cost: add. (142->63), mult. (187->59), div. (0->0), fcn. (168->8), ass. (0->29)
t62 = sin(qJ(3));
t64 = cos(qJ(3));
t87 = -m(4) * pkin(2) - t64 * mrSges(4,1) + t62 * mrSges(4,2) - mrSges(3,1);
t86 = mrSges(3,2) - mrSges(5,3) - mrSges(6,3);
t85 = -mrSges(5,1) - mrSges(6,1);
t84 = mrSges(5,2) + mrSges(6,2);
t83 = -m(3) - m(4) - m(5) - m(6);
t82 = m(4) * pkin(6) + mrSges(4,3);
t60 = sin(pkin(8));
t61 = cos(pkin(8));
t81 = t86 * t60 + t87 * t61 - mrSges(2,1);
t59 = qJ(3) + qJ(4);
t52 = sin(t59);
t79 = t62 * pkin(3);
t80 = -m(5) * t79 - m(6) * (pkin(4) * t52 + t79) + mrSges(2,2) - mrSges(3,3) - t62 * mrSges(4,1) - t64 * mrSges(4,2);
t66 = -pkin(7) - pkin(6);
t51 = t64 * pkin(3) + pkin(2);
t63 = sin(qJ(1));
t76 = t63 * t52;
t53 = cos(t59);
t75 = t63 * t53;
t65 = cos(qJ(1));
t74 = t65 * t52;
t73 = t65 * t53;
t49 = pkin(4) * t53 + t51;
t58 = -qJ(5) + t66;
t69 = t49 * t61 - t58 * t60;
t68 = t51 * t61 - t60 * t66;
t1 = (-mrSges(1,3) + t85 * (t61 * t73 + t76) - t84 * (-t61 * t74 + t75) + t83 * (t65 * pkin(1) + t63 * qJ(2)) + t80 * t63 + (-m(5) * t68 - m(6) * t69 - t82 * t60 + t81) * t65) * g(3) + (-mrSges(1,2) + t85 * (-t61 * t75 + t74) - t84 * (t61 * t76 + t73) + (m(3) * pkin(1) - m(4) * (-pkin(6) * t60 - pkin(1)) + t60 * mrSges(4,3) - m(5) * (-pkin(1) - t68) - m(6) * (-pkin(1) - t69) - t81) * t63 + (t83 * qJ(2) + t80) * t65) * g(2) + (-mrSges(1,1) - mrSges(2,3) + (-m(2) + t83) * pkin(5) + (-m(5) * t66 - m(6) * t58 + t82 - t86) * t61 + (-m(5) * t51 - m(6) * t49 + t84 * t52 + t85 * t53 + t87) * t60) * g(1);
U = t1;
