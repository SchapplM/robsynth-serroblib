% Calculate potential energy for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:40
% EndTime: 2019-12-05 16:47:40
% DurationCPUTime: 0.37s
% Computational Cost: add. (142->55), mult. (187->52), div. (0->0), fcn. (168->8), ass. (0->24)
t65 = cos(qJ(3));
t52 = t65 * pkin(3) + pkin(2);
t60 = qJ(3) + qJ(4);
t54 = cos(t60);
t63 = sin(qJ(3));
t87 = -m(5) * t52 - m(6) * (pkin(4) * t54 + t52) - mrSges(3,1) - m(4) * pkin(2) - t65 * mrSges(4,1) + t63 * mrSges(4,2);
t67 = -pkin(7) - pkin(6);
t86 = m(5) * t67 + m(6) * (-qJ(5) + t67) + mrSges(3,2) - mrSges(5,3) - mrSges(6,3) - m(4) * pkin(6) - mrSges(4,3);
t85 = -mrSges(5,1) - mrSges(6,1);
t84 = mrSges(5,2) + mrSges(6,2);
t83 = -m(3) - m(5) - m(6);
t82 = -m(4) + t83;
t64 = sin(qJ(2));
t66 = cos(qJ(2));
t80 = t86 * t64 + t87 * t66 - mrSges(2,1);
t53 = sin(t60);
t78 = pkin(3) * t63;
t79 = -m(5) * t78 - m(6) * (pkin(4) * t53 + t78) + mrSges(2,2) - mrSges(3,3) - t63 * mrSges(4,1) - t65 * mrSges(4,2);
t61 = sin(pkin(8));
t77 = t61 * t66;
t62 = cos(pkin(8));
t76 = t62 * t66;
t56 = t61 * pkin(1);
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) + t82) * qJ(1) - t86 * t66 + (t84 * t53 + t85 * t54 + t87) * t64) * g(3) + (-m(4) * t56 - mrSges(1,2) + t85 * (-t53 * t62 + t54 * t77) - t84 * (-t53 * t77 - t54 * t62) + t83 * (-t62 * pkin(5) + t56) + (m(4) * pkin(5) - t79) * t62 + t80 * t61) * g(2) + (-mrSges(1,1) + t85 * (t53 * t61 + t54 * t76) - t84 * (-t53 * t76 + t54 * t61) + t82 * (t62 * pkin(1) + t61 * pkin(5)) + t79 * t61 + t80 * t62) * g(1);
U = t1;
