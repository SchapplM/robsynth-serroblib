% Calculate potential energy for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPPR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:43
% EndTime: 2019-12-31 17:36:43
% DurationCPUTime: 0.31s
% Computational Cost: add. (138->53), mult. (138->47), div. (0->0), fcn. (113->8), ass. (0->24)
t84 = mrSges(4,2) - mrSges(5,3);
t58 = sin(pkin(8));
t60 = cos(pkin(8));
t83 = pkin(3) * t60 + qJ(4) * t58;
t82 = -m(6) * pkin(4) - mrSges(4,1) - mrSges(5,1);
t80 = -m(5) - m(6);
t57 = pkin(7) + qJ(2);
t52 = sin(t57);
t79 = t83 * t52;
t63 = sin(qJ(5));
t64 = cos(qJ(5));
t67 = t58 * t63 + t60 * t64;
t68 = t58 * t64 - t60 * t63;
t78 = -t67 * mrSges(6,1) - t68 * mrSges(6,2) + t84 * t58 + t82 * t60 - mrSges(3,1);
t77 = mrSges(5,2) + mrSges(4,3) - mrSges(3,2) - mrSges(6,3);
t59 = sin(pkin(7));
t55 = t59 * pkin(1);
t61 = cos(pkin(7));
t56 = t61 * pkin(1);
t62 = pkin(5) + qJ(1);
t75 = t52 * pkin(2) + t55;
t53 = cos(t57);
t72 = -qJ(3) * t53 + t75;
t1 = (-m(2) * qJ(1) - t68 * mrSges(6,1) + t67 * mrSges(6,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t80 * (t58 * pkin(3) - qJ(4) * t60 + t62) + (-m(3) - m(4)) * t62 - t84 * t60 + t82 * t58) * g(3) + (-mrSges(1,2) - t59 * mrSges(2,1) - t61 * mrSges(2,2) - m(3) * t55 - m(4) * t72 - m(5) * (t72 + t79) - m(6) * (t75 + t79) + (-m(6) * (pkin(6) - qJ(3)) + t77) * t53 + t78 * t52) * g(2) + (-m(3) * t56 - t61 * mrSges(2,1) + t59 * mrSges(2,2) - mrSges(1,1) + (m(6) * pkin(6) - t77) * t52 + (-m(4) + t80) * (t53 * pkin(2) + t52 * qJ(3) + t56) + (t80 * t83 + t78) * t53) * g(1);
U = t1;
