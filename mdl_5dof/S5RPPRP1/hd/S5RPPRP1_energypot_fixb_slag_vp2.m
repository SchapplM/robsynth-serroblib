% Calculate potential energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:08
% EndTime: 2022-01-23 09:12:09
% DurationCPUTime: 0.36s
% Computational Cost: add. (147->53), mult. (148->51), div. (0->0), fcn. (125->8), ass. (0->24)
t62 = cos(qJ(4));
t80 = -m(5) * pkin(3) - m(6) * (pkin(4) * t62 + pkin(3)) - mrSges(4,1);
t79 = m(5) * pkin(6) - m(6) * (-qJ(5) - pkin(6)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t78 = m(6) * pkin(4);
t77 = -mrSges(5,1) - mrSges(6,1);
t76 = mrSges(3,2) - mrSges(4,3);
t75 = mrSges(5,2) + mrSges(6,2);
t74 = -m(4) - m(5) - m(6);
t56 = sin(pkin(8));
t57 = cos(pkin(8));
t73 = -t79 * t56 + t80 * t57 - mrSges(3,1);
t61 = sin(qJ(1));
t53 = t61 * pkin(1);
t63 = cos(qJ(1));
t54 = t63 * pkin(1);
t55 = qJ(1) + pkin(7);
t51 = sin(t55);
t60 = sin(qJ(4));
t72 = t51 * t60;
t52 = cos(t55);
t71 = t52 * t60;
t70 = t57 * t60;
t69 = t57 * t62;
t1 = (-m(2) * pkin(5) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) + t74) * (qJ(2) + pkin(5)) + t79 * t57 + (t75 * t60 + t77 * t62 + t80) * t56) * g(3) + (t71 * t78 - m(3) * t53 - t61 * mrSges(2,1) - t63 * mrSges(2,2) - mrSges(1,2) + t74 * (t51 * pkin(2) - qJ(3) * t52 + t53) - t76 * t52 + t77 * (t51 * t69 - t71) - t75 * (-t51 * t70 - t52 * t62) + t73 * t51) * g(2) + (-t72 * t78 - m(3) * t54 - t63 * mrSges(2,1) + t61 * mrSges(2,2) - mrSges(1,1) + t74 * (t52 * pkin(2) + t51 * qJ(3) + t54) + t76 * t51 + t77 * (t52 * t69 + t72) - t75 * (t51 * t62 - t52 * t70) + t73 * t52) * g(1);
U = t1;
