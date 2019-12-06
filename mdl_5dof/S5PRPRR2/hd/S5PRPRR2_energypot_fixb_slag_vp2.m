% Calculate potential energy for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:42
% EndTime: 2019-12-05 15:44:43
% DurationCPUTime: 0.25s
% Computational Cost: add. (150->50), mult. (137->41), div. (0->0), fcn. (110->10), ass. (0->23)
t60 = sin(qJ(5));
t62 = cos(qJ(5));
t80 = -m(6) * pkin(4) - t62 * mrSges(6,1) + t60 * mrSges(6,2) - mrSges(5,1);
t79 = -m(6) * pkin(7) + mrSges(5,2) - mrSges(6,3);
t78 = m(5) + m(6);
t56 = qJ(2) + pkin(9);
t52 = qJ(4) + t56;
t47 = sin(t52);
t48 = cos(t52);
t63 = cos(qJ(2));
t49 = t63 * pkin(2) + pkin(1);
t50 = sin(t56);
t51 = cos(t56);
t61 = sin(qJ(2));
t76 = -m(3) * pkin(1) - m(4) * t49 - t63 * mrSges(3,1) - t51 * mrSges(4,1) + t61 * mrSges(3,2) + t50 * mrSges(4,2) + t79 * t47 + t80 * t48 - mrSges(2,1);
t59 = -qJ(3) - pkin(5);
t75 = m(3) * pkin(5) - m(4) * t59 + t60 * mrSges(6,1) + t62 * mrSges(6,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t73 = t61 * pkin(2) + qJ(1);
t58 = cos(pkin(8));
t57 = sin(pkin(8));
t55 = -pkin(6) + t59;
t44 = pkin(3) * t51 + t49;
t1 = (-m(4) * t73 - t61 * mrSges(3,1) - t50 * mrSges(4,1) - t63 * mrSges(3,2) - t51 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - t78 * (pkin(3) * t50 + t73) - t79 * t48 + t80 * t47 + (-m(2) - m(3)) * qJ(1)) * g(3) + (-mrSges(1,2) - t78 * (t57 * t44 + t58 * t55) + t75 * t58 + t76 * t57) * g(2) + (-mrSges(1,1) + (t78 * t55 - t75) * t57 + (-t78 * t44 + t76) * t58) * g(1);
U = t1;
