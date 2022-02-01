% Calculate potential energy for
% S5RRPRR4
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
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:00
% EndTime: 2022-01-20 10:48:00
% DurationCPUTime: 0.18s
% Computational Cost: add. (136->49), mult. (92->40), div. (0->0), fcn. (61->10), ass. (0->20)
t60 = -m(4) - m(5) - m(6);
t46 = qJ(4) + qJ(5);
t39 = sin(t46);
t41 = cos(t46);
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t59 = -mrSges(4,1) - m(5) * pkin(3) - t50 * mrSges(5,1) + t48 * mrSges(5,2) - m(6) * (t50 * pkin(4) + pkin(3)) - t41 * mrSges(6,1) + t39 * mrSges(6,2);
t58 = m(5) * pkin(7) - m(6) * (-pkin(8) - pkin(7)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t57 = pkin(6) + pkin(5);
t49 = sin(qJ(1));
t43 = t49 * pkin(1);
t51 = cos(qJ(1));
t44 = t51 * pkin(1);
t47 = qJ(1) + qJ(2);
t42 = cos(t47);
t40 = sin(t47);
t38 = pkin(9) + t47;
t34 = cos(t38);
t33 = sin(t38);
t1 = (-m(2) * pkin(5) - m(3) * t57 - t39 * mrSges(6,1) - t50 * mrSges(5,2) - t41 * mrSges(6,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t48 + t60 * (qJ(3) + t57)) * g(3) + (-m(3) * t43 - t49 * mrSges(2,1) - t40 * mrSges(3,1) - t51 * mrSges(2,2) - t42 * mrSges(3,2) - mrSges(1,2) + t60 * (pkin(2) * t40 + t43) + t58 * t34 + t59 * t33) * g(2) + (-m(3) * t44 - t51 * mrSges(2,1) - t42 * mrSges(3,1) + t49 * mrSges(2,2) + t40 * mrSges(3,2) - mrSges(1,1) + t60 * (pkin(2) * t42 + t44) + t59 * t34 - t58 * t33) * g(1);
U = t1;
