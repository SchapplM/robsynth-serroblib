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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:25:04
% EndTime: 2020-01-03 11:25:05
% DurationCPUTime: 0.40s
% Computational Cost: add. (147->60), mult. (148->59), div. (0->0), fcn. (125->8), ass. (0->27)
t77 = -mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t76 = m(6) * pkin(4);
t75 = -mrSges(5,1) - mrSges(6,1);
t74 = -mrSges(3,2) + mrSges(4,3);
t73 = mrSges(5,2) + mrSges(6,2);
t72 = -m(4) - m(5) - m(6);
t51 = sin(pkin(8));
t52 = cos(pkin(8));
t71 = mrSges(4,1) * t52 + t77 * t51 + mrSges(3,1);
t58 = cos(qJ(1));
t70 = pkin(1) * t58;
t56 = sin(qJ(1));
t49 = t56 * pkin(1);
t50 = qJ(1) + pkin(7);
t47 = sin(t50);
t55 = sin(qJ(4));
t69 = t47 * t55;
t48 = cos(t50);
t68 = t48 * t55;
t65 = t52 * t55;
t57 = cos(qJ(4));
t64 = t52 * t57;
t62 = pkin(3) * t52 + pkin(6) * t51;
t46 = pkin(4) * t57 + pkin(3);
t53 = -qJ(5) - pkin(6);
t59 = t46 * t52 - t51 * t53;
t1 = (t69 * t76 + m(3) * t70 + mrSges(2,1) * t58 - t56 * mrSges(2,2) - mrSges(1,3) + t72 * (-t47 * qJ(3) - t70) + t74 * t47 + t75 * (-t48 * t64 - t69) - t73 * (-t47 * t57 + t48 * t65) + (m(4) * pkin(2) - m(5) * (-pkin(2) - t62) - m(6) * (-pkin(2) - t59) + t71) * t48) * g(3) + (t68 * t76 - m(3) * t49 - t56 * mrSges(2,1) - mrSges(2,2) * t58 - mrSges(1,2) + t72 * (t47 * pkin(2) - qJ(3) * t48 + t49) + t74 * t48 + t75 * (t47 * t64 - t68) - t73 * (-t47 * t65 - t48 * t57) + (-m(5) * t62 - m(6) * t59 - t71) * t47) * g(2) + (-m(2) * pkin(5) - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) + (-m(3) + t72) * (qJ(2) + pkin(5)) + (m(5) * pkin(6) - m(6) * t53 + t77) * t52 + (-m(5) * pkin(3) - m(6) * t46 + t73 * t55 + t75 * t57 - mrSges(4,1)) * t51) * g(1);
U = t1;
