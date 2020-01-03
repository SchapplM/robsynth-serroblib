% Calculate potential energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:20
% EndTime: 2020-01-03 11:30:21
% DurationCPUTime: 0.42s
% Computational Cost: add. (152->49), mult. (187->38), div. (0->0), fcn. (168->10), ass. (0->22)
t52 = pkin(9) + qJ(4);
t48 = qJ(5) + t52;
t43 = sin(t48);
t44 = cos(t48);
t55 = cos(pkin(9));
t45 = t55 * pkin(3) + pkin(2);
t46 = sin(t52);
t47 = cos(t52);
t53 = sin(pkin(9));
t77 = mrSges(3,1) + m(4) * pkin(2) + mrSges(4,1) * t55 - mrSges(4,2) * t53 + m(5) * t45 + mrSges(5,1) * t47 - mrSges(5,2) * t46 + m(6) * (pkin(4) * t47 + t45) + mrSges(6,1) * t44 - mrSges(6,2) * t43;
t57 = -pkin(6) - qJ(3);
t76 = -mrSges(3,2) + m(4) * qJ(3) - m(5) * t57 - m(6) * (-pkin(7) + t57) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t75 = m(3) + m(4);
t54 = sin(pkin(8));
t56 = cos(pkin(8));
t70 = m(5) + m(6) + t75;
t74 = t70 * pkin(1) + t76 * t54 + t77 * t56 + mrSges(2,1);
t72 = -pkin(3) * t53 - qJ(2);
t69 = -m(5) * t72 + t46 * mrSges(5,1) + t47 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3) - m(6) * (-pkin(4) * t46 + t72) + t43 * mrSges(6,1) + t44 * mrSges(6,2) + t53 * mrSges(4,1) + t55 * mrSges(4,2) + t75 * qJ(2);
t59 = cos(qJ(1));
t58 = sin(qJ(1));
t1 = (t69 * t58 + t74 * t59 - mrSges(1,3)) * g(3) + (-t74 * t58 + t69 * t59 - mrSges(1,2)) * g(2) + (-mrSges(1,1) - mrSges(2,3) + (-m(2) - t70) * pkin(5) + t76 * t56 - t77 * t54) * g(1);
U = t1;
