% Calculate potential energy for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP6_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:54
% EndTime: 2019-12-31 17:54:54
% DurationCPUTime: 0.31s
% Computational Cost: add. (100->51), mult. (117->39), div. (0->0), fcn. (86->6), ass. (0->21)
t63 = mrSges(5,1) + mrSges(6,1);
t62 = mrSges(5,2) - mrSges(6,3);
t61 = m(3) + m(4);
t60 = -m(5) - m(6);
t41 = pkin(7) + qJ(4);
t35 = sin(t41);
t36 = cos(t41);
t42 = sin(pkin(7));
t43 = cos(pkin(7));
t59 = t42 * mrSges(4,1) + t43 * mrSges(4,2) + t63 * t35 + t62 * t36 - mrSges(2,2) + mrSges(3,3);
t58 = -m(4) * qJ(3) - mrSges(2,1) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t57 = pkin(2) + pkin(5);
t56 = pkin(3) * t42;
t45 = sin(qJ(1));
t46 = cos(qJ(1));
t55 = t46 * pkin(1) + t45 * qJ(2);
t53 = -qJ(2) - t56;
t48 = pkin(4) * t35 - qJ(5) * t36;
t44 = -pkin(6) - qJ(3);
t39 = t45 * pkin(1);
t1 = (-m(4) * t57 - t43 * mrSges(4,1) + t42 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t60 * (t43 * pkin(3) + t57) + (-m(6) * pkin(4) - t63) * t36 + (-m(6) * qJ(5) + t62) * t35 + (-m(2) - m(3)) * pkin(5)) * g(3) + (-mrSges(1,2) + t60 * (-t45 * t44 + t39) - t61 * t39 + (-m(5) * t53 - m(6) * (-t48 + t53) + t61 * qJ(2) + t59) * t46 + t58 * t45) * g(2) + (-mrSges(1,1) - t61 * t55 + t60 * (-t46 * t44 + t45 * t56 + t55) + t58 * t46 + (-m(6) * t48 - t59) * t45) * g(1);
U = t1;
