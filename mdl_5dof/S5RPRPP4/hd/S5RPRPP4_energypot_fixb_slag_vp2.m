% Calculate potential energy for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:14
% EndTime: 2019-12-31 18:14:14
% DurationCPUTime: 0.30s
% Computational Cost: add. (100->51), mult. (117->39), div. (0->0), fcn. (86->6), ass. (0->21)
t63 = mrSges(5,1) + mrSges(6,1);
t62 = mrSges(5,2) - mrSges(6,3);
t61 = m(3) + m(4);
t60 = -m(5) - m(6);
t41 = qJ(3) + pkin(7);
t35 = sin(t41);
t36 = cos(t41);
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t59 = t43 * mrSges(4,1) + t45 * mrSges(4,2) + t63 * t35 + t62 * t36 - mrSges(2,2) + mrSges(3,3);
t58 = -m(4) * pkin(6) - mrSges(2,1) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t57 = pkin(2) + pkin(5);
t56 = pkin(3) * t43;
t44 = sin(qJ(1));
t46 = cos(qJ(1));
t55 = t46 * pkin(1) + t44 * qJ(2);
t53 = -qJ(2) - t56;
t48 = pkin(4) * t35 - qJ(5) * t36;
t42 = -qJ(4) - pkin(6);
t38 = t44 * pkin(1);
t1 = (-m(4) * t57 - t45 * mrSges(4,1) + t43 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t60 * (t45 * pkin(3) + t57) + (-m(6) * pkin(4) - t63) * t36 + (-m(6) * qJ(5) + t62) * t35 + (-m(2) - m(3)) * pkin(5)) * g(3) + (-mrSges(1,2) + t60 * (-t44 * t42 + t38) - t61 * t38 + (-m(5) * t53 - m(6) * (-t48 + t53) + t61 * qJ(2) + t59) * t46 + t58 * t44) * g(2) + (-mrSges(1,1) - t61 * t55 + t60 * (-t46 * t42 + t44 * t56 + t55) + t58 * t46 + (-m(6) * t48 - t59) * t44) * g(1);
U = t1;
