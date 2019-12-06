% Calculate potential energy for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:31
% EndTime: 2019-12-05 16:15:31
% DurationCPUTime: 0.17s
% Computational Cost: add. (136->49), mult. (92->40), div. (0->0), fcn. (61->10), ass. (0->20)
t60 = -m(4) - m(5) - m(6);
t46 = pkin(9) + qJ(5);
t38 = sin(t46);
t40 = cos(t46);
t48 = sin(pkin(9));
t50 = cos(pkin(9));
t59 = -mrSges(4,1) - m(5) * pkin(3) - t50 * mrSges(5,1) + t48 * mrSges(5,2) - m(6) * (t50 * pkin(4) + pkin(3)) - t40 * mrSges(6,1) + t38 * mrSges(6,2);
t58 = m(5) * qJ(4) - m(6) * (-pkin(7) - qJ(4)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t49 = sin(pkin(8));
t43 = t49 * pkin(1);
t51 = cos(pkin(8));
t44 = t51 * pkin(1);
t57 = pkin(5) + qJ(1);
t47 = pkin(8) + qJ(2);
t42 = qJ(3) + t47;
t41 = cos(t47);
t39 = sin(t47);
t36 = cos(t42);
t35 = sin(t42);
t1 = (-m(2) * qJ(1) - m(3) * t57 - t38 * mrSges(6,1) - t50 * mrSges(5,2) - t40 * mrSges(6,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t48 + t60 * (pkin(6) + t57)) * g(3) + (-m(3) * t43 - t49 * mrSges(2,1) - t39 * mrSges(3,1) - t51 * mrSges(2,2) - t41 * mrSges(3,2) - mrSges(1,2) + t60 * (pkin(2) * t39 + t43) + t58 * t36 + t59 * t35) * g(2) + (-m(3) * t44 - t51 * mrSges(2,1) - t41 * mrSges(3,1) + t49 * mrSges(2,2) + t39 * mrSges(3,2) - mrSges(1,1) + t60 * (pkin(2) * t41 + t44) + t59 * t36 - t58 * t35) * g(1);
U = t1;
