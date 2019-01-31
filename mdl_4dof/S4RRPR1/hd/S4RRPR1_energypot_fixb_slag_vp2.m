% Calculate potential energy for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-01-31 13:16
% Revision: 9ef80adae39e3cd5824e7abdb6e4e1e7895c437e (2019-01-31)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR1_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-31 13:16:33
% EndTime: 2019-01-31 13:16:33
% DurationCPUTime: 0.08s
% Computational Cost: add. (76->41), mult. (53->34), div. (0->0), fcn. (28->8), ass. (0->17)
t41 = pkin(5) + pkin(4);
t34 = qJ(1) + qJ(2);
t30 = sin(t34);
t35 = sin(qJ(1));
t40 = t35 * pkin(1) + pkin(2) * t30;
t31 = cos(t34);
t36 = cos(qJ(1));
t39 = t36 * pkin(1) + pkin(2) * t31;
t38 = qJ(3) + t41;
t37 = -m(3) * pkin(1) - mrSges(2,1);
t29 = pkin(7) + t34;
t28 = qJ(4) + t29;
t25 = cos(t29);
t24 = sin(t29);
t23 = cos(t28);
t22 = sin(t28);
t1 = (-mrSges(1,3) - m(2) * pkin(4) - mrSges(2,3) - m(3) * t41 - mrSges(3,3) - m(4) * t38 - mrSges(4,3) - m(5) * (pkin(6) + t38) - mrSges(5,3)) * g(3) + (-mrSges(1,2) - t36 * mrSges(2,2) - t30 * mrSges(3,1) - t31 * mrSges(3,2) - m(4) * t40 - t24 * mrSges(4,1) - t25 * mrSges(4,2) - m(5) * (pkin(3) * t24 + t40) - t22 * mrSges(5,1) - t23 * mrSges(5,2) + t37 * t35) * g(2) + (-mrSges(1,1) + t35 * mrSges(2,2) - t31 * mrSges(3,1) + t30 * mrSges(3,2) - m(4) * t39 - t25 * mrSges(4,1) + t24 * mrSges(4,2) - m(5) * (pkin(3) * t25 + t39) - t23 * mrSges(5,1) + t22 * mrSges(5,2) + t37 * t36) * g(1);
U  = t1;
