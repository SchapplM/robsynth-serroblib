% Calculate potential energy for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PPRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR1_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:27
% EndTime: 2019-03-08 18:15:27
% DurationCPUTime: 0.12s
% Computational Cost: add. (65->39), mult. (79->36), div. (0->0), fcn. (62->6), ass. (0->15)
t41 = m(3) + m(4) + m(5);
t37 = sin(qJ(3));
t40 = m(5) * pkin(3) * t37 - mrSges(2,2) + mrSges(3,3);
t38 = cos(qJ(3));
t39 = -m(4) * pkin(2) - m(5) * (t38 * pkin(3) + pkin(2)) - mrSges(2,1) - mrSges(3,1);
t36 = cos(pkin(6));
t35 = sin(pkin(6));
t34 = qJ(3) + qJ(4);
t31 = cos(t34);
t30 = sin(t34);
t27 = t35 * t38 - t36 * t37;
t26 = -t35 * t37 - t36 * t38;
t25 = -t36 * t30 + t35 * t31;
t24 = -t35 * t30 - t36 * t31;
t1 = (-mrSges(1,3) - mrSges(2,3) - mrSges(3,2) + m(4) * pkin(4) + mrSges(4,3) - m(5) * (-pkin(5) - pkin(4)) + mrSges(5,3) + (-m(2) - t41) * qJ(1)) * g(3) + (-t27 * mrSges(4,1) - t25 * mrSges(5,1) - t26 * mrSges(4,2) - t24 * mrSges(5,2) - mrSges(1,2) + (t41 * qJ(2) + t40) * t36 + (-t41 * pkin(1) + t39) * t35) * g(2) + (t26 * mrSges(4,1) + t24 * mrSges(5,1) - t27 * mrSges(4,2) - t25 * mrSges(5,2) - t40 * t35 + t39 * t36 - mrSges(1,1) - t41 * (t36 * pkin(1) + t35 * qJ(2))) * g(1);
U  = t1;
