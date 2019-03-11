% Calculate potential energy for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRPR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR1_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:20:51
% EndTime: 2019-03-08 18:20:51
% DurationCPUTime: 0.11s
% Computational Cost: add. (76->38), mult. (65->32), div. (0->0), fcn. (44->6), ass. (0->15)
t45 = m(4) + m(5);
t44 = mrSges(3,2) - mrSges(4,3);
t43 = pkin(4) + qJ(1);
t40 = -m(3) * pkin(1) - mrSges(2,1);
t39 = -m(5) * pkin(3) - mrSges(3,1) - mrSges(4,1);
t38 = cos(qJ(4));
t37 = sin(qJ(4));
t35 = cos(pkin(6));
t34 = sin(pkin(6));
t33 = pkin(6) + qJ(2);
t30 = cos(t33);
t29 = sin(t33);
t25 = t29 * t38 - t30 * t37;
t24 = -t29 * t37 - t30 * t38;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - mrSges(3,3) - mrSges(4,2) - m(5) * (-pkin(5) + t43) + mrSges(5,3) + (-m(3) - m(4)) * t43) * g(3) + (-mrSges(1,2) - t35 * mrSges(2,2) - t25 * mrSges(5,1) - t24 * mrSges(5,2) + t40 * t34 + (qJ(3) * t45 - t44) * t30 + t39 * t29 - t45 * (pkin(1) * t34 + pkin(2) * t29)) * g(2) + (t24 * mrSges(5,1) + t34 * mrSges(2,2) - t25 * mrSges(5,2) + t44 * t29 + t39 * t30 + t40 * t35 - mrSges(1,1) - t45 * (pkin(1) * t35 + pkin(2) * t30 + qJ(3) * t29)) * g(1);
U  = t1;
