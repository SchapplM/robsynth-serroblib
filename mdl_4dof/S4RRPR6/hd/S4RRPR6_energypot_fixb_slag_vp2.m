% Calculate potential energy for
% S4RRPR6
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR6_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:20
% EndTime: 2019-12-31 17:04:20
% DurationCPUTime: 0.14s
% Computational Cost: add. (83->40), mult. (86->32), div. (0->0), fcn. (61->8), ass. (0->16)
t36 = qJ(2) + pkin(7);
t32 = qJ(4) + t36;
t27 = sin(t32);
t28 = cos(t32);
t40 = cos(qJ(2));
t29 = t40 * pkin(2) + pkin(1);
t30 = sin(t36);
t31 = cos(t36);
t38 = sin(qJ(2));
t50 = -mrSges(2,1) - m(5) * (pkin(3) * t31 + t29) - t28 * mrSges(5,1) + t27 * mrSges(5,2) - m(4) * t29 - t31 * mrSges(4,1) + t30 * mrSges(4,2) - m(3) * pkin(1) - t40 * mrSges(3,1) + t38 * mrSges(3,2);
t37 = -qJ(3) - pkin(5);
t49 = mrSges(2,2) + m(5) * (-pkin(6) + t37) - mrSges(5,3) + m(4) * t37 - mrSges(4,3) - m(3) * pkin(5) - mrSges(3,3);
t48 = t38 * pkin(2) + pkin(4);
t41 = cos(qJ(1));
t39 = sin(qJ(1));
t1 = (-mrSges(1,3) - mrSges(2,3) - t38 * mrSges(3,1) - t40 * mrSges(3,2) - m(4) * t48 - t30 * mrSges(4,1) - t31 * mrSges(4,2) - m(5) * (pkin(3) * t30 + t48) - t27 * mrSges(5,1) - t28 * mrSges(5,2) + (-m(2) - m(3)) * pkin(4)) * g(3) + (t50 * t39 - t49 * t41 - mrSges(1,2)) * g(2) + (t49 * t39 + t50 * t41 - mrSges(1,1)) * g(1);
U = t1;
