% Calculate potential energy for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:41
% EndTime: 2019-12-31 16:26:41
% DurationCPUTime: 0.19s
% Computational Cost: add. (78->32), mult. (74->22), div. (0->0), fcn. (49->6), ass. (0->15)
t46 = -m(4) - m(5);
t45 = mrSges(4,2) + mrSges(5,2);
t44 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t42 = -m(3) + t46;
t43 = pkin(1) * t42 - mrSges(2,1);
t35 = sin(qJ(3));
t36 = cos(qJ(3));
t41 = t46 * pkin(2) + t45 * t35 + t44 * t36 - mrSges(3,1);
t39 = m(4) * pkin(5) - m(5) * (-qJ(4) - pkin(5)) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t32 = cos(pkin(6));
t31 = sin(pkin(6));
t30 = pkin(6) + qJ(2);
t27 = cos(t30);
t26 = sin(t30);
t1 = (-m(2) * qJ(1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t45 * t36 + t44 * t35 + t42 * (pkin(4) + qJ(1))) * g(3) + (-t32 * mrSges(2,2) + t41 * t26 + t39 * t27 + t43 * t31 - mrSges(1,2)) * g(2) + (t31 * mrSges(2,2) - t39 * t26 + t41 * t27 + t43 * t32 - mrSges(1,1)) * g(1);
U = t1;
