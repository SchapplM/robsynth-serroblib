% Calculate potential energy for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:45
% EndTime: 2019-12-31 16:20:46
% DurationCPUTime: 0.15s
% Computational Cost: add. (84->36), mult. (74->28), div. (0->0), fcn. (49->8), ass. (0->15)
t45 = -m(3) - m(4) - m(5);
t46 = pkin(1) * t45 - mrSges(2,1);
t32 = pkin(7) + qJ(4);
t26 = sin(t32);
t28 = cos(t32);
t34 = sin(pkin(7));
t36 = cos(pkin(7));
t44 = -mrSges(3,1) - m(4) * pkin(2) - mrSges(4,1) * t36 + mrSges(4,2) * t34 - m(5) * (pkin(3) * t36 + pkin(2)) - mrSges(5,1) * t28 + mrSges(5,2) * t26;
t42 = m(4) * qJ(3) - m(5) * (-pkin(5) - qJ(3)) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t37 = cos(pkin(6));
t35 = sin(pkin(6));
t33 = pkin(6) + qJ(2);
t29 = cos(t33);
t27 = sin(t33);
t1 = (-m(2) * qJ(1) - t26 * mrSges(5,1) - t36 * mrSges(4,2) - t28 * mrSges(5,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(5) * pkin(3) - mrSges(4,1)) * t34 + t45 * (pkin(4) + qJ(1))) * g(3) + (-t37 * mrSges(2,2) + t44 * t27 + t42 * t29 + t46 * t35 - mrSges(1,2)) * g(2) + (t35 * mrSges(2,2) - t42 * t27 + t44 * t29 + t46 * t37 - mrSges(1,1)) * g(1);
U = t1;
