% Calculate potential energy for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRPR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRPR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR6_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:16
% EndTime: 2019-12-31 16:24:17
% DurationCPUTime: 0.30s
% Computational Cost: add. (102->44), mult. (135->35), div. (0->0), fcn. (113->8), ass. (0->16)
t11 = cos(pkin(7));
t8 = pkin(7) + qJ(4);
t2 = sin(t8);
t3 = cos(t8);
t9 = sin(pkin(7));
t33 = -m(4) * pkin(2) - t11 * mrSges(4,1) + t9 * mrSges(4,2) - mrSges(3,1) - m(5) * (pkin(3) * t11 + pkin(2)) - t3 * mrSges(5,1) + t2 * mrSges(5,2);
t32 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) + m(5) * (-pkin(5) - qJ(3)) - mrSges(5,3);
t31 = -m(1) - m(2);
t30 = m(3) + m(4) + m(5);
t14 = sin(qJ(2));
t15 = cos(qJ(2));
t27 = t32 * t14 + t33 * t15 - mrSges(2,1);
t26 = t2 * mrSges(5,1) + t11 * mrSges(4,2) + t3 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3) + (m(5) * pkin(3) + mrSges(4,1)) * t9;
t12 = cos(pkin(6));
t10 = sin(pkin(6));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - t30) * (qJ(1) + r_base(3)) - t32 * t15 + t33 * t14) * g(3) + (-mrSges(1,2) + t31 * r_base(2) - t30 * (t10 * pkin(1) + r_base(2)) + (t30 * pkin(4) + t26) * t12 + t27 * t10) * g(2) + (-mrSges(1,1) + t31 * r_base(1) - t30 * (t12 * pkin(1) + t10 * pkin(4) + r_base(1)) + t27 * t12 - t26 * t10) * g(1);
U = t1;
