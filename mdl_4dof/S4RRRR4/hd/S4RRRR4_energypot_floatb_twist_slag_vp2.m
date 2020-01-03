% Calculate potential energy for
% S4RRRR4
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR4_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:38
% EndTime: 2019-12-31 17:25:38
% DurationCPUTime: 0.27s
% Computational Cost: add. (104->45), mult. (111->36), div. (0->0), fcn. (85->8), ass. (0->20)
t10 = sin(qJ(4));
t13 = cos(qJ(4));
t32 = -m(5) * pkin(3) - t13 * mrSges(5,1) + t10 * mrSges(5,2) - mrSges(4,1);
t31 = -m(5) * pkin(7) + mrSges(4,2) - mrSges(5,3);
t30 = -m(2) - m(3);
t29 = m(4) + m(5);
t28 = -m(1) + t30;
t11 = sin(qJ(2));
t14 = cos(qJ(2));
t9 = qJ(2) + qJ(3);
t5 = sin(t9);
t6 = cos(t9);
t26 = -m(3) * pkin(1) - t14 * mrSges(3,1) + t11 * mrSges(3,2) + t31 * t5 + t32 * t6 - mrSges(2,1);
t25 = m(3) * pkin(5) + t10 * mrSges(5,1) + t13 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t8 = pkin(4) + r_base(3);
t16 = -pkin(6) - pkin(5);
t15 = cos(qJ(1));
t12 = sin(qJ(1));
t3 = pkin(2) * t14 + pkin(1);
t1 = (-m(1) * r_base(3) - mrSges(3,1) * t11 - mrSges(3,2) * t14 - mrSges(1,3) - mrSges(2,3) + t30 * t8 - t29 * (t11 * pkin(2) + t8) - t31 * t6 + t32 * t5) * g(3) + (-mrSges(1,2) - t29 * (t12 * t3 + t15 * t16 + r_base(2)) + t28 * r_base(2) + t25 * t15 + t26 * t12) * g(2) + (-mrSges(1,1) - t29 * (t15 * t3 + r_base(1)) + t28 * r_base(1) + t26 * t15 + (t29 * t16 - t25) * t12) * g(1);
U = t1;
