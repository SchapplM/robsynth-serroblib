% Calculate potential energy for
% S4RRRR3
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
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:14
% EndTime: 2019-12-31 17:24:14
% DurationCPUTime: 0.20s
% Computational Cost: add. (98->45), mult. (91->35), div. (0->0), fcn. (61->8), ass. (0->19)
t25 = -m(2) - m(3);
t13 = sin(qJ(2));
t15 = cos(qJ(2));
t12 = qJ(2) + qJ(3);
t7 = qJ(4) + t12;
t2 = sin(t7);
t3 = cos(t7);
t4 = t15 * pkin(2) + pkin(1);
t5 = sin(t12);
t6 = cos(t12);
t24 = -m(3) * pkin(1) - t15 * mrSges(3,1) + t13 * mrSges(3,2) - mrSges(2,1) - m(5) * (pkin(3) * t6 + t4) - t3 * mrSges(5,1) + t2 * mrSges(5,2) - m(4) * t4 - t6 * mrSges(4,1) + t5 * mrSges(4,2);
t23 = -m(1) - m(4) - m(5) + t25;
t17 = -pkin(6) - pkin(5);
t22 = m(3) * pkin(5) - m(4) * t17 - m(5) * (-pkin(7) + t17) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t10 = pkin(4) + r_base(3);
t21 = t13 * pkin(2) + t10;
t16 = cos(qJ(1));
t14 = sin(qJ(1));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - mrSges(3,1) * t13 - mrSges(3,2) * t15 - m(4) * t21 - t5 * mrSges(4,1) - t6 * mrSges(4,2) - m(5) * (pkin(3) * t5 + t21) - t2 * mrSges(5,1) - t3 * mrSges(5,2) + t25 * t10) * g(3) + (t24 * t14 + t22 * t16 + t23 * r_base(2) - mrSges(1,2)) * g(2) + (-t22 * t14 + t24 * t16 + t23 * r_base(1) - mrSges(1,1)) * g(1);
U = t1;
