% Calculate potential energy for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRPR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRPR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR7_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:33
% EndTime: 2019-12-31 16:25:34
% DurationCPUTime: 0.32s
% Computational Cost: add. (84->51), mult. (122->44), div. (0->0), fcn. (96->6), ass. (0->22)
t37 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t12 = sin(qJ(4));
t14 = cos(qJ(4));
t36 = -t12 * mrSges(5,1) - t14 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t33 = -m(1) - m(2);
t32 = m(4) + m(5);
t10 = sin(pkin(6));
t13 = sin(qJ(2));
t23 = qJ(3) * t13;
t15 = cos(qJ(2));
t28 = t10 * t15;
t31 = pkin(2) * t28 + t10 * t23;
t30 = t14 * mrSges(5,1) - t12 * mrSges(5,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t29 = t36 * t13 + t37 * t15 - mrSges(2,1);
t11 = cos(pkin(6));
t27 = t11 * t15;
t22 = t10 * pkin(1) + r_base(2);
t9 = qJ(1) + r_base(3);
t21 = t11 * pkin(1) + t10 * pkin(4) + r_base(1);
t19 = pkin(2) * t27 + t11 * t23 + t21;
t18 = -t11 * pkin(4) + t22;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t9 - t32 * (t13 * pkin(2) + t9) + (t32 * qJ(3) - t36) * t15 + (-m(5) * pkin(5) + t37) * t13) * g(3) + (-mrSges(1,2) - m(3) * t18 - m(4) * (t18 + t31) - m(5) * (pkin(5) * t28 + t22 + t31) + t33 * r_base(2) + (-m(5) * (-pkin(3) - pkin(4)) + t30) * t11 + t29 * t10) * g(2) + (-mrSges(1,1) - m(3) * t21 - m(4) * t19 - m(5) * (pkin(5) * t27 + t19) + t33 * r_base(1) + (-m(5) * pkin(3) - t30) * t10 + t29 * t11) * g(1);
U = t1;
