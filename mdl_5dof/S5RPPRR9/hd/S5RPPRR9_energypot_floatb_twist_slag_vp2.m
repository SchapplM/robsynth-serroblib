% Calculate potential energy for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:15
% EndTime: 2019-12-31 18:02:15
% DurationCPUTime: 0.39s
% Computational Cost: add. (132->55), mult. (189->45), div. (0->0), fcn. (189->8), ass. (0->24)
t42 = m(5) + m(6);
t15 = sin(qJ(5));
t17 = cos(qJ(5));
t41 = m(6) * pkin(4) + mrSges(6,1) * t17 - mrSges(6,2) * t15 + mrSges(5,1);
t40 = -m(6) * pkin(7) + mrSges(5,2) - mrSges(6,3);
t39 = -m(1) - m(2);
t37 = -mrSges(2,1) - mrSges(3,1);
t36 = mrSges(2,2) - mrSges(3,3);
t16 = sin(qJ(4));
t18 = cos(qJ(4));
t35 = -t16 * t40 + t18 * t41 + mrSges(4,1);
t33 = t15 * mrSges(6,1) + t17 * mrSges(6,2) + pkin(6) * t42 - mrSges(4,2) + mrSges(5,3);
t32 = cos(qJ(1));
t31 = sin(qJ(1));
t30 = cos(pkin(8));
t29 = sin(pkin(8));
t14 = pkin(5) + r_base(3);
t28 = t32 * pkin(1) + t31 * qJ(2) + r_base(1);
t27 = t32 * pkin(2) + t28;
t24 = t31 * pkin(1) - qJ(2) * t32 + r_base(2);
t21 = t31 * pkin(2) + t24;
t4 = t29 * t32 - t30 * t31;
t3 = -t29 * t31 - t30 * t32;
t1 = (-m(1) * r_base(3) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + (-m(4) - t42) * (-qJ(3) + t14) + t40 * t18 + t41 * t16 + (-m(2) - m(3)) * t14) * g(3) + (-m(3) * t24 - m(4) * t21 - mrSges(1,2) + t39 * r_base(2) - t36 * t32 + t37 * t31 - t42 * (-t4 * pkin(3) + t21) + t35 * t4 + t33 * t3) * g(2) + (-m(3) * t28 - m(4) * t27 - mrSges(1,1) + t39 * r_base(1) + t37 * t32 + t36 * t31 - t42 * (-t3 * pkin(3) + t27) - t33 * t4 + t35 * t3) * g(1);
U = t1;
