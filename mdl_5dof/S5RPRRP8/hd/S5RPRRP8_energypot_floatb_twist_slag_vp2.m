% Calculate potential energy for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:04
% EndTime: 2019-12-31 18:47:04
% DurationCPUTime: 0.31s
% Computational Cost: add. (123->53), mult. (168->42), div. (0->0), fcn. (160->6), ass. (0->22)
t36 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t35 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t34 = -m(1) - m(2);
t33 = -m(5) - m(6);
t32 = -mrSges(2,1) - mrSges(3,1);
t31 = mrSges(2,2) - mrSges(3,3);
t15 = sin(qJ(4));
t16 = cos(qJ(4));
t30 = t35 * t15 + t36 * t16 + mrSges(4,1);
t29 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t28 = cos(qJ(3));
t27 = sin(qJ(1));
t26 = sin(qJ(3));
t14 = pkin(5) + r_base(3);
t17 = cos(qJ(1));
t25 = t17 * pkin(1) + t27 * qJ(2) + r_base(1);
t24 = t17 * pkin(2) + t25;
t22 = t27 * pkin(1) - t17 * qJ(2) + r_base(2);
t21 = t27 * pkin(2) + t22;
t4 = t17 * t26 - t27 * t28;
t3 = -t17 * t28 - t27 * t26;
t1 = (-m(1) * r_base(3) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + (-m(4) + t33) * (-pkin(6) + t14) - t35 * t16 + t36 * t15 + (-m(2) - m(3)) * t14) * g(3) + (-m(3) * t22 - m(4) * t21 - mrSges(1,2) + t34 * r_base(2) + t32 * t27 + t33 * (-t4 * pkin(3) - t3 * pkin(7) + t21) - t31 * t17 + t30 * t4 - t29 * t3) * g(2) + (-m(3) * t25 - m(4) * t24 - mrSges(1,1) + t34 * r_base(1) + t31 * t27 + t33 * (-t3 * pkin(3) + t4 * pkin(7) + t24) + t32 * t17 + t29 * t4 + t30 * t3) * g(1);
U = t1;
