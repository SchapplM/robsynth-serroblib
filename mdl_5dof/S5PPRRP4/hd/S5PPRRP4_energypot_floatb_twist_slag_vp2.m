% Calculate potential energy for
% S5PPRRP4
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPRRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:23
% EndTime: 2019-12-31 17:34:23
% DurationCPUTime: 0.30s
% Computational Cost: add. (118->51), mult. (153->38), div. (0->0), fcn. (141->6), ass. (0->21)
t35 = m(5) + m(6);
t34 = mrSges(5,2) + mrSges(6,2);
t33 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t32 = -m(1) - m(2);
t31 = -mrSges(2,1) - mrSges(3,1);
t30 = mrSges(2,2) - mrSges(3,3);
t29 = -m(4) - t35;
t16 = sin(qJ(4));
t17 = cos(qJ(4));
t28 = t35 * pkin(3) - t34 * t16 + t33 * t17 + mrSges(4,1);
t27 = m(5) * pkin(6) - m(6) * (-qJ(5) - pkin(6)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t26 = cos(qJ(3));
t25 = sin(qJ(3));
t24 = sin(pkin(7));
t13 = qJ(1) + r_base(3);
t14 = cos(pkin(7));
t23 = t14 * pkin(1) + t24 * qJ(2) + r_base(1);
t21 = t24 * pkin(1) - qJ(2) * t14 + r_base(2);
t2 = t14 * t25 - t24 * t26;
t1 = -t14 * t26 - t24 * t25;
t3 = (-m(1) * r_base(3) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + t34 * t17 + t29 * (-pkin(5) + t13) + t33 * t16 + (-m(2) - m(3)) * t13) * g(3) + (-m(3) * t21 - mrSges(1,2) + t32 * r_base(2) + t31 * t24 - t30 * t14 + t29 * (t24 * pkin(2) + t21) + t28 * t2 + t27 * t1) * g(2) + (-m(3) * t23 - mrSges(1,1) + t32 * r_base(1) + t30 * t24 + t31 * t14 + t29 * (t14 * pkin(2) + t23) - t27 * t2 + t28 * t1) * g(1);
U = t3;
