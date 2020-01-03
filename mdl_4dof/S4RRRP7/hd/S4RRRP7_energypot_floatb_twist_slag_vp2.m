% Calculate potential energy for
% S4RRRP7
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRP7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:06
% EndTime: 2019-12-31 17:20:06
% DurationCPUTime: 0.31s
% Computational Cost: add. (94->53), mult. (148->52), div. (0->0), fcn. (130->6), ass. (0->21)
t37 = -m(1) - m(2);
t36 = -m(4) - m(5);
t16 = sin(qJ(2));
t19 = cos(qJ(2));
t35 = -t19 * mrSges(3,1) + t16 * mrSges(3,2) - mrSges(2,1);
t34 = mrSges(2,2) - mrSges(3,3);
t33 = -mrSges(4,3) - mrSges(5,2);
t32 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t31 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3);
t17 = sin(qJ(1));
t30 = t17 * t16;
t29 = t17 * t19;
t20 = cos(qJ(1));
t28 = t20 * t16;
t27 = t20 * t19;
t14 = pkin(4) + r_base(3);
t26 = t20 * pkin(1) + t17 * pkin(5) + r_base(1);
t24 = t17 * pkin(1) - pkin(5) * t20 + r_base(2);
t18 = cos(qJ(3));
t15 = sin(qJ(3));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + t36 * (t16 * pkin(2) - pkin(6) * t19 + t14) + (-m(2) - m(3)) * t14 + (-mrSges(3,2) - t33) * t19 + (t31 * t15 + t32 * t18 - mrSges(3,1)) * t16) * g(3) + (-m(3) * t24 - mrSges(1,2) + t37 * r_base(2) + t33 * t30 + t36 * (pkin(2) * t29 + pkin(6) * t30 + t24) - t34 * t20 + t32 * (-t15 * t20 + t18 * t29) + t35 * t17 + t31 * (t15 * t29 + t18 * t20)) * g(2) + (-m(3) * t26 - mrSges(1,1) + t37 * r_base(1) + t36 * (pkin(2) * t27 + pkin(6) * t28 + t26) + t32 * (t17 * t15 + t18 * t27) + t31 * (t15 * t27 - t17 * t18) + t33 * t28 + t35 * t20 + t34 * t17) * g(1);
U = t1;
