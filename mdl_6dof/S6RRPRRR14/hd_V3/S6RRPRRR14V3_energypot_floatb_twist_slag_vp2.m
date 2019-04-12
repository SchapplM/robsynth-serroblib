% Calculate potential energy for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR14V3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(1,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR14V3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_energypot_floatb_twist_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14V3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:17
% EndTime: 2019-04-12 15:03:18
% DurationCPUTime: 0.26s
% Computational Cost: add. (124->53), mult. (241->61), div. (0->0), fcn. (240->10), ass. (0->29)
t17 = sin(qJ(2));
t18 = sin(qJ(1));
t39 = t17 * t18;
t21 = cos(qJ(4));
t38 = t17 * t21;
t23 = cos(qJ(1));
t37 = t17 * t23;
t22 = cos(qJ(2));
t36 = t18 * t22;
t16 = sin(qJ(4));
t35 = t23 * t16;
t34 = t23 * t21;
t33 = -mrSges(3,1) - mrSges(4,1);
t32 = mrSges(6,2) - mrSges(7,3);
t31 = qJ(3) * t17;
t30 = -m(1) - m(2) - m(3);
t29 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2);
t28 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t27 = -m(4) - m(5) - m(6) - m(7);
t14 = sin(qJ(6));
t19 = cos(qJ(6));
t26 = -t19 * mrSges(7,1) + t14 * mrSges(7,2) - mrSges(6,1);
t25 = -t14 * mrSges(7,1) - t19 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t24 = t17 * t28 + t22 * t33 - mrSges(2,1);
t20 = cos(qJ(5));
t15 = sin(qJ(5));
t10 = t18 * t16 + t22 * t34;
t8 = t21 * t36 - t35;
t1 = (-mrSges(1,3) - mrSges(2,3) + t32 * (t15 * t38 + t22 * t20) + t30 * r_base(3) + t26 * (-t22 * t15 + t20 * t38) - t28 * t22 + t27 * (-t22 * qJ(3) + r_base(3)) + (-mrSges(5,1) * t21 + t25 * t16 + t33) * t17) * g(3) + (-t8 * mrSges(5,1) - mrSges(1,2) + t26 * (t15 * t39 + t8 * t20) + t25 * (t16 * t36 + t34) + t32 * (t8 * t15 - t20 * t39) + t30 * r_base(2) - t29 * t23 + t27 * (t18 * t31 + r_base(2)) + t24 * t18) * g(2) + (-t10 * mrSges(5,1) - mrSges(1,1) + t26 * (t10 * t20 + t15 * t37) + t25 * (-t18 * t21 + t22 * t35) + t32 * (t10 * t15 - t20 * t37) + t30 * r_base(1) + t29 * t18 + t27 * (t23 * t31 + r_base(1)) + t24 * t23) * g(1);
U  = t1;
