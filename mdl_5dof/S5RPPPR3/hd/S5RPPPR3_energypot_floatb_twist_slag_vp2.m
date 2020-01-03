% Calculate potential energy for
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPPR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR3_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:47
% EndTime: 2019-12-31 17:43:48
% DurationCPUTime: 0.38s
% Computational Cost: add. (156->60), mult. (143->50), div. (0->0), fcn. (113->8), ass. (0->26)
t44 = mrSges(4,2) - mrSges(5,3);
t15 = sin(pkin(8));
t16 = cos(pkin(8));
t43 = pkin(3) * t16 + qJ(4) * t15;
t42 = -m(6) * pkin(4) - mrSges(4,1) - mrSges(5,1);
t40 = -m(1) - m(2);
t39 = -m(5) - m(6);
t14 = qJ(1) + pkin(7);
t8 = sin(t14);
t38 = t43 * t8;
t17 = sin(qJ(5));
t19 = cos(qJ(5));
t25 = t15 * t17 + t16 * t19;
t26 = t15 * t19 - t16 * t17;
t37 = -t25 * mrSges(6,1) - t26 * mrSges(6,2) + t44 * t15 + t42 * t16 - mrSges(3,1);
t36 = mrSges(5,2) + mrSges(4,3) - mrSges(3,2) - mrSges(6,3);
t33 = pkin(5) + r_base(3);
t18 = sin(qJ(1));
t32 = t18 * pkin(1) + r_base(2);
t20 = cos(qJ(1));
t31 = t20 * pkin(1) + r_base(1);
t30 = t8 * pkin(2) + t32;
t10 = qJ(2) + t33;
t9 = cos(t14);
t23 = -qJ(3) * t9 + t30;
t1 = (-m(1) * r_base(3) - m(2) * t33 - t26 * mrSges(6,1) + t25 * mrSges(6,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t39 * (t15 * pkin(3) - qJ(4) * t16 + t10) - t44 * t16 + t42 * t15 + (-m(3) - m(4)) * t10) * g(3) + (-mrSges(1,2) - t18 * mrSges(2,1) - mrSges(2,2) * t20 - m(3) * t32 - m(4) * t23 - m(5) * (t23 + t38) - m(6) * (t30 + t38) + t40 * r_base(2) + (-m(6) * (pkin(6) - qJ(3)) + t36) * t9 + t37 * t8) * g(2) + (-m(3) * t31 - mrSges(2,1) * t20 + t18 * mrSges(2,2) - mrSges(1,1) + t40 * r_base(1) + (m(6) * pkin(6) - t36) * t8 + (-m(4) + t39) * (t9 * pkin(2) + t8 * qJ(3) + t31) + (t39 * t43 + t37) * t9) * g(1);
U = t1;
