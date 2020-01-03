% Calculate potential energy for
% S5RRPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:34
% EndTime: 2019-12-31 19:52:34
% DurationCPUTime: 0.30s
% Computational Cost: add. (132->51), mult. (106->37), div. (0->0), fcn. (70->6), ass. (0->20)
t32 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t31 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t30 = -m(1) - m(2);
t29 = -m(5) - m(6);
t28 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t12 = sin(qJ(4));
t14 = cos(qJ(4));
t27 = t32 * t12 + t31 * t14 - mrSges(3,2) + mrSges(4,3);
t26 = pkin(5) + r_base(3);
t13 = sin(qJ(1));
t25 = t13 * pkin(1) + r_base(2);
t15 = cos(qJ(1));
t24 = t15 * pkin(1) + r_base(1);
t8 = pkin(6) + t26;
t11 = qJ(1) + qJ(2);
t6 = sin(t11);
t23 = t6 * pkin(2) + t25;
t7 = cos(t11);
t20 = t7 * pkin(2) + t6 * qJ(3) + t24;
t1 = (-m(1) * r_base(3) - m(2) * t26 - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - m(4)) * t8 + t29 * (pkin(3) + t8) - t32 * t14 + t31 * t12) * g(3) + (-m(3) * t25 - m(4) * t23 - t13 * mrSges(2,1) - mrSges(2,2) * t15 - mrSges(1,2) + t30 * r_base(2) + t29 * (t6 * pkin(7) + t23) + ((m(4) - t29) * qJ(3) + t27) * t7 + t28 * t6) * g(2) + (-m(3) * t24 - m(4) * t20 - mrSges(2,1) * t15 + t13 * mrSges(2,2) - mrSges(1,1) + t30 * r_base(1) + t29 * (t7 * pkin(7) + t20) + t28 * t7 - t27 * t6) * g(1);
U = t1;
