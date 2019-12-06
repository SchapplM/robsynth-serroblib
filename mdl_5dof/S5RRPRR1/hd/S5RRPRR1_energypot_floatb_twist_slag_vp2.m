% Calculate potential energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energypot_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:23:46
% EndTime: 2019-12-05 18:23:46
% DurationCPUTime: 0.36s
% Computational Cost: add. (108->50), mult. (127->44), div. (0->0), fcn. (95->8), ass. (0->24)
t37 = -m(4) * pkin(1) - mrSges(3,1) - mrSges(4,1);
t36 = mrSges(3,2) + mrSges(4,2);
t35 = -m(6) * pkin(4) + mrSges(5,2) - mrSges(6,3);
t33 = -m(5) - m(6);
t32 = -m(1) - m(2) - m(3) - m(4);
t30 = m(4) * qJ(3) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t10 = sin(qJ(2));
t13 = cos(qJ(2));
t7 = qJ(2) + qJ(4);
t5 = sin(t7);
t6 = cos(t7);
t29 = -t6 * mrSges(5,1) + t36 * t10 + t37 * t13 + t35 * t5 - mrSges(2,1);
t11 = sin(qJ(1));
t9 = sin(qJ(5));
t27 = t11 * t9;
t14 = cos(qJ(1));
t26 = t14 * t9;
t12 = cos(qJ(5));
t24 = t11 * t12;
t23 = t12 * t14;
t15 = pkin(2) + pkin(1);
t22 = t13 * t15;
t8 = -pkin(3) - qJ(3);
t1 = (-mrSges(1,3) - mrSges(2,3) + t33 * (t10 * t15 + r_base(3)) - t35 * t6 + (-t12 * mrSges(6,1) + t9 * mrSges(6,2) - mrSges(5,1)) * t5 - t36 * t13 + t37 * t10 + t32 * r_base(3)) * g(3) + (-mrSges(1,2) - (t6 * t24 - t26) * mrSges(6,1) - (-t6 * t27 - t23) * mrSges(6,2) + t33 * (t11 * t22 + t14 * t8 + r_base(2)) + t32 * r_base(2) + t30 * t14 + t29 * t11) * g(2) + (-mrSges(1,1) - (t6 * t23 + t27) * mrSges(6,1) - (-t6 * t26 + t24) * mrSges(6,2) + t33 * (-t11 * t8 + t14 * t22 + r_base(1)) + t32 * r_base(1) - t30 * t11 + t29 * t14) * g(1);
U = t1;
