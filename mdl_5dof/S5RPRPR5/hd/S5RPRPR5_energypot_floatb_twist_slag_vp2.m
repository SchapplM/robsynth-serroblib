% Calculate potential energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:41:18
% EndTime: 2020-01-03 11:41:19
% DurationCPUTime: 0.58s
% Computational Cost: add. (170->57), mult. (192->43), div. (0->0), fcn. (168->10), ass. (0->23)
t17 = sin(qJ(3));
t19 = cos(qJ(3));
t12 = qJ(3) + pkin(9);
t8 = qJ(5) + t12;
t3 = sin(t8);
t4 = cos(t8);
t5 = t19 * pkin(3) + pkin(2);
t6 = sin(t12);
t7 = cos(t12);
t40 = -mrSges(3,1) - m(4) * pkin(2) - mrSges(4,1) * t19 + mrSges(4,2) * t17 - m(6) * (pkin(4) * t7 + t5) - mrSges(6,1) * t4 + mrSges(6,2) * t3 - m(5) * t5 - mrSges(5,1) * t7 + mrSges(5,2) * t6;
t16 = -qJ(4) - pkin(6);
t39 = mrSges(3,2) - m(4) * pkin(6) + m(5) * t16 + m(6) * (-pkin(7) + t16) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t38 = m(3) + m(4);
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t37 = t39 * t14 + t40 * t15 - mrSges(2,1);
t35 = -pkin(3) * t17 - qJ(2);
t33 = m(5) + m(6) + t38;
t32 = -m(2) - t33;
t31 = -m(5) * t35 + t6 * mrSges(5,1) + t7 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3) - m(6) * (-pkin(4) * t6 + t35) + t3 * mrSges(6,1) + t4 * mrSges(6,2) + t17 * mrSges(4,1) + t19 * mrSges(4,2) + t38 * qJ(2);
t20 = cos(qJ(1));
t18 = sin(qJ(1));
t1 = (-mrSges(1,3) + (t33 * pkin(1) - t37) * t20 + (-m(1) + t32) * r_base(3) + t31 * t18) * g(3) + (-mrSges(1,2) + (-m(1) - m(2)) * r_base(2) - t33 * (t18 * pkin(1) + r_base(2)) + t31 * t20 + t37 * t18) * g(2) + (-m(1) * r_base(1) - mrSges(1,1) - mrSges(2,3) + t32 * (pkin(5) + r_base(1)) - t39 * t15 + t40 * t14) * g(1);
U = t1;
