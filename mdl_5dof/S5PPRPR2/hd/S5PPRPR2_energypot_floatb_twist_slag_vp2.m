% Calculate potential energy for
% S5PPRPR2
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRPR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPRPR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:02:53
% EndTime: 2019-12-05 15:02:54
% DurationCPUTime: 0.45s
% Computational Cost: add. (151->72), mult. (153->65), div. (0->0), fcn. (121->8), ass. (0->33)
t45 = -mrSges(4,1) + mrSges(5,2);
t44 = mrSges(4,2) - mrSges(5,3);
t43 = -m(2) - m(3);
t42 = m(5) + m(6);
t18 = cos(pkin(7));
t14 = pkin(8) + qJ(3);
t10 = sin(t14);
t31 = qJ(4) * t10;
t11 = cos(t14);
t34 = t18 * t11;
t41 = pkin(3) * t34 + t18 * t31;
t40 = -m(1) + t43;
t15 = sin(pkin(8));
t17 = cos(pkin(8));
t39 = -m(3) * pkin(1) - t17 * mrSges(3,1) + t15 * mrSges(3,2) + t44 * t10 + t45 * t11 - mrSges(2,1);
t38 = -m(3) * qJ(2) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t16 = sin(pkin(7));
t37 = t16 * t11;
t20 = sin(qJ(5));
t36 = t16 * t20;
t21 = cos(qJ(5));
t35 = t16 * t21;
t33 = t18 * t20;
t32 = t18 * t21;
t9 = pkin(2) * t17 + pkin(1);
t30 = t18 * t9 + r_base(1);
t13 = qJ(1) + r_base(3);
t19 = -pkin(5) - qJ(2);
t29 = t16 * t9 + t18 * t19 + r_base(2);
t28 = t15 * pkin(2) + t13;
t26 = pkin(3) * t37 + t16 * t31 + t29;
t23 = -t16 * t19 + t30;
t1 = (-m(1) * r_base(3) - m(4) * t28 - mrSges(3,1) * t15 - t17 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) - t42 * (t10 * pkin(3) + t28) + t43 * t13 + (t20 * mrSges(6,1) + t21 * mrSges(6,2) + t42 * qJ(4) - t44) * t11 + (-m(6) * pkin(6) - mrSges(6,3) + t45) * t10) * g(3) + (-mrSges(1,2) - m(4) * t29 - m(5) * t26 - m(6) * (pkin(6) * t37 + t26) - (t10 * t36 - t32) * mrSges(6,1) - (t10 * t35 + t33) * mrSges(6,2) - mrSges(6,3) * t37 + t40 * r_base(2) + (m(6) * pkin(4) - t38) * t18 + t39 * t16) * g(2) + (-mrSges(1,1) - m(4) * t23 - m(5) * (t23 + t41) - m(6) * (pkin(6) * t34 + t30 + t41) - (t10 * t33 + t35) * mrSges(6,1) - (t10 * t32 - t36) * mrSges(6,2) - mrSges(6,3) * t34 + t40 * r_base(1) + t39 * t18 + (-m(6) * (pkin(4) - t19) + t38) * t16) * g(1);
U = t1;
