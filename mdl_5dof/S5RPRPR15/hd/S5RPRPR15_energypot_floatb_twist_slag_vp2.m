% Calculate potential energy for
% S5RPRPR15
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR15_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR15_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR15_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:33
% EndTime: 2019-12-31 18:36:33
% DurationCPUTime: 0.38s
% Computational Cost: add. (124->54), mult. (155->40), div. (0->0), fcn. (127->8), ass. (0->19)
t11 = sin(pkin(8));
t12 = cos(pkin(8));
t9 = pkin(8) + qJ(5);
t2 = sin(t9);
t3 = cos(t9);
t36 = mrSges(4,1) + m(5) * pkin(3) + t12 * mrSges(5,1) - t11 * mrSges(5,2) + m(6) * (pkin(4) * t12 + pkin(3)) + t3 * mrSges(6,1) - t2 * mrSges(6,2);
t35 = mrSges(4,2) - m(5) * qJ(4) + m(6) * (-pkin(7) - qJ(4)) - mrSges(5,3) - mrSges(6,3);
t14 = sin(qJ(3));
t16 = cos(qJ(3));
t34 = t36 * t14 + t35 * t16 - mrSges(2,2) + mrSges(3,3);
t33 = -m(1) - m(2);
t31 = -m(4) - m(5) - m(6);
t29 = -t2 * mrSges(6,1) - t12 * mrSges(5,2) - t3 * mrSges(6,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t11;
t10 = pkin(5) + r_base(3);
t15 = sin(qJ(1));
t27 = t15 * pkin(1) + r_base(2);
t17 = cos(qJ(1));
t24 = t17 * pkin(1) + t15 * qJ(2) + r_base(1);
t1 = (-m(1) * r_base(3) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t10 + t31 * (pkin(2) + t10) - t36 * t16 + t35 * t14) * g(3) + (-m(3) * t27 - mrSges(1,2) + t33 * r_base(2) + t31 * (t15 * pkin(6) + t27) + ((m(3) - t31) * qJ(2) + t34) * t17 + t29 * t15) * g(2) + (-m(3) * t24 - mrSges(1,1) + t33 * r_base(1) + t31 * (t17 * pkin(6) + t24) + t29 * t17 - t34 * t15) * g(1);
U = t1;
