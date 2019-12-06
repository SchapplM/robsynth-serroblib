% Calculate potential energy for
% S5RPPRR1
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:56
% EndTime: 2019-12-05 17:37:56
% DurationCPUTime: 0.28s
% Computational Cost: add. (98->56), mult. (103->39), div. (0->0), fcn. (67->6), ass. (0->20)
t27 = -m(1) - m(2);
t26 = -m(5) - m(6);
t9 = qJ(4) + qJ(5);
t1 = sin(t9);
t10 = sin(qJ(4));
t12 = cos(qJ(4));
t2 = cos(t9);
t25 = -t1 * mrSges(6,1) - mrSges(5,2) * t12 - t2 * mrSges(6,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t10;
t24 = mrSges(4,2) + mrSges(3,3) - mrSges(2,2) - mrSges(5,3) - mrSges(6,3);
t8 = pkin(5) + r_base(3);
t11 = sin(qJ(1));
t23 = t11 * pkin(1) + r_base(2);
t22 = pkin(2) + t8;
t13 = cos(qJ(1));
t20 = t13 * pkin(1) + t11 * qJ(2) + r_base(1);
t19 = pkin(3) + t22;
t16 = -t13 * qJ(2) + t23;
t14 = -pkin(7) - pkin(6);
t3 = t11 * qJ(3);
t4 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - mrSges(3,1) - m(4) * t22 - mrSges(4,1) - m(5) * t19 - t12 * mrSges(5,1) + t10 * mrSges(5,2) - m(6) * (pkin(4) * t12 + t19) - t2 * mrSges(6,1) + t1 * mrSges(6,2) + (-m(2) - m(3)) * t8) * g(3) + (-mrSges(1,2) - m(3) * t16 - m(4) * (t16 + t3) + t27 * r_base(2) + t26 * (t3 + t23) + (-m(5) * (pkin(6) - qJ(2)) - m(6) * (-qJ(2) - t14) + t24) * t13 + t25 * t11) * g(2) + (-m(3) * t20 - mrSges(1,1) + t27 * r_base(1) + (-m(4) + t26) * (t13 * qJ(3) + t20) + t25 * t13 + (m(5) * pkin(6) - m(6) * t14 - t24) * t11) * g(1);
U = t4;
