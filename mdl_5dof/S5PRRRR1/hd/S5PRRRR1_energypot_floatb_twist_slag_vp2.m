% Calculate potential energy for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_energypot_floatb_twist_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:17
% EndTime: 2019-07-18 13:28:17
% DurationCPUTime: 0.16s
% Computational Cost: add. (87->44), mult. (99->31), div. (0->0), fcn. (67->8), ass. (0->20)
t24 = -m(3) - m(4);
t23 = -m(5) - m(6);
t13 = cos(qJ(3));
t22 = pkin(2) * t13;
t21 = mrSges(5,2) - mrSges(6,3);
t20 = -m(2) + t24;
t19 = pkin(1) + r_base(1);
t18 = qJ(1) + r_base(3);
t12 = cos(qJ(5));
t9 = sin(qJ(5));
t17 = mrSges(6,1) * t12 - mrSges(6,2) * t9 + mrSges(5,1);
t16 = t9 * mrSges(6,1) + t12 * mrSges(6,2) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t10 = sin(qJ(3));
t8 = qJ(3) + qJ(4);
t4 = sin(t8);
t5 = cos(t8);
t15 = -mrSges(4,1) * t13 + mrSges(4,2) * t10 - t17 * t5 + t21 * t4 - mrSges(3,1);
t14 = cos(qJ(2));
t11 = sin(qJ(2));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + t23 * (t11 * t22 + t18) + t20 * t18 + t16 * t14 + t15 * t11) * g(3) + (mrSges(4,1) * t10 + mrSges(4,2) * t13 - mrSges(1,2) - mrSges(2,2) + mrSges(3,3) + t21 * t5 + t23 * (-pkin(2) * t10 + r_base(2)) + t17 * t4 + (-m(1) + t20) * r_base(2)) * g(2) + (-mrSges(1,1) - mrSges(2,1) + (-m(1) - m(2)) * r_base(1) + t24 * t19 + t23 * (t14 * t22 + t19) - t16 * t11 + t15 * t14) * g(1);
U  = t1;
