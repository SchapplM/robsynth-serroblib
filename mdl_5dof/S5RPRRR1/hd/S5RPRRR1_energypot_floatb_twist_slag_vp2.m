% Calculate potential energy for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(1,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_energypot_floatb_twist_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:24:50
% EndTime: 2019-07-18 13:24:50
% DurationCPUTime: 0.18s
% Computational Cost: add. (80->40), mult. (135->38), div. (0->0), fcn. (115->8), ass. (0->18)
t23 = -m(1) - m(2);
t10 = sin(qJ(1));
t13 = cos(qJ(3));
t22 = t10 * t13;
t14 = cos(qJ(1));
t21 = t13 * t14;
t20 = mrSges(5,2) - mrSges(6,3);
t19 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t18 = -m(3) - m(4) - m(5) - m(6);
t11 = cos(qJ(5));
t7 = sin(qJ(5));
t17 = -t11 * mrSges(6,1) + t7 * mrSges(6,2) - mrSges(5,1);
t16 = t7 * mrSges(6,1) + t11 * mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t9 = sin(qJ(3));
t15 = -mrSges(4,1) * t13 - t16 * t9 - mrSges(2,1) - mrSges(3,1);
t12 = cos(qJ(4));
t8 = sin(qJ(4));
t1 = (-mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t16 * t13 + (t17 * t12 + t20 * t8 - mrSges(4,1)) * t9 + (t18 + t23) * r_base(3)) * g(3) + (-mrSges(1,2) + t23 * r_base(2) + t20 * (t12 * t14 + t8 * t22) + t17 * (t12 * t22 - t14 * t8) - t19 * t14 + t18 * (-qJ(2) * t14 + r_base(2)) + t15 * t10) * g(2) + (-mrSges(1,1) + t23 * r_base(1) + t20 * (-t10 * t12 + t8 * t21) + t17 * (t10 * t8 + t12 * t21) + t19 * t10 + t18 * (qJ(2) * t10 + r_base(1)) + t15 * t14) * g(1);
U  = t1;
