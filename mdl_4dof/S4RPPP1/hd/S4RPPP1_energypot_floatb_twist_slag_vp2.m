% Calculate potential energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPPP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:18
% EndTime: 2019-03-08 18:26:19
% DurationCPUTime: 0.31s
% Computational Cost: add. (108->54), mult. (185->53), div. (0->0), fcn. (177->6), ass. (0->27)
t40 = -m(1) - m(2);
t39 = m(4) + m(5);
t16 = cos(pkin(6));
t17 = cos(pkin(4));
t19 = cos(qJ(1));
t32 = t17 * t19;
t14 = sin(pkin(6));
t18 = sin(qJ(1));
t33 = t14 * t18;
t3 = -t16 * t32 + t33;
t30 = t18 * t16;
t4 = t14 * t32 + t30;
t38 = t4 * pkin(2) + t3 * qJ(3);
t37 = mrSges(4,1) + mrSges(5,1) + mrSges(3,3);
t36 = mrSges(3,2) - mrSges(4,3) - mrSges(5,2);
t35 = -m(5) * pkin(3) - t37;
t34 = -m(5) * qJ(4) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t15 = sin(pkin(4));
t27 = qJ(2) * t15;
t26 = pkin(5) + r_base(3);
t25 = t18 * pkin(1) + r_base(2);
t24 = t17 * qJ(2) + t26;
t23 = t19 * pkin(1) + t18 * t27 + r_base(1);
t21 = -t19 * t27 + t25;
t6 = t16 * t19 - t17 * t33;
t5 = t14 * t19 + t17 * t30;
t1 = (-m(1) * r_base(3) - m(2) * t26 - m(3) * t24 - mrSges(1,3) - mrSges(2,3) - t39 * (t15 * t14 * pkin(2) + t24) + t35 * t17 + ((t39 * qJ(3) - t36) * t16 + t34 * t14) * t15) * g(3) + (-mrSges(1,2) - t18 * mrSges(2,1) - m(3) * t21 - m(4) * (t21 + t38) - m(5) * (t25 + t38) + t40 * r_base(2) + t34 * t4 + t36 * t3 + (-mrSges(2,2) + (-m(5) * (-pkin(3) - qJ(2)) + t37) * t15) * t19) * g(2) + (-m(3) * t23 - mrSges(2,1) * t19 - mrSges(1,1) + t40 * r_base(1) - t39 * (t6 * pkin(2) + qJ(3) * t5 + t23) + t34 * t6 + t36 * t5 + (t35 * t15 + mrSges(2,2)) * t18) * g(1);
U  = t1;
