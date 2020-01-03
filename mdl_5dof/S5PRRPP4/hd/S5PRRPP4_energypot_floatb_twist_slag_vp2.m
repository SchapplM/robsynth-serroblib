% Calculate potential energy for
% S5PRRPP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRPP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:43
% EndTime: 2019-12-31 17:40:43
% DurationCPUTime: 0.34s
% Computational Cost: add. (146->55), mult. (127->42), div. (0->0), fcn. (91->6), ass. (0->22)
t17 = sin(qJ(3));
t18 = cos(qJ(3));
t40 = pkin(3) * t18 + qJ(4) * t17;
t39 = -m(6) * pkin(4) - mrSges(4,1) - mrSges(5,1) - mrSges(6,1);
t38 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t37 = -m(1) - m(2);
t36 = -m(5) - m(6);
t14 = pkin(7) + qJ(2);
t8 = sin(t14);
t35 = t40 * t8;
t33 = t38 * t17 + t39 * t18 - mrSges(3,1);
t32 = mrSges(5,2) + mrSges(4,3) - mrSges(3,2) - mrSges(6,3);
t15 = sin(pkin(7));
t29 = t15 * pkin(1) + r_base(2);
t16 = cos(pkin(7));
t28 = t16 * pkin(1) + r_base(1);
t27 = qJ(1) + r_base(3);
t26 = t8 * pkin(2) + t29;
t10 = pkin(5) + t27;
t9 = cos(t14);
t22 = -t9 * pkin(6) + t26;
t1 = (-m(1) * r_base(3) - m(2) * t27 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t36 * (t17 * pkin(3) - t18 * qJ(4) + t10) + (-m(3) - m(4)) * t10 - t38 * t18 + t39 * t17) * g(3) + (-mrSges(1,2) - mrSges(2,1) * t15 - mrSges(2,2) * t16 - m(3) * t29 - m(4) * t22 - m(5) * (t22 + t35) - m(6) * (t26 + t35) + t37 * r_base(2) + (-m(6) * (-pkin(6) + qJ(5)) + t32) * t9 + t33 * t8) * g(2) + (-m(3) * t28 - mrSges(2,1) * t16 + mrSges(2,2) * t15 - mrSges(1,1) + t37 * r_base(1) + (m(6) * qJ(5) - t32) * t8 + (-m(4) + t36) * (t9 * pkin(2) + t8 * pkin(6) + t28) + (t36 * t40 + t33) * t9) * g(1);
U = t1;
