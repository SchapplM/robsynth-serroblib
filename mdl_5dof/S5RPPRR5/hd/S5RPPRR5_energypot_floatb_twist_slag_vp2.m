% Calculate potential energy for
% S5RPPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:26
% EndTime: 2019-12-31 17:56:26
% DurationCPUTime: 0.27s
% Computational Cost: add. (153->56), mult. (120->45), div. (0->0), fcn. (98->8), ass. (0->24)
t37 = -m(1) - m(2);
t36 = -m(5) - m(6);
t35 = -mrSges(3,1) - mrSges(4,1);
t16 = sin(qJ(5));
t18 = cos(qJ(5));
t34 = m(6) * pkin(4) + t18 * mrSges(6,1) - t16 * mrSges(6,2) + mrSges(5,1);
t33 = mrSges(3,2) - mrSges(4,3);
t32 = m(6) * pkin(7) - mrSges(5,2) + mrSges(6,3);
t31 = cos(qJ(4));
t30 = sin(qJ(4));
t29 = pkin(5) + r_base(3);
t28 = qJ(1) + pkin(8);
t17 = sin(qJ(1));
t27 = t17 * pkin(1) + r_base(2);
t19 = cos(qJ(1));
t26 = t19 * pkin(1) + r_base(1);
t13 = qJ(2) + t29;
t25 = sin(t28);
t12 = cos(t28);
t24 = t12 * pkin(2) + t25 * qJ(3) + t26;
t22 = t25 * pkin(2) - t12 * qJ(3) + t27;
t2 = t12 * t30 - t25 * t31;
t1 = -t12 * t31 - t25 * t30;
t3 = (-m(1) * r_base(3) - m(2) * t29 + mrSges(6,1) * t16 + t18 * mrSges(6,2) - mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + mrSges(5,3) + (-m(3) - m(4)) * t13 + t36 * (-pkin(6) + t13)) * g(3) + (-m(3) * t27 - m(4) * t22 - t17 * mrSges(2,1) - mrSges(2,2) * t19 - mrSges(1,2) + t37 * r_base(2) + t35 * t25 + t36 * (t25 * pkin(3) + t22) + t34 * t2 - t33 * t12 + t32 * t1) * g(2) + (-m(3) * t26 - m(4) * t24 - mrSges(2,1) * t19 + t17 * mrSges(2,2) - mrSges(1,1) + t37 * r_base(1) + t33 * t25 + t36 * (t12 * pkin(3) + t24) - t32 * t2 + t35 * t12 + t34 * t1) * g(1);
U = t3;
