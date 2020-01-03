% Calculate potential energy for
% S5RRRPP7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPP7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPP7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP7_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:52
% EndTime: 2019-12-31 21:03:53
% DurationCPUTime: 0.42s
% Computational Cost: add. (139->67), mult. (228->64), div. (0->0), fcn. (214->6), ass. (0->30)
t46 = -m(1) - m(2);
t45 = -m(5) - m(6);
t19 = sin(qJ(3));
t20 = sin(qJ(2));
t22 = cos(qJ(3));
t44 = (pkin(3) * t22 + qJ(4) * t19) * t20;
t43 = mrSges(2,2) - mrSges(3,3);
t42 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t41 = -mrSges(4,3) - mrSges(5,2) + mrSges(6,3);
t23 = cos(qJ(2));
t40 = -t23 * mrSges(3,1) - mrSges(2,1) + (m(6) * qJ(5) + mrSges(3,2)) * t20;
t39 = -m(6) * pkin(4) - mrSges(4,1) - mrSges(5,1) - mrSges(6,1);
t21 = sin(qJ(1));
t38 = t21 * t20;
t37 = t21 * t23;
t24 = cos(qJ(1));
t36 = t24 * t20;
t35 = t24 * t23;
t18 = pkin(5) + r_base(3);
t33 = t20 * pkin(2) + t18;
t32 = t24 * pkin(1) + t21 * pkin(6) + r_base(1);
t30 = t21 * pkin(1) - t24 * pkin(6) + r_base(2);
t29 = pkin(2) * t35 + pkin(7) * t36 + t32;
t28 = -t23 * pkin(7) + t33;
t27 = pkin(2) * t37 + pkin(7) * t38 + t30;
t6 = t21 * t19 + t22 * t35;
t5 = t19 * t35 - t21 * t22;
t4 = -t24 * t19 + t22 * t37;
t3 = t19 * t37 + t24 * t22;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - m(4) * t28 - m(5) * (t28 + t44) - m(6) * (t33 + t44) + (-m(2) - m(3)) * t18 + (-mrSges(3,2) - m(6) * (-pkin(7) + qJ(5)) - t41) * t23 + (t42 * t19 + t39 * t22 - mrSges(3,1)) * t20) * g(3) + (-m(3) * t30 - m(4) * t27 - mrSges(1,2) + t46 * r_base(2) + t45 * (t4 * pkin(3) + t3 * qJ(4) + t27) - t43 * t24 + t40 * t21 + t39 * t4 + t41 * t38 + t42 * t3) * g(2) + (-m(3) * t32 - m(4) * t29 - mrSges(1,1) + t46 * r_base(1) + t45 * (t6 * pkin(3) + t5 * qJ(4) + t29) + t40 * t24 + t43 * t21 + t39 * t6 + t42 * t5 + t41 * t36) * g(1);
U = t1;
