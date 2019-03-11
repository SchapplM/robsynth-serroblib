% Calculate potential energy for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPPR8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:41
% EndTime: 2019-03-09 02:58:41
% DurationCPUTime: 0.52s
% Computational Cost: add. (142->73), mult. (193->54), div. (0->0), fcn. (155->6), ass. (0->26)
t45 = m(6) + m(7);
t44 = -mrSges(4,1) - mrSges(5,1);
t43 = -m(1) - m(2);
t41 = mrSges(7,3) - mrSges(6,2);
t40 = m(5) + t45;
t13 = sin(qJ(6));
t16 = cos(qJ(6));
t39 = -m(7) * pkin(5) - t16 * mrSges(7,1) + t13 * mrSges(7,2) - mrSges(6,1) - mrSges(5,3);
t14 = sin(qJ(3));
t17 = cos(qJ(3));
t38 = -mrSges(4,2) * t17 + t44 * t14 + mrSges(2,2) - mrSges(3,3);
t37 = -m(7) * pkin(8) - t41;
t36 = mrSges(7,1) * t13 + mrSges(7,2) * t16 + t45 * qJ(5) - mrSges(2,1) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t35 = -pkin(3) - pkin(4);
t15 = sin(qJ(1));
t33 = t14 * t15;
t12 = pkin(6) + r_base(3);
t32 = t15 * pkin(1) + r_base(2);
t31 = pkin(2) + t12;
t30 = t15 * pkin(7) + t32;
t18 = cos(qJ(1));
t29 = t18 * pkin(1) + t15 * qJ(2) + r_base(1);
t27 = t18 * pkin(7) + t29;
t26 = t17 * pkin(3) + t14 * qJ(4) + t31;
t23 = pkin(3) * t33 + t27;
t1 = (-m(1) * r_base(3) - m(4) * t31 - m(5) * t26 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) - t45 * (t17 * pkin(4) + t26) + (-m(2) - m(3)) * t12 + (t37 + t44) * t17 + (mrSges(4,2) + t39) * t14) * g(3) + (-m(3) * t32 - m(4) * t30 - mrSges(1,2) + t43 * r_base(2) - t40 * (t18 * t17 * qJ(4) + t30) + (t39 * t17 + (m(5) * pkin(3) - m(6) * t35 - m(7) * (-pkin(8) + t35) + t41) * t14 + (m(3) + m(4) + t40) * qJ(2) - t38) * t18 + t36 * t15) * g(2) + (-m(3) * t29 - m(4) * t27 - m(5) * t23 - mrSges(1,1) + t43 * r_base(1) - t45 * (pkin(4) * t33 + t23) + t36 * t18 + (t37 * t14 + (t40 * qJ(4) - t39) * t17 + t38) * t15) * g(1);
U  = t1;
