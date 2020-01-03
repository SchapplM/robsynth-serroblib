% Calculate potential energy for
% S5RRPRP7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRP7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP7_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:51
% EndTime: 2019-12-31 19:59:51
% DurationCPUTime: 0.41s
% Computational Cost: add. (167->65), mult. (179->60), div. (0->0), fcn. (155->8), ass. (0->31)
t49 = -m(2) - m(3);
t48 = -m(5) - m(6);
t47 = -mrSges(5,3) - mrSges(6,2);
t46 = -m(1) + t49;
t18 = qJ(2) + pkin(8);
t15 = sin(t18);
t16 = cos(t18);
t22 = sin(qJ(2));
t25 = cos(qJ(2));
t45 = -m(3) * pkin(1) - t25 * mrSges(3,1) - t16 * mrSges(4,1) + t22 * mrSges(3,2) + t15 * mrSges(4,2) - mrSges(2,1);
t44 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t43 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t42 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t41 = pkin(3) * t16;
t23 = sin(qJ(1));
t40 = t23 * t15;
t21 = sin(qJ(4));
t39 = t23 * t21;
t24 = cos(qJ(4));
t38 = t23 * t24;
t26 = cos(qJ(1));
t37 = t26 * t15;
t36 = t26 * t21;
t35 = t26 * t24;
t19 = pkin(5) + r_base(3);
t34 = t22 * pkin(2) + t19;
t14 = t25 * pkin(2) + pkin(1);
t20 = -qJ(3) - pkin(6);
t33 = t23 * t14 + t26 * t20 + r_base(2);
t30 = t26 * t14 - t23 * t20 + r_base(1);
t1 = (-m(1) * r_base(3) - m(4) * t34 - t22 * mrSges(3,1) - t25 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t48 * (t15 * pkin(3) - t16 * pkin(7) + t34) + t49 * t19 + (-mrSges(4,2) - t47) * t16 + (t43 * t21 + t44 * t24 - mrSges(4,1)) * t15) * g(3) + (-m(4) * t33 - mrSges(1,2) + t47 * t40 + t48 * (pkin(7) * t40 + t23 * t41 + t33) + t44 * (t16 * t38 - t36) + t43 * (t16 * t39 + t35) + t46 * r_base(2) + t42 * t26 + t45 * t23) * g(2) + (-m(4) * t30 - mrSges(1,1) + t48 * (pkin(7) * t37 + t26 * t41 + t30) + t44 * (t16 * t35 + t39) + t47 * t37 + t43 * (t16 * t36 - t38) + t46 * r_base(1) + t45 * t26 - t42 * t23) * g(1);
U = t1;
