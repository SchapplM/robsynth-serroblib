% Calculate Gravitation load on the joints for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:20
% EndTime: 2019-12-31 18:12:21
% DurationCPUTime: 0.33s
% Computational Cost: add. (177->61), mult. (216->60), div. (0->0), fcn. (166->6), ass. (0->27)
t44 = -mrSges(6,2) - mrSges(5,3);
t15 = sin(qJ(1));
t43 = g(2) * t15;
t35 = m(5) + m(6);
t16 = cos(qJ(1));
t30 = qJ(4) * t16;
t11 = pkin(7) + qJ(3);
t10 = cos(t11);
t7 = t10 * pkin(3);
t9 = sin(t11);
t42 = t16 * t7 + t9 * t30;
t23 = m(6) * (-pkin(3) - qJ(5)) - mrSges(6,3);
t41 = (m(5) * pkin(3) - mrSges(5,2) - t23) * t9 + t44 * t10;
t40 = g(1) * t16 + t43;
t39 = (-mrSges(4,1) + mrSges(5,2)) * t10 + (mrSges(4,2) + t44) * t9;
t38 = t10 * t35;
t13 = cos(pkin(7));
t37 = -mrSges(2,1) - m(3) * pkin(1) - t13 * mrSges(3,1) + sin(pkin(7)) * mrSges(3,2) + t39;
t14 = -pkin(6) - qJ(2);
t36 = mrSges(2,2) - m(3) * qJ(2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,1) - m(6) * (pkin(4) - t14) - mrSges(6,1);
t6 = t9 * qJ(4);
t34 = t7 + t6;
t8 = pkin(2) * t13 + pkin(1);
t28 = -t8 - t6;
t4 = t16 * t8;
t27 = -t15 * t14 + t4;
t1 = [(-m(4) * t27 - m(5) * (t27 + t42) - m(6) * (t4 + t42) + (-(m(6) * qJ(5) + mrSges(6,3)) * t10 + t37) * t16 + t36 * t15) * g(2) + (((m(4) + m(5)) * t14 + t36) * t16 + (m(4) * t8 - m(5) * (t28 - t7) - m(6) * t28 - t23 * t10 - t37) * t15) * g(1), (-g(1) * t15 + g(2) * t16) * (m(3) + m(4) + t35), t40 * (mrSges(4,1) * t9 + mrSges(4,2) * t10) + (t41 * t16 - t30 * t38) * g(1) + (-m(5) * t34 - m(6) * (qJ(5) * t10 + t34) - t10 * mrSges(6,3) + t39) * g(3) + (-qJ(4) * t38 + t41) * t43, (g(3) * t10 - t40 * t9) * t35, (-g(3) * t9 - t10 * t40) * m(6)];
taug = t1(:);
