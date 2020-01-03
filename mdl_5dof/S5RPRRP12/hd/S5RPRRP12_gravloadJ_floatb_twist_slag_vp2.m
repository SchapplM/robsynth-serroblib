% Calculate Gravitation load on the joints for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:16
% EndTime: 2019-12-31 18:56:18
% DurationCPUTime: 0.51s
% Computational Cost: add. (143->51), mult. (295->60), div. (0->0), fcn. (258->6), ass. (0->29)
t53 = mrSges(4,2) - m(5) * pkin(7) + m(6) * (-qJ(5) - pkin(7)) - mrSges(5,3) - mrSges(6,3);
t12 = sin(qJ(3));
t15 = cos(qJ(3));
t52 = t12 * mrSges(4,1) + t53 * t15;
t14 = cos(qJ(4));
t51 = -m(5) * pkin(3) - m(6) * (t14 * pkin(4) + pkin(3));
t49 = mrSges(5,1) + mrSges(6,1);
t41 = mrSges(5,2) + mrSges(6,2);
t13 = sin(qJ(1));
t16 = cos(qJ(1));
t38 = -g(1) * t13 + g(2) * t16;
t48 = m(6) * pkin(4);
t46 = -m(4) - m(5) - m(6);
t45 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t44 = m(3) - t46;
t43 = t51 * t12 + mrSges(2,2) - mrSges(3,3) - t52;
t11 = sin(qJ(4));
t40 = -t41 * t11 + t49 * t14 - t51;
t37 = t49 + t48;
t34 = t16 * pkin(1) + t13 * qJ(2);
t30 = t13 * t11;
t29 = t13 * t14;
t28 = t16 * t11;
t27 = t16 * t14;
t3 = t12 * t28 + t29;
t1 = -t12 * t30 + t27;
t4 = t12 * t27 - t30;
t2 = t12 * t29 + t28;
t5 = [(-t28 * t48 - m(3) * t34 + t46 * (t16 * pkin(6) + t34) - t49 * t2 - t41 * t1 + t45 * t16 + t43 * t13) * g(2) + (-t49 * t4 + t41 * t3 + (t11 * t48 + m(3) * pkin(1) + t46 * (-pkin(1) - pkin(6)) - t45) * t13 + (-t44 * qJ(2) + t43) * t16) * g(1), t38 * t44, (t40 * t12 + t52) * g(3) - t38 * ((-mrSges(4,1) - t40) * t15 + t53 * t12), (t37 * t11 + t41 * t14) * g(3) * t15 + (-t37 * t3 - t41 * t4) * g(2) + (-t37 * t1 + t41 * t2) * g(1), (-g(3) * t12 - t38 * t15) * m(6)];
taug = t5(:);
