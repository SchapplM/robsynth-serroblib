% Calculate Gravitation load on the joints for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:39
% EndTime: 2019-12-05 15:32:41
% DurationCPUTime: 0.31s
% Computational Cost: add. (145->36), mult. (207->44), div. (0->0), fcn. (171->8), ass. (0->22)
t41 = m(5) + m(6);
t35 = mrSges(5,2) + mrSges(6,2);
t34 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t23 = m(4) + t41;
t40 = t23 * pkin(2) + mrSges(3,1);
t13 = sin(qJ(4));
t15 = cos(qJ(4));
t38 = pkin(3) * t41 - t13 * t35 + t15 * t34 + mrSges(4,1);
t37 = -m(5) * pkin(6) + m(6) * (-qJ(5) - pkin(6)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t10 = sin(pkin(7));
t11 = cos(pkin(7));
t36 = g(1) * t11 + g(2) * t10;
t27 = t10 * t13;
t26 = t10 * t15;
t25 = t11 * t13;
t24 = t11 * t15;
t16 = cos(qJ(2));
t14 = sin(qJ(2));
t9 = qJ(2) + pkin(8);
t7 = cos(t9);
t6 = sin(t9);
t1 = [(-m(2) - m(3) - t23) * g(3), (t14 * mrSges(3,2) - t16 * t40 + t37 * t6 - t38 * t7) * g(3) + (mrSges(3,2) * t16 + t14 * t40 + t37 * t7 + t38 * t6) * t36, (-g(1) * t10 + g(2) * t11) * t23, (t13 * t34 + t15 * t35) * g(3) * t6 + (-t35 * (-t26 * t7 + t25) - t34 * (-t27 * t7 - t24)) * g(2) + (-t35 * (-t24 * t7 - t27) - t34 * (-t25 * t7 + t26)) * g(1), (g(3) * t7 - t36 * t6) * m(6)];
taug = t1(:);
