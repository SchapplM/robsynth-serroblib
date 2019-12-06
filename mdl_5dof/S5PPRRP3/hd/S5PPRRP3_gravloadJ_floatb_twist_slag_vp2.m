% Calculate Gravitation load on the joints for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:34
% EndTime: 2019-12-05 15:10:35
% DurationCPUTime: 0.28s
% Computational Cost: add. (134->43), mult. (332->67), div. (0->0), fcn. (355->8), ass. (0->33)
t48 = mrSges(5,1) + mrSges(6,1);
t47 = -mrSges(5,2) + mrSges(6,3);
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t46 = t47 * t22 + t48 * t24 + mrSges(4,1);
t45 = -m(5) - m(6);
t27 = pkin(4) * t24 + qJ(5) * t22;
t43 = -m(6) * t27 - t46;
t42 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t18 = sin(pkin(8));
t41 = t18 * t22;
t40 = t18 * t24;
t25 = cos(qJ(3));
t39 = t18 * t25;
t19 = sin(pkin(7));
t23 = sin(qJ(3));
t38 = t19 * t23;
t37 = t19 * t25;
t21 = cos(pkin(7));
t36 = t21 * t23;
t35 = t21 * t25;
t33 = m(3) + m(4) - t45;
t31 = m(6) * pkin(4) + t48;
t30 = -m(6) * qJ(5) - t47;
t20 = cos(pkin(8));
t11 = t20 * t24 + t22 * t39;
t10 = t20 * t35 + t38;
t9 = -t20 * t36 + t37;
t8 = t20 * t37 - t36;
t7 = -t20 * t38 - t35;
t3 = t10 * t22 - t21 * t40;
t1 = -t19 * t40 + t8 * t22;
t2 = [(-m(2) - t33) * g(3), (-g(1) * t19 + g(2) * t21) * t33, (t45 * pkin(6) * t39 + (t42 * t25 + (m(5) * pkin(3) - m(6) * (-pkin(3) - t27) + t46) * t23) * t18) * g(3) + (t45 * (t7 * pkin(3) + t8 * pkin(6)) + t42 * t8 + t43 * t7) * g(2) + (t45 * (t9 * pkin(3) + t10 * pkin(6)) + t43 * t9 + t42 * t10) * g(1), (t30 * (-t20 * t22 + t24 * t39) + t31 * t11) * g(3) + (t30 * (t19 * t41 + t8 * t24) + t31 * t1) * g(2) + (t30 * (t10 * t24 + t21 * t41) + t31 * t3) * g(1), (-g(1) * t3 - g(2) * t1 - g(3) * t11) * m(6)];
taug = t2(:);
