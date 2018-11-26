% Calculate Gravitation load on the joints for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:47
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:47:24
% EndTime: 2018-11-23 15:47:25
% DurationCPUTime: 0.52s
% Computational Cost: add. (276->77), mult. (383->96), div. (0->0), fcn. (342->8), ass. (0->45)
t65 = mrSges(6,1) + mrSges(7,1);
t64 = -mrSges(6,2) + mrSges(7,3);
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t63 = t64 * t23 + t65 * t25;
t60 = -mrSges(6,3) - mrSges(7,2);
t61 = m(6) + m(7);
t62 = t61 * pkin(8) - t60;
t33 = pkin(5) * t25 + qJ(6) * t23;
t59 = m(6) * pkin(4) - m(7) * (-pkin(4) - t33) + t63;
t58 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t20 = sin(pkin(9));
t19 = pkin(9) + qJ(4);
t14 = sin(t19);
t15 = cos(t19);
t36 = t14 * mrSges(5,1) + t15 * mrSges(5,2);
t57 = mrSges(2,2) - mrSges(3,3) - t20 * mrSges(4,1) - cos(pkin(9)) * mrSges(4,2) - t36;
t56 = m(7) * pkin(5) + t65;
t55 = m(7) * qJ(6) + t64;
t53 = pkin(3) * t20;
t24 = sin(qJ(1));
t52 = g(1) * t24;
t26 = cos(qJ(1));
t51 = g(2) * t26;
t50 = g(3) * t15;
t49 = t14 * t24;
t48 = t15 * t24;
t47 = t15 * t26;
t46 = t24 * t23;
t45 = t24 * t25;
t44 = t26 * t23;
t43 = t26 * t25;
t42 = t26 * pkin(1) + t24 * qJ(2);
t41 = -m(4) - m(5) - t61;
t17 = t26 * qJ(2);
t40 = -t24 * pkin(1) + t17;
t37 = mrSges(5,1) * t15 - mrSges(5,2) * t14;
t22 = -pkin(7) - qJ(3);
t32 = t24 * t22 + t26 * t53 + t40;
t31 = -t26 * t22 + t24 * t53 + t42;
t4 = t14 * t43 - t46;
t3 = t14 * t44 + t45;
t2 = t14 * t45 + t44;
t1 = t14 * t46 - t43;
t5 = [(-m(5) * t31 - t60 * t48 + (-m(3) - m(4)) * t42 - t61 * (pkin(4) * t49 - pkin(8) * t48 + t31) - t56 * t2 - t55 * t1 + (-m(4) * qJ(3) - t58) * t26 + t57 * t24) * g(2) + (-m(3) * t40 - m(4) * t17 - m(5) * t32 - t60 * t47 - t61 * (t26 * t14 * pkin(4) - pkin(8) * t47 + t32) - t56 * t4 - t55 * t3 + t57 * t26 + (-m(4) * (-pkin(1) - qJ(3)) + t58) * t24) * g(1) (-t52 + t51) * (m(3) - t41) (g(1) * t26 + g(2) * t24) * t41, -t37 * t52 + (-t61 * (pkin(4) * t48 + pkin(8) * t49) + ((-m(7) * t33 - t63) * t15 + t60 * t14) * t24) * g(1) + (t62 * t14 + t59 * t15 + t37) * t51 + (t59 * t14 - t62 * t15 + t36) * g(3) (t56 * t23 - t55 * t25) * t50 + (-t56 * t3 + t55 * t4) * g(2) + (t56 * t1 - t55 * t2) * g(1) (-g(1) * t1 + g(2) * t3 - t23 * t50) * m(7)];
taug  = t5(:);
