% Calculate Gravitation load on the joints for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:52:53
% EndTime: 2019-03-09 01:52:55
% DurationCPUTime: 0.48s
% Computational Cost: add. (260->73), mult. (317->79), div. (0->0), fcn. (265->10), ass. (0->37)
t65 = -mrSges(5,2) - m(7) * (-pkin(8) - qJ(5)) + mrSges(7,3) + m(6) * qJ(5) + mrSges(6,3);
t17 = pkin(9) + qJ(4);
t10 = sin(t17);
t12 = cos(t17);
t64 = -t10 * mrSges(5,1) + t65 * t12;
t24 = sin(qJ(1));
t25 = cos(qJ(1));
t55 = -g(1) * t24 + g(2) * t25;
t18 = sin(pkin(10));
t20 = cos(pkin(10));
t62 = -m(7) * (t20 * pkin(5) + pkin(4)) - m(6) * pkin(4) - t20 * mrSges(6,1) + t18 * mrSges(6,2);
t60 = -m(5) - m(7);
t59 = m(6) + m(7);
t23 = -pkin(7) - qJ(3);
t19 = sin(pkin(9));
t51 = pkin(3) * t19;
t58 = t24 * t23 + t25 * t51;
t16 = pkin(10) + qJ(6);
t11 = cos(t16);
t9 = sin(t16);
t57 = t11 * mrSges(7,1) - t9 * mrSges(7,2) - t62;
t54 = -t19 * mrSges(4,1) - cos(pkin(9)) * mrSges(4,2) + mrSges(2,2) - mrSges(3,3) + t62 * t10 + t64;
t53 = t20 * mrSges(6,2) + mrSges(2,1) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + (m(7) * pkin(5) + mrSges(6,1)) * t18;
t47 = t24 * t9;
t46 = t25 * t9;
t45 = t24 * t11;
t44 = t25 * t11;
t43 = t25 * pkin(1) + t24 * qJ(2);
t42 = t24 * t51 + t43;
t41 = -m(4) - m(5) - t59;
t14 = t25 * qJ(2);
t39 = -t24 * pkin(1) + t14;
t4 = t10 * t44 - t47;
t3 = t10 * t46 + t45;
t2 = t10 * t45 + t46;
t1 = -t10 * t47 + t44;
t5 = [(-m(6) * t42 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (-m(3) - m(4)) * t43 + t60 * (-t25 * t23 + t42) + (-m(4) * qJ(3) + m(6) * t23 - t53) * t25 + t54 * t24) * g(2) + (-m(3) * t39 - m(4) * t14 - m(6) * (t14 + t58) - t4 * mrSges(7,1) + t3 * mrSges(7,2) + t60 * (t39 + t58) + (-m(4) * (-pkin(1) - qJ(3)) + m(6) * pkin(1) + t53) * t24 + t54 * t25) * g(1), t55 * (m(3) - t41) (g(1) * t25 + g(2) * t24) * t41 (t57 * t10 - t64) * g(3) + t55 * ((mrSges(5,1) + t57) * t12 + t65 * t10) (-t10 * g(3) - t55 * t12) * t59, -g(1) * (t1 * mrSges(7,1) - t2 * mrSges(7,2)) - g(2) * (t3 * mrSges(7,1) + t4 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t9 - mrSges(7,2) * t11) * t12];
taug  = t5(:);
