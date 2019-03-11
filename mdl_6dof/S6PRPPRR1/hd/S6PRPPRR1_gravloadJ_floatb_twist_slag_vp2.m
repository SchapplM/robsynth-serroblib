% Calculate Gravitation load on the joints for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:13:59
% EndTime: 2019-03-08 19:14:01
% DurationCPUTime: 0.64s
% Computational Cost: add. (506->82), mult. (1082->123), div. (0->0), fcn. (1304->14), ass. (0->45)
t44 = sin(qJ(6));
t46 = cos(qJ(6));
t84 = -m(7) * pkin(5) - t46 * mrSges(7,1) + t44 * mrSges(7,2) - mrSges(6,1);
t83 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t85 = -m(6) - m(7);
t82 = m(5) - t85;
t64 = m(4) + t82;
t37 = sin(pkin(11));
t45 = sin(qJ(2));
t47 = cos(qJ(2));
t68 = cos(pkin(11));
t24 = -t47 * t37 - t45 * t68;
t38 = sin(pkin(10));
t41 = cos(pkin(10));
t42 = cos(pkin(6));
t72 = t42 * t47;
t87 = -t38 * t45 + t41 * t72;
t86 = -m(4) - m(5);
t55 = -t45 * t37 + t47 * t68;
t69 = t24 * t42;
t14 = t38 * t69 + t41 * t55;
t9 = -t38 * t55 + t41 * t69;
t35 = pkin(12) + qJ(5);
t33 = sin(t35);
t34 = cos(t35);
t40 = cos(pkin(12));
t81 = -m(5) * pkin(3) - t40 * mrSges(5,1) + sin(pkin(12)) * mrSges(5,2) - mrSges(4,1) + t84 * t34 + t83 * t33;
t80 = m(5) * qJ(4) + mrSges(7,1) * t44 + mrSges(7,2) * t46 - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t39 = sin(pkin(6));
t78 = t38 * t39;
t76 = t39 * t41;
t73 = t42 * t45;
t30 = t39 * t47 * pkin(2);
t60 = t87 * pkin(2);
t50 = t42 * t55;
t43 = -pkin(8) - qJ(4);
t32 = pkin(4) * t40 + pkin(3);
t22 = t24 * t39;
t21 = t55 * t39;
t16 = -t22 * t34 + t33 * t42;
t13 = t24 * t41 - t38 * t50;
t10 = t38 * t24 + t41 * t50;
t4 = t14 * t34 + t33 * t78;
t2 = -t33 * t76 - t34 * t9;
t1 = [(-m(2) - m(3) - t64) * g(3) (-(mrSges(3,1) * t47 - mrSges(3,2) * t45) * t39 + t85 * (t21 * t32 + t22 * t43 + t30) + t86 * t30 + t80 * t22 + t81 * t21) * g(3) + (-t87 * mrSges(3,1) - (-t38 * t47 - t41 * t73) * mrSges(3,2) + t86 * t60 + t85 * (t10 * t32 + t9 * t43 + t60) + t80 * t9 + t81 * t10) * g(2) + (-(t38 * t73 - t41 * t47) * mrSges(3,2) + t85 * (t13 * t32 - t14 * t43) + t81 * t13 - t80 * t14 + (-t64 * pkin(2) - mrSges(3,1)) * (-t38 * t72 - t41 * t45)) * g(1) (-g(3) * t42 + (-g(1) * t38 + g(2) * t41) * t39) * t64, t82 * (g(1) * t13 + g(2) * t10 + g(3) * t21) (t83 * t16 + t84 * (t22 * t33 + t34 * t42)) * g(3) + (t83 * t2 + t84 * (t33 * t9 - t34 * t76)) * g(2) + (t83 * t4 + t84 * (-t14 * t33 + t34 * t78)) * g(1), -g(1) * ((-t13 * t46 - t4 * t44) * mrSges(7,1) + (t13 * t44 - t4 * t46) * mrSges(7,2)) - g(2) * ((-t10 * t46 - t2 * t44) * mrSges(7,1) + (t10 * t44 - t2 * t46) * mrSges(7,2)) - g(3) * ((-t16 * t44 - t21 * t46) * mrSges(7,1) + (-t16 * t46 + t21 * t44) * mrSges(7,2))];
taug  = t1(:);
