% Calculate Gravitation load on the joints for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:32
% EndTime: 2019-12-31 19:14:34
% DurationCPUTime: 0.54s
% Computational Cost: add. (184->69), mult. (310->80), div. (0->0), fcn. (277->8), ass. (0->43)
t70 = mrSges(4,2) - m(5) * pkin(7) + m(6) * (-pkin(8) - pkin(7)) - mrSges(5,3) - mrSges(6,3);
t22 = sin(qJ(3));
t25 = cos(qJ(3));
t69 = t22 * mrSges(4,1) + t70 * t25;
t24 = cos(qJ(4));
t68 = -m(5) * pkin(3) - m(6) * (t24 * pkin(4) + pkin(3));
t65 = m(6) * pkin(4);
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t66 = g(1) * t23 - g(2) * t26;
t21 = sin(qJ(4));
t64 = t21 * t65;
t63 = -m(4) - m(5) - m(6);
t62 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t61 = m(3) - t63;
t60 = t68 * t22 + mrSges(2,2) - mrSges(3,3) - t69;
t58 = -mrSges(5,1) - t65;
t20 = qJ(4) + qJ(5);
t14 = sin(t20);
t15 = cos(t20);
t57 = t24 * mrSges(5,1) + t15 * mrSges(6,1) - t21 * mrSges(5,2) - t14 * mrSges(6,2) - t68;
t42 = t26 * t15;
t47 = t23 * t14;
t5 = -t22 * t47 + t42;
t43 = t26 * t14;
t46 = t23 * t15;
t6 = t22 * t46 + t43;
t53 = t5 * mrSges(6,1) - t6 * mrSges(6,2);
t7 = t22 * t43 + t46;
t8 = t22 * t42 - t47;
t52 = t7 * mrSges(6,1) + t8 * mrSges(6,2);
t48 = g(3) * t25;
t45 = t23 * t21;
t44 = t23 * t24;
t41 = t26 * t21;
t40 = t26 * t24;
t39 = t26 * pkin(1) + t23 * qJ(2);
t32 = -mrSges(6,1) * t14 - mrSges(6,2) * t15;
t11 = t22 * t41 + t44;
t9 = -t22 * t45 + t40;
t12 = t22 * t40 - t45;
t10 = t22 * t44 + t41;
t1 = [(-t41 * t65 - m(3) * t39 - t10 * mrSges(5,1) - t6 * mrSges(6,1) - t9 * mrSges(5,2) - t5 * mrSges(6,2) + t63 * (t26 * pkin(6) + t39) + t62 * t26 + t60 * t23) * g(2) + (-t12 * mrSges(5,1) - t8 * mrSges(6,1) + t11 * mrSges(5,2) + t7 * mrSges(6,2) + (m(3) * pkin(1) + t64 + t63 * (-pkin(1) - pkin(6)) - t62) * t23 + (-t61 * qJ(2) + t60) * t26) * g(1), -t66 * t61, (t57 * t22 + t69) * g(3) + t66 * ((-mrSges(4,1) - t57) * t25 + t70 * t22), (mrSges(5,1) * t21 + mrSges(5,2) * t24 - t32 + t64) * t48 + (-t12 * mrSges(5,2) + t58 * t11 - t52) * g(2) + (t10 * mrSges(5,2) + t58 * t9 - t53) * g(1), -g(1) * t53 - g(2) * t52 - t32 * t48];
taug = t1(:);
