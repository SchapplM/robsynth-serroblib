% Calculate Gravitation load on the joints for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:53:47
% EndTime: 2020-01-03 11:53:49
% DurationCPUTime: 0.35s
% Computational Cost: add. (252->56), mult. (169->60), div. (0->0), fcn. (125->10), ass. (0->34)
t32 = qJ(4) + qJ(5);
t27 = sin(t32);
t28 = cos(t32);
t61 = t28 * mrSges(6,1) - mrSges(6,2) * t27;
t33 = sin(qJ(4));
t60 = mrSges(5,2) * t33 - t61;
t59 = mrSges(6,1) * t27 + mrSges(6,2) * t28;
t58 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t35 = cos(qJ(4));
t57 = -mrSges(5,1) * t35 - mrSges(4,1) + t60;
t43 = m(6) * pkin(4) + mrSges(5,1);
t56 = -mrSges(5,2) * t35 - t43 * t33;
t31 = qJ(1) + pkin(9);
t26 = qJ(3) + t31;
t21 = sin(t26);
t22 = cos(t26);
t23 = pkin(4) * t35 + pkin(3);
t37 = -pkin(8) - pkin(7);
t55 = t21 * t23 + t22 * t37;
t54 = t59 * t22;
t53 = g(2) * t21;
t47 = t22 * pkin(3) + t21 * pkin(7);
t24 = sin(t31);
t34 = sin(qJ(1));
t46 = t34 * pkin(1) + pkin(2) * t24;
t25 = cos(t31);
t36 = cos(qJ(1));
t45 = t36 * pkin(1) + pkin(2) * t25;
t44 = -m(3) * pkin(1) - mrSges(2,1);
t42 = -t21 * t37 + t22 * t23;
t39 = t58 * t21 + t57 * t22;
t38 = (m(5) * pkin(7) - t58) * t22 + t57 * t21;
t16 = t21 * pkin(3);
t1 = [(-mrSges(2,2) * t36 - mrSges(3,1) * t24 - mrSges(3,2) * t25 - m(4) * t46 - m(5) * (t16 + t46) - m(6) * (t46 + t55) + t44 * t34 + t38) * g(3) + (mrSges(2,2) * t34 - mrSges(3,1) * t25 + mrSges(3,2) * t24 - m(4) * t45 - m(5) * (t45 + t47) - m(6) * (t42 + t45) + t44 * t36 + t39) * g(2), (-m(3) - m(4) - m(5) - m(6)) * g(1), (-m(5) * t16 - m(6) * t55 + t38) * g(3) + (-m(5) * t47 - m(6) * t42 + t39) * g(2), (t56 * t22 - t54) * g(3) + (-t43 * t35 + t60) * g(1) + (t59 - t56) * t53, -g(1) * t61 - g(3) * t54 + t53 * t59];
taug = t1(:);
