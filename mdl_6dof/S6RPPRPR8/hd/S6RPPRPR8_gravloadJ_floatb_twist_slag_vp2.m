% Calculate Gravitation load on the joints for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2018-11-23 15:43
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:43:27
% EndTime: 2018-11-23 15:43:28
% DurationCPUTime: 0.48s
% Computational Cost: add. (224->81), mult. (298->90), div. (0->0), fcn. (242->8), ass. (0->40)
t64 = mrSges(5,1) - mrSges(6,2);
t18 = pkin(9) + qJ(4);
t13 = sin(t18);
t14 = cos(t18);
t63 = -t64 * t14 + (mrSges(5,2) - mrSges(6,3)) * t13;
t23 = sin(qJ(1));
t52 = g(1) * t23;
t62 = t14 * mrSges(5,2) + t64 * t13;
t61 = pkin(4) * t14 + qJ(5) * t13;
t60 = m(6) + m(7);
t25 = cos(qJ(1));
t51 = g(2) * t25;
t59 = -t52 + t51;
t11 = t14 * qJ(5);
t19 = sin(pkin(9));
t58 = mrSges(2,2) - mrSges(3,3) - t19 * mrSges(4,1) - cos(pkin(9)) * mrSges(4,2) - (-m(6) * qJ(5) - mrSges(6,3)) * t14 - m(7) * (t13 * pkin(8) - t11) - t13 * mrSges(7,3) - t62;
t57 = m(7) * pkin(5) + mrSges(2,1) + mrSges(6,1) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t54 = pkin(3) * t19;
t50 = g(3) * t13;
t49 = t13 * pkin(4);
t22 = sin(qJ(6));
t47 = t23 * t22;
t24 = cos(qJ(6));
t46 = t23 * t24;
t45 = t25 * t22;
t44 = t25 * t24;
t43 = t25 * pkin(1) + t23 * qJ(2);
t41 = -m(4) - m(5) - t60;
t16 = t25 * qJ(2);
t40 = -t23 * pkin(1) + t16;
t38 = m(7) * (-pkin(4) - pkin(8)) - mrSges(7,3);
t33 = t22 * mrSges(7,1) + t24 * mrSges(7,2);
t21 = -pkin(7) - qJ(3);
t31 = t23 * t21 + t25 * t54 + t40;
t30 = -t25 * t21 + t23 * t54 + t43;
t4 = -t14 * t47 + t44;
t3 = -t14 * t46 - t45;
t2 = -t14 * t45 - t46;
t1 = -t14 * t44 + t47;
t5 = [(-m(5) * t30 - t4 * mrSges(7,1) - t3 * mrSges(7,2) + (-m(3) - m(4)) * t43 - t60 * (t23 * t49 + t30) + (-m(4) * qJ(3) - t57) * t25 + t58 * t23) * g(2) + (-m(3) * t40 - m(4) * t16 - m(5) * t31 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t60 * (t25 * t49 + t31) + (-m(4) * (-pkin(1) - qJ(3)) + t57) * t23 + t58 * t25) * g(1), t59 * (m(3) - t41) (g(1) * t25 + g(2) * t23) * t41 (m(6) * t61 - t38 * t14 - (-m(7) * qJ(5) - t33) * t13 - t63) * t51 + (-m(6) * (t11 - t49) - m(7) * t11 - t38 * t13 + (-mrSges(6,3) - t33) * t14 + t62) * g(3) + (-t60 * t61 - (m(7) * pkin(8) + mrSges(7,3)) * t14 - t33 * t13 + t63) * t52 (-t14 * t59 - t50) * t60, -g(1) * (t3 * mrSges(7,1) - t4 * mrSges(7,2)) - g(2) * (-t1 * mrSges(7,1) + t2 * mrSges(7,2)) - (mrSges(7,1) * t24 - mrSges(7,2) * t22) * t50];
taug  = t5(:);
