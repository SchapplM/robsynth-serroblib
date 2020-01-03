% Calculate Gravitation load on the joints for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:26
% EndTime: 2019-12-31 19:46:28
% DurationCPUTime: 0.65s
% Computational Cost: add. (182->78), mult. (336->85), div. (0->0), fcn. (291->8), ass. (0->40)
t72 = mrSges(3,1) - mrSges(4,2);
t71 = -mrSges(3,2) + mrSges(4,3);
t64 = m(5) + m(6);
t70 = m(4) + t64;
t17 = pkin(8) + qJ(5);
t10 = sin(t17);
t11 = cos(t17);
t18 = sin(pkin(8));
t19 = cos(pkin(8));
t30 = t18 * mrSges(5,1) + t19 * mrSges(5,2);
t53 = pkin(4) * t18;
t69 = -m(6) * t53 - t10 * mrSges(6,1) - t11 * mrSges(6,2) - t30;
t20 = -pkin(7) - qJ(4);
t68 = -m(5) * (-pkin(2) - qJ(4)) + mrSges(5,3) - m(6) * (-pkin(2) + t20) + mrSges(6,3);
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t61 = g(1) * t24 + g(2) * t22;
t21 = sin(qJ(2));
t12 = t21 * qJ(3);
t23 = cos(qJ(2));
t41 = t23 * pkin(2) + t12;
t63 = m(4) * t41;
t62 = t71 * t21 + t72 * t23;
t58 = -t23 * mrSges(6,3) - t62;
t57 = m(3) + t70;
t55 = -m(5) * pkin(3) - m(6) * (t19 * pkin(4) + pkin(3)) - t19 * mrSges(5,1) + t18 * mrSges(5,2) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t50 = g(3) * t23;
t48 = t22 * t10;
t47 = t22 * t11;
t45 = t23 * t20;
t44 = t24 * t10;
t43 = t24 * t11;
t42 = t24 * t23;
t40 = t24 * pkin(1) + t22 * pkin(6);
t38 = t23 * qJ(4);
t4 = -t21 * t48 + t43;
t3 = t21 * t47 + t44;
t2 = t21 * t44 + t47;
t1 = t21 * t43 - t48;
t5 = [(-m(3) * t40 - t2 * mrSges(6,1) - t1 * mrSges(6,2) - mrSges(5,3) * t42 - t70 * (pkin(2) * t42 + t40) + t55 * t22) * g(2) + (-t4 * mrSges(6,1) + t3 * mrSges(6,2) + (mrSges(2,1) + t63 + t68 * t23 + (m(5) * qJ(3) + t30 - m(6) * (-qJ(3) - t53)) * t21 + t57 * pkin(1) + t62) * t22) * g(1) + ((-t30 * t21 - t70 * t12 - mrSges(2,1) - m(5) * t38 - m(6) * (t21 * t53 - t45) + t58) * g(2) + (-t57 * pkin(6) + t55) * g(1)) * t24, (-t63 - m(5) * (t38 + t41) - t23 * mrSges(5,3) - m(6) * (t41 - t45) + t69 * t21 + t58) * g(3) + ((m(4) * pkin(2) + t68 + t72) * t21 + (-qJ(3) * t70 + t69 - t71) * t23) * t61, (-t61 * t21 + t50) * t70, (-t21 * g(3) - t61 * t23) * t64, -g(1) * (t1 * mrSges(6,1) - t2 * mrSges(6,2)) - g(2) * (t3 * mrSges(6,1) + t4 * mrSges(6,2)) - (-mrSges(6,1) * t11 + mrSges(6,2) * t10) * t50];
taug = t5(:);
