% Calculate Gravitation load on the joints for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2018-11-23 15:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:56:02
% EndTime: 2018-11-23 15:56:03
% DurationCPUTime: 0.52s
% Computational Cost: add. (172->90), mult. (349->101), div. (0->0), fcn. (288->6), ass. (0->44)
t65 = mrSges(7,3) - mrSges(6,2);
t68 = -m(7) * pkin(8) - t65;
t67 = m(7) * pkin(5);
t19 = sin(qJ(3));
t22 = cos(qJ(3));
t64 = t22 * mrSges(4,2) + (mrSges(4,1) + mrSges(5,1)) * t19;
t63 = m(6) + m(7);
t18 = sin(qJ(6));
t21 = cos(qJ(6));
t30 = -t21 * mrSges(7,1) + t18 * mrSges(7,2);
t62 = -mrSges(6,1) + t30 - t67;
t23 = cos(qJ(1));
t53 = g(2) * t23;
t20 = sin(qJ(1));
t54 = g(1) * t20;
t61 = -t54 + t53;
t42 = m(5) + t63;
t56 = -pkin(3) - pkin(4);
t60 = m(5) * pkin(3) - m(7) * (-pkin(8) + t56) - m(6) * t56 + t65;
t59 = -mrSges(2,1) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t14 = t22 * qJ(4);
t39 = -m(6) * qJ(4) - mrSges(6,1);
t58 = m(7) * t14 + mrSges(2,2) - mrSges(3,3) - t64 + (m(5) * qJ(4) + mrSges(5,3) - t39 + t67) * t22 + t68 * t19;
t43 = qJ(4) * t19;
t48 = t20 * t22;
t55 = pkin(3) * t48 + t20 * t43;
t52 = g(3) * t19;
t51 = t20 * pkin(1);
t49 = t20 * t19;
t47 = t23 * t19;
t46 = t23 * t22;
t15 = t23 * qJ(2);
t45 = pkin(3) * t47 + t15;
t44 = t23 * pkin(1) + t20 * qJ(2);
t40 = t23 * pkin(7) + t44;
t38 = pkin(3) * t49 + t40;
t37 = pkin(4) * t47 + t20 * qJ(5) + t45;
t33 = mrSges(4,1) * t22 - mrSges(4,2) * t19;
t31 = t22 * mrSges(5,1) + t19 * mrSges(5,3);
t4 = t20 * t18 - t21 * t46;
t3 = t18 * t46 + t20 * t21;
t2 = t23 * t18 + t21 * t48;
t1 = t18 * t48 - t23 * t21;
t5 = [(-m(3) * t44 - m(4) * t40 - m(5) * t38 + t2 * mrSges(7,1) - t1 * mrSges(7,2) - t63 * (pkin(4) * t49 - t23 * qJ(5) + t38) + t59 * t23 + t58 * t20) * g(2) + (-m(3) * (t15 - t51) - m(4) * t15 - m(5) * t45 - m(6) * t37 - m(7) * (t37 - t51) - t4 * mrSges(7,1) - t3 * mrSges(7,2) + (m(7) * pkin(7) + (-m(4) - m(5) - m(6)) * (-pkin(1) - pkin(7)) - t59) * t20 + t58 * t23) * g(1), t61 * (m(3) + m(4) + t42) -t33 * t54 + (-m(5) * t55 - t63 * (pkin(4) * t48 + t55) + (t62 * t19 + t68 * t22 - t31) * t20) * g(1) + (m(5) * t43 + t31 + t33 + t60 * t22 + (-t39 - m(7) * (-pkin(5) - qJ(4)) - t30) * t19) * t53 + (t60 * t19 - t42 * t14 + (-mrSges(5,3) + t62) * t22 + t64) * g(3) (-t61 * t22 - t52) * t42, t63 * (g(1) * t23 + g(2) * t20) -g(1) * (t1 * mrSges(7,1) + t2 * mrSges(7,2)) - g(2) * (-t3 * mrSges(7,1) + t4 * mrSges(7,2)) - (-mrSges(7,1) * t18 - mrSges(7,2) * t21) * t52];
taug  = t5(:);
