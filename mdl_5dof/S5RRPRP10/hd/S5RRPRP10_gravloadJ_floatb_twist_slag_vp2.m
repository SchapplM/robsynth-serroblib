% Calculate Gravitation load on the joints for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:25
% EndTime: 2019-12-31 20:09:27
% DurationCPUTime: 0.65s
% Computational Cost: add. (170->75), mult. (366->80), div. (0->0), fcn. (322->6), ass. (0->37)
t61 = -mrSges(5,1) - mrSges(6,1);
t60 = m(6) * pkin(4) - t61;
t59 = mrSges(5,2) + mrSges(6,2);
t68 = mrSges(3,1) - mrSges(4,2);
t67 = -mrSges(3,2) + mrSges(4,3);
t16 = sin(qJ(4));
t19 = cos(qJ(4));
t66 = -t60 * t16 - t59 * t19;
t15 = -qJ(5) - pkin(7);
t65 = -m(5) * (-pkin(2) - pkin(7)) + mrSges(5,3) - m(6) * (-pkin(2) + t15) + mrSges(6,3);
t18 = sin(qJ(1));
t21 = cos(qJ(1));
t57 = g(1) * t21 + g(2) * t18;
t17 = sin(qJ(2));
t20 = cos(qJ(2));
t58 = t67 * t17 + t20 * t68;
t56 = m(4) + m(5) + m(6);
t53 = -t20 * mrSges(6,3) - t58;
t51 = -m(5) * pkin(3) - m(6) * (t19 * pkin(4) + pkin(3)) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t49 = pkin(4) * t16;
t46 = g(3) * t20;
t12 = t20 * pkin(2);
t45 = t18 * t16;
t44 = t18 * t19;
t42 = t20 * t15;
t41 = t21 * t16;
t40 = t21 * t19;
t39 = t21 * t20;
t10 = t17 * qJ(3);
t38 = t12 + t10;
t37 = t21 * pkin(1) + t18 * pkin(6);
t32 = -pkin(1) - t10;
t1 = t17 * t40 - t45;
t3 = t17 * t44 + t41;
t4 = -t17 * t45 + t40;
t2 = t17 * t41 + t44;
t5 = [(-m(3) * t37 + (-m(5) * pkin(7) - mrSges(5,3)) * t39 - t56 * (pkin(2) * t39 + t21 * t10 + t37) + t61 * t2 - t59 * t1 + t51 * t18 + (-mrSges(2,1) - m(6) * (t17 * t49 - t42) + t53) * t21) * g(2) + (t61 * t4 + t59 * t3 + (mrSges(2,1) + m(3) * pkin(1) - m(4) * (t32 - t12) - m(5) * t32 - m(6) * (-pkin(1) + (-qJ(3) - t49) * t17) + t65 * t20 + t58) * t18 + ((-m(3) - t56) * pkin(6) + t51) * t21) * g(1), (-m(4) * t38 - m(5) * (t20 * pkin(7) + t38) - t20 * mrSges(5,3) - m(6) * (t38 - t42) + t66 * t17 + t53) * g(3) + ((m(4) * pkin(2) + t65 + t68) * t17 + (-qJ(3) * t56 + t66 - t67) * t20) * t57, (-t57 * t17 + t46) * t56, (-t59 * t16 + t60 * t19) * t46 + (-t3 * t60 - t59 * t4) * g(2) + (-t1 * t60 + t59 * t2) * g(1), (-g(3) * t17 - t57 * t20) * m(6)];
taug = t5(:);
