% Calculate Gravitation load on the joints for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:28
% EndTime: 2019-12-31 20:15:29
% DurationCPUTime: 0.28s
% Computational Cost: add. (230->66), mult. (199->74), div. (0->0), fcn. (149->8), ass. (0->40)
t33 = cos(qJ(4));
t58 = mrSges(5,2) * t33;
t68 = -t58 + mrSges(3,2) - mrSges(4,3);
t67 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t30 = qJ(1) + qJ(2);
t25 = sin(t30);
t35 = -pkin(8) - pkin(7);
t27 = cos(t30);
t31 = sin(qJ(4));
t52 = t27 * t31;
t66 = pkin(4) * t52 + t25 * t35;
t45 = m(6) * pkin(4) + mrSges(5,1);
t65 = -mrSges(5,2) * t31 + t45 * t33;
t29 = qJ(4) + qJ(5);
t26 = cos(t29);
t53 = t26 * t27;
t24 = sin(t29);
t56 = t24 * t27;
t64 = -mrSges(5,1) * t52 - mrSges(6,1) * t56 - mrSges(6,2) * t53 + t68 * t27 + (-m(5) * (-pkin(2) - pkin(7)) - t67) * t25;
t54 = t25 * t31;
t55 = t25 * t26;
t57 = t24 * t25;
t63 = -mrSges(5,1) * t54 - mrSges(6,1) * t57 - mrSges(6,2) * t55 + t68 * t25 + t67 * t27;
t32 = sin(qJ(1));
t62 = pkin(1) * t32;
t61 = g(1) * t25;
t60 = g(2) * t27;
t34 = cos(qJ(1));
t28 = t34 * pkin(1);
t50 = t27 * pkin(2) + t25 * qJ(3);
t46 = t28 + t50;
t15 = t27 * qJ(3);
t44 = -pkin(2) * t25 + t15;
t39 = -mrSges(6,1) * t24 - mrSges(6,2) * t26;
t38 = pkin(4) * t54 - t27 * t35 + t50;
t37 = t44 - t62;
t22 = t27 * pkin(7);
t4 = mrSges(6,2) * t56;
t3 = mrSges(6,1) * t55;
t1 = [(-mrSges(2,1) * t34 + mrSges(2,2) * t32 - m(3) * t28 - m(4) * t46 - m(5) * (t22 + t46) - m(6) * (t28 + t38) + t63) * g(2) + (mrSges(2,1) * t32 + mrSges(2,2) * t34 + m(3) * t62 - m(4) * t37 - m(5) * (t15 - t62) - m(6) * (t37 + t66) + t64) * g(1), (-m(4) * t50 - m(5) * (t22 + t50) - m(6) * t38 + t63) * g(2) + (-m(4) * t44 - m(5) * t15 - m(6) * (t44 + t66) + t64) * g(1), (t60 - t61) * (m(4) + m(5) + m(6)), -g(1) * t3 - g(2) * t4 + (t45 * t31 - t39 + t58) * g(3) + (mrSges(6,1) * t26 + t65) * t60 + (mrSges(6,2) * t24 - t65) * t61, -g(1) * (-mrSges(6,2) * t57 + t3) - g(2) * (-mrSges(6,1) * t53 + t4) - g(3) * t39];
taug = t1(:);
