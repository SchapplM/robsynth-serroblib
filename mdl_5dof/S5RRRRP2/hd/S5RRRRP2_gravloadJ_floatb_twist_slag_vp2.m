% Calculate Gravitation load on the joints for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:05
% EndTime: 2020-01-03 12:11:06
% DurationCPUTime: 0.30s
% Computational Cost: add. (289->69), mult. (244->75), div. (0->0), fcn. (186->8), ass. (0->44)
t61 = mrSges(5,1) + mrSges(6,1);
t60 = mrSges(5,2) + mrSges(6,2);
t40 = qJ(3) + qJ(4);
t34 = cos(t40);
t70 = t61 * t34;
t44 = cos(qJ(3));
t42 = sin(qJ(3));
t32 = sin(t40);
t53 = t60 * t32;
t49 = t42 * mrSges(4,2) + t53;
t72 = -t44 * mrSges(4,1) - mrSges(3,1) + t49;
t71 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t56 = m(5) * pkin(3) + mrSges(4,1);
t69 = -t56 * t42 + m(6) * (-t42 * pkin(3) - pkin(4) * t32) - mrSges(4,2) * t44;
t68 = m(6) * pkin(4);
t46 = -pkin(8) - pkin(7);
t30 = pkin(4) * t34;
t41 = qJ(1) + qJ(2);
t33 = sin(t41);
t67 = g(2) * t33;
t37 = t44 * pkin(3);
t58 = t30 + t37;
t3 = pkin(2) + t58;
t35 = cos(t41);
t39 = -qJ(5) + t46;
t66 = t33 * t3 + t35 * t39;
t31 = t37 + pkin(2);
t65 = t33 * t31 + t35 * t46;
t64 = t34 * t35;
t63 = t35 * t32;
t59 = t35 * pkin(2) + t33 * pkin(7);
t57 = -m(3) * pkin(1) - mrSges(2,1);
t55 = t35 * t3 - t33 * t39;
t54 = t35 * t31 - t33 * t46;
t52 = t60 * t34;
t51 = -t60 * t64 - t61 * t63;
t48 = t71 * t33 + t72 * t35 - t61 * t64;
t47 = (m(4) * pkin(7) - t71) * t35 + (t72 - t70) * t33;
t45 = cos(qJ(1));
t43 = sin(qJ(1));
t38 = t45 * pkin(1);
t36 = t43 * pkin(1);
t28 = t33 * pkin(2);
t1 = [(-t45 * mrSges(2,2) - m(4) * (t28 + t36) - m(5) * (t36 + t65) - m(6) * (t36 + t66) + t57 * t43 + t47) * g(3) + (t43 * mrSges(2,2) - m(4) * (t38 + t59) - m(5) * (t38 + t54) - m(6) * (t38 + t55) + t57 * t45 + t48) * g(2), (-m(4) * t28 - m(5) * t65 - m(6) * t66 + t47) * g(3) + (-m(4) * t59 - m(5) * t54 - m(6) * t55 + t48) * g(2), (t69 * t35 + t51) * g(3) + (-m(6) * t58 - t56 * t44 + t49 - t70) * g(1) + (t61 * t32 + t52 - t69) * t67, (-t63 * t68 + t51) * g(3) + (-m(6) * t30 + t53 - t70) * g(1) + (t52 + (t61 + t68) * t32) * t67, (g(2) * t35 + g(3) * t33) * m(6)];
taug = t1(:);
