% Calculate Gravitation load on the joints for
% S5RRPRP11
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:23
% EndTime: 2019-12-31 20:12:25
% DurationCPUTime: 0.59s
% Computational Cost: add. (179->65), mult. (395->75), div. (0->0), fcn. (361->6), ass. (0->36)
t63 = -m(5) - m(6);
t54 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t53 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t70 = -mrSges(3,1) + mrSges(4,2);
t69 = mrSges(3,2) - mrSges(4,3);
t62 = -mrSges(5,3) - mrSges(6,2);
t19 = sin(qJ(4));
t22 = cos(qJ(4));
t68 = -t54 * t19 + t53 * t22;
t67 = -t62 + t63 * (-pkin(2) - pkin(7));
t21 = sin(qJ(1));
t24 = cos(qJ(1));
t60 = g(1) * t24 + g(2) * t21;
t20 = sin(qJ(2));
t23 = cos(qJ(2));
t61 = t69 * t20 + t70 * t23;
t59 = m(4) - t63;
t58 = -mrSges(2,1) + t61;
t57 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t49 = g(3) * t23;
t15 = t23 * pkin(2);
t48 = t21 * t19;
t47 = t21 * t22;
t46 = t24 * t19;
t45 = t24 * t22;
t44 = t24 * t23;
t11 = t20 * qJ(3);
t43 = t15 + t11;
t41 = t24 * pkin(1) + t21 * pkin(6);
t36 = pkin(2) * t44 + t24 * t11 + t41;
t16 = t24 * pkin(6);
t4 = -t20 * t48 + t45;
t3 = t20 * t47 + t46;
t2 = t20 * t46 + t47;
t1 = -t20 * t45 + t48;
t5 = [(-m(3) * t41 - m(4) * t36 + t62 * t44 + t63 * (t21 * pkin(3) + pkin(7) * t44 + t36) - t54 * t2 - t53 * t1 + t58 * t24 + t57 * t21) * g(2) + (t63 * (t24 * pkin(3) + t16) - t54 * t4 - t53 * t3 + (-m(3) - m(4)) * t16 + t57 * t24 + (m(3) * pkin(1) + m(4) * t15 - t59 * (-pkin(1) - t11) + t67 * t23 - t58) * t21) * g(1), (-m(4) * t43 + t63 * (t23 * pkin(7) + t43) + t62 * t23 + t68 * t20 + t61) * g(3) + ((m(4) * pkin(2) + t67 - t70) * t20 + (-qJ(3) * t59 + t68 + t69) * t23) * t60, (-t60 * t20 + t49) * t59, (t53 * t19 + t54 * t22) * t49 + (-t54 * t3 + t53 * t4) * g(2) + (t54 * t1 - t53 * t2) * g(1), (-g(1) * t1 + g(2) * t3 - t22 * t49) * m(6)];
taug = t5(:);
