% Calculate Gravitation load on the joints for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:34
% EndTime: 2019-12-05 17:09:36
% DurationCPUTime: 0.44s
% Computational Cost: add. (259->72), mult. (300->89), div. (0->0), fcn. (259->10), ass. (0->46)
t29 = qJ(4) + qJ(5);
t24 = sin(t29);
t33 = sin(qJ(4));
t97 = -t33 * mrSges(5,2) - t24 * mrSges(6,2);
t26 = cos(t29);
t35 = cos(qJ(4));
t96 = t35 * mrSges(5,1) + t26 * mrSges(6,1);
t30 = qJ(2) + qJ(3);
t25 = sin(t30);
t27 = cos(t30);
t94 = -mrSges(5,3) - mrSges(6,3);
t95 = t27 * (-m(5) * pkin(7) + t94) + t97 * t25;
t90 = t96 * t25;
t89 = m(6) * pkin(4) + mrSges(5,1);
t23 = t35 * pkin(4) + pkin(3);
t37 = -pkin(8) - pkin(7);
t46 = -t23 * t25 - t27 * t37;
t78 = pkin(3) * t25;
t34 = sin(qJ(2));
t79 = pkin(2) * t34;
t88 = -m(6) * (t46 - t79) - m(5) * (-t78 - t79) + t90;
t31 = sin(pkin(9));
t32 = cos(pkin(9));
t87 = g(1) * t32 + g(2) * t31;
t86 = (-mrSges(4,1) - t96 - t97) * t27 + (mrSges(4,2) + t94) * t25;
t85 = t95 * t31;
t84 = t95 * t32;
t83 = m(5) * t78 - m(6) * t46 + t90;
t61 = t32 * t26;
t62 = t32 * t24;
t65 = t31 * t26;
t66 = t31 * t24;
t81 = (-t27 * t66 - t61) * mrSges(6,1) + (-t27 * t65 + t62) * mrSges(6,2);
t80 = (-t27 * t62 + t65) * mrSges(6,1) + (-t27 * t61 - t66) * mrSges(6,2);
t75 = g(3) * t25;
t36 = cos(qJ(2));
t28 = t36 * pkin(2);
t64 = t31 * t33;
t63 = t31 * t35;
t60 = t32 * t33;
t59 = t32 * t35;
t56 = t27 * pkin(3) + t25 * pkin(7);
t51 = t27 * t23 - t25 * t37;
t48 = mrSges(4,1) * t25 + mrSges(4,2) * t27;
t47 = -mrSges(6,1) * t24 - mrSges(6,2) * t26;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), (t31 * t88 + t85) * g(2) + (t32 * t88 + t84) * g(1) + (-t36 * mrSges(3,1) + t34 * mrSges(3,2) - m(4) * t28 - m(5) * (t28 + t56) - m(6) * (t28 + t51) + t86) * g(3) + (m(4) * t79 + mrSges(3,1) * t34 + mrSges(3,2) * t36 + t48) * t87, t87 * t48 + (t31 * t83 + t85) * g(2) + (t32 * t83 + t84) * g(1) + (-m(5) * t56 - m(6) * t51 + t86) * g(3), (mrSges(5,2) * t35 + t33 * t89 - t47) * t75 + (-(-t27 * t63 + t60) * mrSges(5,2) - t81 - t89 * (-t27 * t64 - t59)) * g(2) + (-(-t27 * t59 - t64) * mrSges(5,2) - t80 - t89 * (-t27 * t60 + t63)) * g(1), -g(1) * t80 - g(2) * t81 - t47 * t75];
taug = t1(:);
