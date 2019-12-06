% Calculate Gravitation load on the joints for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:01
% EndTime: 2019-12-05 16:35:05
% DurationCPUTime: 0.77s
% Computational Cost: add. (346->79), mult. (900->124), div. (0->0), fcn. (1073->12), ass. (0->50)
t41 = sin(pkin(10));
t43 = cos(pkin(10));
t84 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t90 = m(5) + m(6);
t44 = sin(qJ(5));
t47 = cos(qJ(5));
t98 = m(6) * pkin(4) + t47 * mrSges(6,1) - t44 * mrSges(6,2) + mrSges(5,1);
t101 = pkin(3) * t90 - t84 * t41 + t98 * t43;
t100 = t44 * mrSges(6,1) + t47 * mrSges(6,2) + t90 * qJ(4) + mrSges(5,3);
t96 = -m(4) - t90;
t93 = mrSges(4,1) + t101;
t88 = mrSges(3,2) - mrSges(4,3);
t92 = -t41 * t98 - t84 * t43 + t88;
t45 = sin(qJ(3));
t48 = cos(qJ(3));
t89 = -mrSges(4,1) * t48 + mrSges(4,2) * t45 - mrSges(3,1);
t91 = t100 * t45 + t101 * t48 - t89;
t83 = mrSges(4,2) - t100;
t42 = sin(pkin(5));
t46 = sin(qJ(2));
t76 = t42 * t46;
t49 = cos(qJ(2));
t75 = t42 * t49;
t71 = t48 * t49;
t70 = pkin(2) * t75 + pkin(7) * t76;
t68 = cos(pkin(5));
t67 = cos(pkin(9));
t66 = sin(pkin(9));
t65 = t45 * t75;
t64 = t41 * t75;
t61 = t42 * t67;
t60 = t42 * t66;
t59 = t42 * pkin(3) * t71 + qJ(4) * t65 + t70;
t54 = t68 * t67;
t53 = t68 * t66;
t30 = t68 * t45 + t48 * t76;
t29 = t45 * t76 - t68 * t48;
t28 = -t46 * t53 + t67 * t49;
t27 = t67 * t46 + t49 * t53;
t26 = t46 * t54 + t66 * t49;
t25 = t66 * t46 - t49 * t54;
t15 = (t41 * t46 + t43 * t71) * t42;
t13 = t28 * t48 + t45 * t60;
t12 = t28 * t45 - t48 * t60;
t11 = t26 * t48 - t45 * t61;
t10 = t26 * t45 + t48 * t61;
t9 = t30 * t43 - t64;
t2 = t13 * t43 + t27 * t41;
t1 = t11 * t43 + t25 * t41;
t3 = [(-m(2) - m(3) + t96) * g(3), (-m(4) * t70 - m(5) * t59 - t15 * mrSges(5,1) - mrSges(5,3) * t65 - m(6) * (t15 * pkin(4) + t59) - (t15 * t47 + t44 * t65) * mrSges(6,1) - (-t15 * t44 + t47 * t65) * mrSges(6,2) + (t88 * t46 + t89 * t49) * t42 + t84 * (-t43 * t76 + t48 * t64)) * g(3) + (t96 * (-t25 * pkin(2) + t26 * pkin(7)) + t92 * t26 + t91 * t25) * g(2) + (t96 * (-t27 * pkin(2) + t28 * pkin(7)) + t92 * t28 + t91 * t27) * g(1), (t93 * t29 + t83 * t30) * g(3) + (t93 * t10 + t83 * t11) * g(2) + (t93 * t12 + t83 * t13) * g(1), t90 * (-g(1) * t12 - g(2) * t10 - g(3) * t29), -g(1) * ((t12 * t47 - t2 * t44) * mrSges(6,1) + (-t12 * t44 - t2 * t47) * mrSges(6,2)) - g(2) * ((-t1 * t44 + t10 * t47) * mrSges(6,1) + (-t1 * t47 - t10 * t44) * mrSges(6,2)) - g(3) * ((t29 * t47 - t9 * t44) * mrSges(6,1) + (-t29 * t44 - t9 * t47) * mrSges(6,2))];
taug = t3(:);
