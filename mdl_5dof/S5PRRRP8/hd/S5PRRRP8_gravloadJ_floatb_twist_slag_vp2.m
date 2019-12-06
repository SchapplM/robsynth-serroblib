% Calculate Gravitation load on the joints for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:14
% EndTime: 2019-12-05 16:58:17
% DurationCPUTime: 0.62s
% Computational Cost: add. (348->74), mult. (890->115), div. (0->0), fcn. (1053->10), ass. (0->46)
t47 = sin(qJ(4));
t50 = cos(qJ(4));
t59 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t60 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t95 = t59 * t47 - t60 * t50;
t94 = mrSges(5,3) + mrSges(6,2);
t90 = -m(5) - m(6);
t93 = -m(4) + t90;
t88 = mrSges(3,2) - mrSges(4,3);
t92 = -t60 * t47 - t59 * t50 + t88;
t48 = sin(qJ(3));
t51 = cos(qJ(3));
t89 = -mrSges(4,1) * t51 + mrSges(4,2) * t48 - mrSges(3,1);
t91 = -t89 + (-t90 * pkin(8) + t94) * t48 + (-t90 * pkin(3) - t95) * t51;
t86 = -mrSges(4,1) + t95;
t85 = mrSges(4,2) - t94;
t46 = sin(pkin(5));
t49 = sin(qJ(2));
t81 = t46 * t49;
t80 = t46 * t51;
t52 = cos(qJ(2));
t79 = t46 * t52;
t74 = t51 * t52;
t73 = pkin(2) * t79 + pkin(7) * t81;
t72 = cos(pkin(5));
t71 = cos(pkin(9));
t70 = t48 * t79;
t69 = t46 * t74;
t45 = sin(pkin(9));
t63 = t45 * t72;
t62 = t46 * t71;
t54 = t72 * t71;
t35 = t72 * t48 + t49 * t80;
t34 = -t48 * t81 + t72 * t51;
t33 = -t49 * t63 + t71 * t52;
t32 = t71 * t49 + t52 * t63;
t31 = t45 * t52 + t49 * t54;
t30 = t45 * t49 - t52 * t54;
t15 = t35 * t47 + t50 * t79;
t14 = t45 * t46 * t48 + t33 * t51;
t13 = -t33 * t48 + t45 * t80;
t12 = t31 * t51 - t48 * t62;
t11 = -t31 * t48 - t51 * t62;
t3 = t14 * t47 - t32 * t50;
t1 = t12 * t47 - t30 * t50;
t2 = [(-m(2) - m(3) + t93) * g(3), (-m(4) * t73 - t94 * t70 + t90 * (pkin(3) * t69 + pkin(8) * t70 + t73) + t59 * (t47 * t69 - t50 * t81) + (t88 * t49 + t89 * t52 - t60 * (t47 * t49 + t50 * t74)) * t46) * g(3) + (t93 * (-t30 * pkin(2) + t31 * pkin(7)) + t92 * t31 + t91 * t30) * g(2) + (t93 * (-t32 * pkin(2) + t33 * pkin(7)) + t92 * t33 + t91 * t32) * g(1), (t90 * (t34 * pkin(3) + t35 * pkin(8)) + t85 * t35 + t86 * t34) * g(3) + (t90 * (t11 * pkin(3) + t12 * pkin(8)) + t85 * t12 + t86 * t11) * g(2) + (t90 * (t13 * pkin(3) + t14 * pkin(8)) + t85 * t14 + t86 * t13) * g(1), (t59 * (t35 * t50 - t47 * t79) + t60 * t15) * g(3) + (t59 * (t12 * t50 + t30 * t47) + t60 * t1) * g(2) + (t59 * (t14 * t50 + t32 * t47) + t60 * t3) * g(1), (-g(1) * t3 - g(2) * t1 - g(3) * t15) * m(6)];
taug = t2(:);
