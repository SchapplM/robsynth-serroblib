% Calculate Gravitation load on the joints for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:14
% EndTime: 2019-03-08 20:18:16
% DurationCPUTime: 0.69s
% Computational Cost: add. (403->92), mult. (1019->139), div. (0->0), fcn. (1185->10), ass. (0->53)
t94 = -m(4) - m(5);
t92 = -m(6) - m(7);
t90 = mrSges(6,3) + mrSges(7,2);
t63 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t60 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t93 = -mrSges(5,1) * t48 - mrSges(5,2) * t51 + mrSges(3,2) - mrSges(4,3);
t89 = mrSges(3,1) - mrSges(4,2) + mrSges(5,3);
t47 = sin(qJ(5));
t50 = cos(qJ(5));
t88 = t60 * t47 - t63 * t50 - mrSges(5,1);
t87 = mrSges(5,2) - t90;
t86 = t93 + t92 * (-pkin(9) * t51 + qJ(3)) + t90 * t51 + t94 * qJ(3);
t85 = pkin(4) * t48;
t45 = sin(pkin(6));
t84 = t45 * t48;
t49 = sin(qJ(2));
t83 = t45 * t49;
t82 = t45 * t51;
t52 = cos(qJ(2));
t81 = t45 * t52;
t80 = t47 * t48;
t79 = t48 * t50;
t78 = t49 * t50;
t77 = pkin(2) * t81 + qJ(3) * t83;
t76 = cos(pkin(6));
t75 = t48 * t83;
t74 = t49 * t82;
t73 = pkin(8) * t81 + t77;
t72 = -t92 - t94;
t44 = sin(pkin(10));
t46 = cos(pkin(10));
t64 = t52 * t76;
t29 = t44 * t49 - t46 * t64;
t26 = t29 * pkin(2);
t70 = -pkin(8) * t29 - t26;
t31 = t44 * t64 + t46 * t49;
t27 = t31 * pkin(2);
t69 = -pkin(8) * t31 - t27;
t65 = t49 * t76;
t34 = -t48 * t81 + t76 * t51;
t33 = -t76 * t48 - t51 * t81;
t32 = -t44 * t65 + t46 * t52;
t30 = t44 * t52 + t46 * t65;
t15 = t34 * t47 - t45 * t78;
t14 = -t29 * t48 + t46 * t82;
t13 = t29 * t51 + t46 * t84;
t12 = t31 * t48 + t44 * t82;
t11 = t31 * t51 - t44 * t84;
t3 = -t14 * t47 - t30 * t50;
t1 = t12 * t47 - t32 * t50;
t2 = [(-m(2) - m(3) - t72) * g(3) (-m(4) * t77 - m(5) * t73 + t90 * t74 + t92 * (pkin(4) * t75 - pkin(9) * t74 + t73) + t60 * (t47 * t75 - t50 * t81) + (-t63 * (t47 * t52 + t48 * t78) - t89 * t52 + t93 * t49) * t45) * g(3) + (m(4) * t26 - m(5) * t70 + t92 * (t30 * t85 + t70) - t63 * (-t29 * t47 + t30 * t79) + t60 * (t29 * t50 + t30 * t80) + t89 * t29 + t86 * t30) * g(2) + (m(4) * t27 - m(5) * t69 + t92 * (t32 * t85 + t69) - t63 * (-t31 * t47 + t32 * t79) + t60 * (t31 * t50 + t32 * t80) + t89 * t31 + t86 * t32) * g(1) (-g(1) * t31 - g(2) * t29 + g(3) * t81) * t72 (t92 * (t33 * pkin(4) + pkin(9) * t34) + t87 * t34 + t88 * t33) * g(3) + (t92 * (t13 * pkin(4) - pkin(9) * t14) - t87 * t14 + t88 * t13) * g(2) + (t92 * (t11 * pkin(4) + pkin(9) * t12) + t87 * t12 + t88 * t11) * g(1) (t60 * (t34 * t50 + t47 * t83) + t63 * t15) * g(3) + (t60 * (-t14 * t50 + t30 * t47) + t63 * t3) * g(2) + (t60 * (t12 * t50 + t32 * t47) + t63 * t1) * g(1) (-g(1) * t1 - g(2) * t3 - g(3) * t15) * m(7)];
taug  = t2(:);
