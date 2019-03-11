% Calculate Gravitation load on the joints for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:44
% EndTime: 2019-03-08 20:04:46
% DurationCPUTime: 0.76s
% Computational Cost: add. (518->83), mult. (898->122), div. (0->0), fcn. (1020->12), ass. (0->48)
t45 = cos(qJ(5));
t93 = -m(6) * pkin(4) - m(7) * (pkin(5) * t45 + pkin(4)) - mrSges(5,1);
t81 = -m(6) * pkin(9) + m(7) * (-qJ(6) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t92 = mrSges(6,1) + mrSges(7,1);
t86 = -mrSges(6,2) - mrSges(7,2);
t82 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t79 = m(7) * pkin(5);
t83 = -t79 - t92;
t37 = pkin(11) + qJ(4);
t35 = sin(t37);
t36 = cos(t37);
t40 = cos(pkin(11));
t91 = -mrSges(3,1) - m(4) * pkin(2) - t40 * mrSges(4,1) + sin(pkin(11)) * mrSges(4,2) + t93 * t36 + t81 * t35;
t85 = m(5) + m(6) + m(7);
t43 = sin(qJ(5));
t84 = t86 * t43 + t92 * t45 - t93;
t44 = sin(qJ(2));
t46 = cos(qJ(2));
t64 = cos(pkin(10));
t65 = cos(pkin(6));
t53 = t65 * t64;
t63 = sin(pkin(10));
t24 = t44 * t53 + t46 * t63;
t77 = t24 * t43;
t52 = t65 * t63;
t26 = -t44 * t52 + t46 * t64;
t76 = t26 * t43;
t73 = t36 * t43;
t72 = t36 * t45;
t39 = sin(pkin(6));
t71 = t39 * t44;
t70 = t39 * t46;
t69 = t43 * t46;
t68 = t45 * t46;
t62 = m(4) + t85;
t58 = t39 * t64;
t57 = t39 * t63;
t42 = -pkin(8) - qJ(3);
t33 = pkin(3) * t40 + pkin(2);
t25 = t44 * t64 + t46 * t52;
t23 = t44 * t63 - t46 * t53;
t20 = t35 * t65 + t36 * t71;
t19 = t35 * t71 - t36 * t65;
t12 = t26 * t36 + t35 * t57;
t11 = t26 * t35 - t36 * t57;
t10 = t24 * t36 - t35 * t58;
t9 = t24 * t35 + t36 * t58;
t1 = [(-m(2) - m(3) - t62) * g(3) (-t77 * t79 - t85 * (-t23 * t33 - t24 * t42) - t92 * (-t23 * t72 + t77) + t86 * (t23 * t73 + t24 * t45) + t82 * t24 - t91 * t23) * g(2) + (-t76 * t79 - t92 * (-t25 * t72 + t76) + t86 * (t25 * t73 + t26 * t45) - t85 * (-t25 * t33 - t26 * t42) + t82 * t26 - t91 * t25) * g(1) + (-t85 * t33 * t70 + (t91 * t46 + (-t68 * t92 - t86 * t69) * t36 + (t85 * t42 + t83 * t43 + t86 * t45 + t82) * t44) * t39) * g(3) (-g(1) * t25 - g(2) * t23 + g(3) * t70) * t62 (t84 * t19 + t81 * t20) * g(3) + (t81 * t10 + t84 * t9) * g(2) + (t84 * t11 + t81 * t12) * g(1) (t86 * (-t20 * t45 + t39 * t69) + t83 * (-t20 * t43 - t39 * t68)) * g(3) + (t86 * (-t10 * t45 - t23 * t43) + t83 * (-t10 * t43 + t23 * t45)) * g(2) + (t86 * (-t12 * t45 - t25 * t43) + t83 * (-t12 * t43 + t25 * t45)) * g(1) (-g(1) * t11 - g(2) * t9 - g(3) * t19) * m(7)];
taug  = t1(:);
