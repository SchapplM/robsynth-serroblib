% Calculate Gravitation load on the joints for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:46
% EndTime: 2019-03-08 20:26:49
% DurationCPUTime: 0.99s
% Computational Cost: add. (673->97), mult. (1601->144), div. (0->0), fcn. (1984->14), ass. (0->49)
t38 = qJ(5) + qJ(6);
t36 = sin(t38);
t37 = cos(t38);
t43 = sin(qJ(5));
t46 = cos(qJ(5));
t97 = -m(7) * (pkin(5) * t46 + pkin(4)) - t37 * mrSges(7,1) + t36 * mrSges(7,2) - m(6) * pkin(4) - t46 * mrSges(6,1) + t43 * mrSges(6,2) - mrSges(5,1);
t96 = m(7) * (-pkin(10) - pkin(9)) - mrSges(7,3) - m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t107 = m(5) + m(6);
t99 = t43 * mrSges(6,1) + t46 * mrSges(6,2) - mrSges(4,2) + mrSges(5,3) + t107 * pkin(8) + m(7) * (pkin(5) * t43 + pkin(8)) + t36 * mrSges(7,1) + t37 * mrSges(7,2);
t44 = sin(qJ(4));
t47 = cos(qJ(4));
t101 = t96 * t44 + t97 * t47 - mrSges(4,1);
t45 = sin(qJ(2));
t77 = sin(pkin(12));
t78 = cos(pkin(12));
t89 = cos(qJ(2));
t54 = t45 * t78 + t89 * t77;
t39 = sin(pkin(11));
t41 = cos(pkin(11));
t42 = cos(pkin(6));
t73 = t42 * t89;
t104 = -t39 * t45 + t41 * t73;
t103 = -m(7) - t107;
t98 = -m(7) * pkin(5) - mrSges(6,1);
t26 = t45 * t77 - t89 * t78;
t79 = t54 * t42;
t18 = -t41 * t26 - t39 * t79;
t13 = t39 * t26 - t41 * t79;
t51 = t42 * t26;
t14 = -t39 * t54 - t41 * t51;
t40 = sin(pkin(6));
t86 = t40 * t44;
t8 = -t13 * t47 - t41 * t86;
t92 = (-t14 * t37 - t36 * t8) * mrSges(7,1) + (t14 * t36 - t37 * t8) * mrSges(7,2);
t10 = t18 * t47 + t39 * t86;
t17 = t39 * t51 - t41 * t54;
t91 = (-t10 * t36 - t17 * t37) * mrSges(7,1) + (-t10 * t37 + t17 * t36) * mrSges(7,2);
t25 = t54 * t40;
t20 = t25 * t47 + t42 * t44;
t24 = t26 * t40;
t90 = (-t20 * t36 + t24 * t37) * mrSges(7,1) + (-t20 * t37 - t24 * t36) * mrSges(7,2);
t85 = t40 * t47;
t83 = t42 * t45;
t76 = m(4) - t103;
t33 = t40 * t89 * pkin(2);
t68 = t104 * pkin(2);
t57 = -t39 * t73 - t41 * t45;
t55 = t57 * pkin(2);
t1 = [(-m(2) - m(3) - t76) * g(3) (-(mrSges(3,1) * t89 - mrSges(3,2) * t45) * t40 - m(4) * t33 + t103 * (-t24 * pkin(3) + t33) - t99 * t25 - t101 * t24) * g(3) + (-t104 * mrSges(3,1) - (-t39 * t89 - t41 * t83) * mrSges(3,2) - m(4) * t68 + t103 * (t14 * pkin(3) + t68) + t101 * t14 + t99 * t13) * g(2) + (-t57 * mrSges(3,1) - (t39 * t83 - t41 * t89) * mrSges(3,2) - m(4) * t55 + t103 * (t17 * pkin(3) + t55) + t101 * t17 - t99 * t18) * g(1) (-t42 * g(3) + (-g(1) * t39 + g(2) * t41) * t40) * t76 (t96 * t20 + t97 * (-t25 * t44 + t42 * t47)) * g(3) + (t96 * t8 + t97 * (t13 * t44 - t41 * t85)) * g(2) + (t97 * (-t18 * t44 + t39 * t85) + t96 * t10) * g(1) (-(-t20 * t46 - t24 * t43) * mrSges(6,2) - t90 + t98 * (-t20 * t43 + t24 * t46)) * g(3) + (-(t14 * t43 - t46 * t8) * mrSges(6,2) - t92 + t98 * (-t14 * t46 - t43 * t8)) * g(2) + (-(-t10 * t46 + t17 * t43) * mrSges(6,2) - t91 + t98 * (-t10 * t43 - t17 * t46)) * g(1), -g(1) * t91 - g(2) * t92 - g(3) * t90];
taug  = t1(:);
