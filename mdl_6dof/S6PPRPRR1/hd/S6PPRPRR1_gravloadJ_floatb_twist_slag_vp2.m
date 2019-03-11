% Calculate Gravitation load on the joints for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:42:03
% EndTime: 2019-03-08 18:42:06
% DurationCPUTime: 0.84s
% Computational Cost: add. (770->97), mult. (2124->168), div. (0->0), fcn. (2705->16), ass. (0->64)
t96 = m(6) + m(7);
t50 = sin(qJ(6));
t53 = cos(qJ(6));
t93 = -m(7) * pkin(5) - t53 * mrSges(7,1) + t50 * mrSges(7,2) - mrSges(6,1);
t92 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t41 = sin(pkin(12));
t46 = cos(pkin(12));
t52 = sin(qJ(3));
t48 = cos(pkin(7));
t55 = cos(qJ(3));
t76 = t48 * t55;
t95 = -t41 * t52 + t46 * t76;
t91 = m(5) + t96;
t51 = sin(qJ(5));
t54 = cos(qJ(5));
t90 = t92 * t51 + t93 * t54 - mrSges(5,1);
t43 = sin(pkin(7));
t40 = sin(pkin(13));
t45 = cos(pkin(13));
t68 = t40 * t55 + t52 * t45;
t26 = t68 * t43;
t28 = t68 * t48;
t34 = t40 * t52 - t55 * t45;
t44 = sin(pkin(6));
t49 = cos(pkin(6));
t18 = t26 * t49 + (t28 * t46 - t34 * t41) * t44;
t89 = t50 * mrSges(7,1) + t53 * mrSges(7,2) + t96 * pkin(9) - mrSges(5,2) + mrSges(6,3);
t42 = sin(pkin(11));
t47 = cos(pkin(11));
t77 = t47 * t49;
t31 = t41 * t77 + t42 * t46;
t86 = t31 * t52;
t82 = t42 * t49;
t33 = -t41 * t82 + t46 * t47;
t85 = t33 * t52;
t83 = t42 * t44;
t81 = t43 * t44;
t80 = t43 * t55;
t79 = t44 * t47;
t78 = t44 * t48;
t75 = pkin(3) * t76;
t74 = t44 * t80;
t71 = m(3) + m(4) + t91;
t32 = -t41 * t47 - t46 * t82;
t67 = t32 * t75 + (t42 * t74 - t85) * pkin(3);
t30 = -t41 * t42 + t46 * t77;
t66 = -t30 * t48 + t43 * t79;
t65 = t32 * t48 + t42 * t81;
t63 = (t95 * t44 + t49 * t80) * pkin(3);
t7 = t26 * t79 - t28 * t30 + t31 * t34;
t12 = t26 * t83 + t28 * t32 - t33 * t34;
t58 = t30 * t75 + (-t47 * t74 - t86) * pkin(3);
t29 = -t46 * t81 + t48 * t49;
t27 = t34 * t48;
t25 = t34 * t43;
t20 = -t32 * t43 + t42 * t78;
t19 = -t30 * t43 - t47 * t78;
t17 = -t25 * t49 + (-t27 * t46 - t41 * t68) * t44;
t14 = t18 * t54 + t29 * t51;
t11 = -t25 * t83 - t27 * t32 - t33 * t68;
t8 = t25 * t79 - t27 * t30 - t31 * t68;
t4 = t12 * t54 + t20 * t51;
t2 = t19 * t51 - t54 * t7;
t1 = [(-m(2) - t71) * g(3) (-g(3) * t49 + (-g(1) * t42 + g(2) * t47) * t44) * t71 (-(mrSges(4,1) * t55 - mrSges(4,2) * t52) * t49 * t43 - (t95 * mrSges(4,1) + (-t46 * t48 * t52 - t41 * t55) * mrSges(4,2)) * t44 - m(5) * t63 - t96 * (t17 * pkin(4) + t63) + t90 * t17 - t89 * t18) * g(3) + (-(-t55 * t66 - t86) * mrSges(4,1) - (-t31 * t55 + t52 * t66) * mrSges(4,2) - m(5) * t58 - t96 * (t8 * pkin(4) + t58) + t90 * t8 + t89 * t7) * g(2) + (-(t55 * t65 - t85) * mrSges(4,1) - (-t33 * t55 - t52 * t65) * mrSges(4,2) - m(5) * t67 - t96 * (t11 * pkin(4) + t67) + t90 * t11 - t89 * t12) * g(1), t91 * (-g(1) * t20 - g(2) * t19 - g(3) * t29) (t92 * t14 + t93 * (-t18 * t51 + t29 * t54)) * g(3) + (t92 * t2 + t93 * (t19 * t54 + t51 * t7)) * g(2) + (t92 * t4 + t93 * (-t12 * t51 + t20 * t54)) * g(1), -g(1) * ((-t11 * t53 - t4 * t50) * mrSges(7,1) + (t11 * t50 - t4 * t53) * mrSges(7,2)) - g(2) * ((-t2 * t50 - t53 * t8) * mrSges(7,1) + (-t2 * t53 + t50 * t8) * mrSges(7,2)) - g(3) * ((-t14 * t50 - t17 * t53) * mrSges(7,1) + (-t14 * t53 + t17 * t50) * mrSges(7,2))];
taug  = t1(:);
