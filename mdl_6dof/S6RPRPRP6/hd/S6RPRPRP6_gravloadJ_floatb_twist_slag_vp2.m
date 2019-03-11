% Calculate Gravitation load on the joints for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:18:04
% EndTime: 2019-03-09 03:18:05
% DurationCPUTime: 0.67s
% Computational Cost: add. (327->90), mult. (412->89), div. (0->0), fcn. (354->8), ass. (0->45)
t67 = -mrSges(6,1) - mrSges(7,1);
t66 = m(7) * pkin(5) - t67;
t65 = mrSges(6,2) + mrSges(7,2);
t76 = -mrSges(4,1) + mrSges(5,2);
t75 = mrSges(4,2) - mrSges(5,3);
t21 = sin(qJ(5));
t23 = cos(qJ(5));
t74 = -t66 * t21 - t65 * t23;
t19 = -qJ(6) - pkin(8);
t73 = -m(6) * (-pkin(3) - pkin(8)) + mrSges(6,3) - m(7) * (-pkin(3) + t19) + mrSges(7,3);
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t63 = g(1) * t24 + g(2) * t22;
t69 = m(5) + m(7);
t16 = pkin(9) + qJ(3);
t14 = sin(t16);
t10 = t14 * qJ(4);
t15 = cos(t16);
t49 = t15 * t24;
t68 = pkin(3) * t49 + t24 * t10;
t64 = t75 * t14 + t76 * t15;
t42 = m(6) + t69;
t18 = cos(pkin(9));
t60 = -mrSges(2,1) - m(3) * pkin(1) - t18 * mrSges(3,1) + sin(pkin(9)) * mrSges(3,2) + t64;
t20 = -pkin(7) - qJ(2);
t58 = -m(6) * (pkin(4) - t20) - m(7) * (pkin(5) * t23 + pkin(4)) - mrSges(5,1) + mrSges(2,2) - mrSges(4,3) - m(3) * qJ(2) - mrSges(3,3);
t56 = pkin(5) * t21;
t53 = g(3) * t15;
t11 = t15 * pkin(3);
t51 = t15 * mrSges(7,3);
t50 = t15 * t19;
t48 = t21 * t24;
t47 = t22 * t21;
t46 = t22 * t23;
t45 = t23 * t24;
t44 = t11 + t10;
t12 = pkin(2) * t18 + pkin(1);
t8 = t24 * t12;
t40 = -t22 * t20 + t8;
t38 = -t12 - t10;
t1 = t14 * t45 - t47;
t3 = t14 * t46 + t48;
t4 = -t14 * t47 + t45;
t2 = t14 * t48 + t46;
t5 = [(-m(4) * t40 - m(6) * (pkin(8) * t49 + t68 + t8) - mrSges(6,3) * t49 - t69 * (t40 + t68) + t67 * t2 - t65 * t1 + t58 * t22 + (-m(7) * (t14 * t56 - t50) - t51 + t60) * t24) * g(2) + (t67 * t4 + t65 * t3 + ((m(4) + t69) * t20 + t58) * t24 + (m(4) * t12 - m(5) * (t38 - t11) - m(6) * t38 - m(7) * (-t12 + (-qJ(4) - t56) * t14) + t73 * t15 - t60) * t22) * g(1) (-g(1) * t22 + g(2) * t24) * (m(3) + m(4) + t42) (-m(5) * t44 - m(6) * (pkin(8) * t15 + t44) - t15 * mrSges(6,3) - m(7) * (t44 - t50) - t51 + t74 * t14 + t64) * g(3) + ((m(5) * pkin(3) + t73 - t76) * t14 + (-qJ(4) * t42 + t74 + t75) * t15) * t63 (-t14 * t63 + t53) * t42 (-t21 * t65 + t23 * t66) * t53 + (-t3 * t66 - t4 * t65) * g(2) + (-t1 * t66 + t2 * t65) * g(1) (-g(3) * t14 - t15 * t63) * m(7)];
taug  = t5(:);
