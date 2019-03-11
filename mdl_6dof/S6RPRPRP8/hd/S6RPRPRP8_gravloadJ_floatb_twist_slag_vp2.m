% Calculate Gravitation load on the joints for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:24:26
% EndTime: 2019-03-09 03:24:28
% DurationCPUTime: 0.59s
% Computational Cost: add. (289->84), mult. (410->104), div. (0->0), fcn. (365->8), ass. (0->48)
t75 = mrSges(6,1) + mrSges(7,1);
t74 = -mrSges(6,2) + mrSges(7,3);
t70 = mrSges(6,3) + mrSges(7,2);
t22 = sin(qJ(5));
t25 = cos(qJ(5));
t73 = t74 * t22 + t75 * t25;
t72 = -m(3) - m(4);
t71 = -m(6) - m(7);
t35 = pkin(5) * t25 + qJ(6) * t22;
t69 = -m(7) * (-pkin(4) - t35) + m(6) * pkin(4) + t73;
t20 = qJ(3) + pkin(9);
t15 = sin(t20);
t16 = cos(t20);
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t68 = -t23 * mrSges(4,1) - t15 * mrSges(5,1) - t26 * mrSges(4,2) - t16 * mrSges(5,2);
t61 = pkin(3) * t26;
t67 = -m(5) * t61 - mrSges(4,1) * t26 - mrSges(5,1) * t16 + mrSges(4,2) * t23 + mrSges(5,2) * t15;
t66 = t70 * t15;
t65 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t64 = mrSges(2,2) - mrSges(3,3) + t68;
t63 = m(7) * pkin(5) + t75;
t62 = m(7) * qJ(6) + t74;
t24 = sin(qJ(1));
t60 = g(1) * t24;
t27 = cos(qJ(1));
t59 = g(2) * t27;
t58 = g(3) * t16;
t57 = t23 * pkin(3);
t54 = t15 * t24;
t53 = t16 * t24;
t52 = t16 * t27;
t51 = t22 * t27;
t50 = t24 * t22;
t49 = t24 * t25;
t48 = t25 * t27;
t47 = t27 * pkin(1) + t24 * qJ(2);
t46 = -m(5) + t71;
t18 = t27 * qJ(2);
t44 = -t24 * pkin(1) + t18;
t21 = -qJ(4) - pkin(7);
t34 = t24 * t21 + t27 * t57 + t44;
t33 = -t21 * t27 + t24 * t57 + t47;
t4 = t15 * t48 - t50;
t3 = t15 * t51 + t49;
t2 = t15 * t49 + t51;
t1 = t15 * t50 - t48;
t5 = [(-m(5) * t33 + t70 * t53 + t72 * t47 + t71 * (pkin(4) * t54 - pkin(8) * t53 + t33) - t63 * t2 - t62 * t1 + (-m(4) * pkin(7) - t65) * t27 + t64 * t24) * g(2) + (-m(3) * t44 - m(4) * t18 - m(5) * t34 + t70 * t52 + t71 * (t27 * t15 * pkin(4) - pkin(8) * t52 + t34) - t63 * t4 - t62 * t3 + t64 * t27 + (-m(4) * (-pkin(1) - pkin(7)) + t65) * t24) * g(1) (-t60 + t59) * (-t46 - t72) t67 * t60 + (t71 * (pkin(4) * t53 + pkin(8) * t54 + t24 * t61) + ((-m(7) * t35 - t73) * t16 - t66) * t24) * g(1) + (t71 * (-pkin(8) * t15 - t61) + t69 * t16 + t66 - t67) * t59 + (m(5) * t57 + t71 * (t16 * pkin(8) - t57) - t70 * t16 + t69 * t15 - t68) * g(3) (g(1) * t27 + g(2) * t24) * t46 (t63 * t22 - t62 * t25) * t58 + (-t63 * t3 + t62 * t4) * g(2) + (t63 * t1 - t62 * t2) * g(1) (-g(1) * t1 + g(2) * t3 - t22 * t58) * m(7)];
taug  = t5(:);
