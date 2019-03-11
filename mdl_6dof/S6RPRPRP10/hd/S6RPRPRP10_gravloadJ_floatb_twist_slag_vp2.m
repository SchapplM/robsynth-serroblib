% Calculate Gravitation load on the joints for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:31:03
% EndTime: 2019-03-09 03:31:05
% DurationCPUTime: 0.68s
% Computational Cost: add. (204->88), mult. (435->103), div. (0->0), fcn. (387->6), ass. (0->45)
t73 = mrSges(6,1) + mrSges(7,1);
t72 = mrSges(6,2) - mrSges(7,3);
t64 = -m(6) - m(7);
t21 = sin(qJ(3));
t24 = cos(qJ(3));
t71 = -t24 * mrSges(4,2) + (-mrSges(4,1) + mrSges(5,2)) * t21;
t70 = -m(3) - m(4);
t69 = -mrSges(6,3) - mrSges(7,2);
t20 = sin(qJ(5));
t23 = cos(qJ(5));
t68 = t73 * t20 + t72 * t23;
t67 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,1);
t66 = mrSges(2,2) - mrSges(3,3) - (-m(5) * qJ(4) - mrSges(5,3)) * t24 + t71;
t59 = m(7) * pkin(5) + t73;
t65 = m(7) * qJ(6) - t72;
t30 = pkin(5) * t20 - qJ(6) * t23;
t63 = -m(7) * t30 - t68;
t25 = cos(qJ(1));
t54 = g(2) * t25;
t22 = sin(qJ(1));
t55 = g(1) * t22;
t61 = -t55 + t54;
t43 = m(5) - t64;
t60 = m(5) * pkin(3) - t69 + t64 * (-pkin(3) - pkin(8));
t57 = -pkin(1) - pkin(7);
t53 = g(3) * t21;
t44 = qJ(4) * t21;
t48 = t22 * t24;
t52 = pkin(3) * t48 + t22 * t44;
t50 = t21 * t22;
t49 = t21 * t25;
t47 = t24 * t25;
t16 = t25 * qJ(2);
t46 = pkin(3) * t49 + t16;
t45 = t25 * pkin(1) + t22 * qJ(2);
t15 = t24 * qJ(4);
t41 = t25 * pkin(7) + t45;
t40 = pkin(3) * t50 + t41;
t35 = mrSges(4,1) * t24 - mrSges(4,2) * t21;
t31 = -t24 * mrSges(5,2) + t21 * mrSges(5,3);
t4 = -t20 * t48 + t23 * t25;
t3 = t20 * t25 + t23 * t48;
t2 = t20 * t47 + t22 * t23;
t1 = t20 * t22 - t23 * t47;
t5 = [(-m(3) * t45 - m(4) * t41 - m(5) * t40 + t69 * t50 + t64 * (t25 * pkin(4) + pkin(8) * t50 - t22 * t15 + t40) - t59 * t4 - t65 * t3 + t67 * t25 + t66 * t22) * g(2) + (-m(5) * t46 + t64 * (pkin(8) * t49 - t25 * t15 + (-pkin(4) + t57) * t22 + t46) + t69 * t49 + t59 * t2 + t70 * t16 + t65 * t1 + t66 * t25 + (m(3) * pkin(1) + (-m(4) - m(5)) * t57 - t67) * t22) * g(1), t61 * (t43 - t70) -t35 * t55 + (-m(5) * t52 + t64 * (pkin(8) * t48 + t52) + (t63 * t21 + t69 * t24 - t31) * t22) * g(1) + (m(5) * t44 + t31 + t35 + t60 * t24 + (m(6) * qJ(4) - m(7) * (-qJ(4) - t30) + t68) * t21) * t54 + (t60 * t21 - t43 * t15 + (-mrSges(5,3) + t63) * t24 - t71) * g(3) (-t24 * t61 - t53) * t43 (-t20 * t65 - t23 * t59) * t53 + (t1 * t59 - t2 * t65) * g(2) + (t3 * t59 - t4 * t65) * g(1) (-g(1) * t3 - g(2) * t1 + t23 * t53) * m(7)];
taug  = t5(:);
