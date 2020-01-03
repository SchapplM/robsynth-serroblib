% Calculate Gravitation load on the joints for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:01:48
% EndTime: 2020-01-03 12:01:49
% DurationCPUTime: 0.24s
% Computational Cost: add. (271->61), mult. (182->62), div. (0->0), fcn. (135->10), ass. (0->37)
t34 = qJ(4) + qJ(5);
t28 = sin(t34);
t30 = cos(t34);
t63 = t30 * mrSges(6,1) - mrSges(6,2) * t28;
t36 = sin(qJ(4));
t62 = mrSges(5,2) * t36 - t63;
t61 = mrSges(6,1) * t28 + mrSges(6,2) * t30;
t60 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t38 = cos(qJ(4));
t59 = -mrSges(5,1) * t38 - mrSges(4,1) + t62;
t46 = m(6) * pkin(4) + mrSges(5,1);
t58 = -mrSges(5,2) * t38 - t46 * t36;
t35 = qJ(1) + qJ(2);
t27 = pkin(9) + t35;
t23 = cos(t27);
t57 = t61 * t23;
t29 = sin(t35);
t24 = pkin(2) * t29;
t31 = cos(t35);
t25 = pkin(2) * t31;
t22 = sin(t27);
t56 = g(2) * t22;
t37 = sin(qJ(1));
t32 = t37 * pkin(1);
t50 = t24 + t32;
t26 = pkin(4) * t38 + pkin(3);
t40 = -pkin(8) - pkin(7);
t49 = t22 * t26 + t23 * t40 + t24;
t48 = t23 * pkin(3) + t22 * pkin(7) + t25;
t47 = -m(3) * pkin(1) - mrSges(2,1);
t45 = -t22 * t40 + t23 * t26 + t25;
t42 = -t31 * mrSges(3,1) + mrSges(3,2) * t29 + t60 * t22 + t59 * t23;
t41 = -t29 * mrSges(3,1) - t31 * mrSges(3,2) + (m(5) * pkin(7) - t60) * t23 + t59 * t22;
t39 = cos(qJ(1));
t33 = t39 * pkin(1);
t16 = t22 * pkin(3);
t1 = [(-mrSges(2,2) * t39 - m(4) * t50 - m(5) * (t16 + t50) - m(6) * (t32 + t49) + t47 * t37 + t41) * g(3) + (mrSges(2,2) * t37 - m(4) * (t25 + t33) - m(5) * (t33 + t48) - m(6) * (t33 + t45) + t47 * t39 + t42) * g(2), (-m(4) * t24 - m(5) * (t16 + t24) - m(6) * t49 + t41) * g(3) + (-m(4) * t25 - m(5) * t48 - m(6) * t45 + t42) * g(2), (-m(4) - m(5) - m(6)) * g(1), (t58 * t23 - t57) * g(3) + (-t38 * t46 + t62) * g(1) + (t61 - t58) * t56, -g(1) * t63 - g(3) * t57 + t56 * t61];
taug = t1(:);
