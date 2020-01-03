% Calculate Gravitation load on the joints for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:55
% EndTime: 2019-12-31 19:50:56
% DurationCPUTime: 0.33s
% Computational Cost: add. (267->60), mult. (219->63), div. (0->0), fcn. (169->8), ass. (0->32)
t26 = pkin(8) + qJ(4);
t21 = sin(t26);
t69 = mrSges(5,2) - mrSges(6,3);
t73 = t21 * t69;
t72 = mrSges(5,1) + mrSges(6,1);
t29 = cos(pkin(8));
t71 = -mrSges(3,1) - t29 * mrSges(4,1) + sin(pkin(8)) * mrSges(4,2);
t22 = cos(t26);
t68 = t22 * t72 - t73;
t67 = -mrSges(6,2) - mrSges(5,3) - mrSges(4,3) + mrSges(3,2);
t27 = qJ(1) + qJ(2);
t23 = sin(t27);
t24 = cos(t27);
t64 = g(1) * t24 + g(2) * t23;
t19 = pkin(3) * t29 + pkin(2);
t50 = qJ(5) * t21;
t40 = pkin(4) * t22 + t50;
t63 = t67 * t24 + (-m(6) * (-t19 - t40) + t68 - t71) * t23;
t55 = t22 * t24;
t62 = t67 * t23 - t72 * t55 + (t71 + t73) * t24;
t31 = sin(qJ(1));
t59 = t31 * pkin(1);
t32 = cos(qJ(1));
t25 = t32 * pkin(1);
t30 = -pkin(7) - qJ(3);
t54 = t24 * t30;
t51 = pkin(2) * t24 + qJ(3) * t23;
t47 = -pkin(2) * t23 + qJ(3) * t24;
t46 = t19 * t24 - t23 * t30;
t44 = pkin(4) * t55 + t24 * t50 + t46;
t39 = -t19 * t23 - t54;
t1 = [(-t32 * mrSges(2,1) + t31 * mrSges(2,2) - m(3) * t25 - m(4) * (t25 + t51) - m(5) * (t25 + t46) - m(6) * (t25 + t44) + t62) * g(2) + (t31 * mrSges(2,1) + t32 * mrSges(2,2) + m(3) * t59 - m(4) * (t47 - t59) - m(5) * (t39 - t59) - m(6) * (-t54 - t59) + t63) * g(1), (-m(4) * t51 - m(5) * t46 - m(6) * t44 + t62) * g(2) + (-m(4) * t47 - m(5) * t39 + m(6) * t54 + t63) * g(1), (-g(1) * t23 + g(2) * t24) * (m(4) + m(5) + m(6)), (-m(6) * t40 - t68) * g(3) + ((-m(6) * qJ(5) + t69) * t22 + (m(6) * pkin(4) + t72) * t21) * t64, (g(3) * t22 - t21 * t64) * m(6)];
taug = t1(:);
