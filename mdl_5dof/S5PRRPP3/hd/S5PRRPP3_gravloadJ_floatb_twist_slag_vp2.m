% Calculate Gravitation load on the joints for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:00
% EndTime: 2019-12-05 16:12:03
% DurationCPUTime: 0.64s
% Computational Cost: add. (175->66), mult. (451->90), div. (0->0), fcn. (456->8), ass. (0->34)
t72 = mrSges(5,1) + mrSges(6,1);
t71 = -mrSges(5,2) + mrSges(6,3);
t65 = m(5) + m(6);
t27 = sin(qJ(3));
t29 = cos(qJ(3));
t70 = -t29 * mrSges(4,1) + t27 * mrSges(4,2) - mrSges(3,1);
t69 = mrSges(3,2) - mrSges(4,3);
t68 = mrSges(5,3) + mrSges(6,2);
t23 = sin(pkin(8));
t25 = cos(pkin(8));
t67 = t71 * t23 + t72 * t25;
t63 = -m(4) - t65;
t39 = -pkin(4) * t25 - qJ(5) * t23;
t62 = -m(6) * t39 + mrSges(4,1) + t67;
t61 = mrSges(4,2) - t68;
t58 = -m(6) * pkin(4) - t72;
t57 = m(6) * qJ(5) + t71;
t28 = sin(qJ(2));
t53 = g(3) * t28;
t52 = t28 * t25;
t51 = t28 * t29;
t30 = cos(qJ(2));
t50 = t30 * t27;
t49 = t30 * t29;
t48 = t30 * pkin(2) + t28 * pkin(6);
t47 = qJ(4) * t27;
t36 = t23 * t51 + t30 * t25;
t26 = cos(pkin(7));
t24 = sin(pkin(7));
t10 = t24 * t27 + t26 * t49;
t9 = -t24 * t29 + t26 * t50;
t8 = t24 * t49 - t26 * t27;
t7 = t24 * t50 + t26 * t29;
t1 = [(-m(2) - m(3) + t63) * g(3), (-m(4) * t48 - t68 * t50 - t65 * (pkin(3) * t49 + t30 * t47 + t48) + t70 * t30 + t69 * t28 + t58 * (t28 * t23 + t25 * t49) - t57 * (t23 * t49 - t52)) * g(3) + (t26 * g(1) + t24 * g(2)) * (t36 * t57 - t25 * t51 * t58 + (m(4) * pkin(2) - t65 * (-pkin(3) * t29 - pkin(2) - t47) + t68 * t27 - t70) * t28 + (pkin(6) * t63 + t58 * t23 + t69) * t30), -(-mrSges(4,1) * t27 - mrSges(4,2) * t29) * t53 + ((-t68 * t29 + (m(5) * pkin(3) - m(6) * (-pkin(3) + t39) + t67) * t27) * t28 - t65 * qJ(4) * t51) * g(3) + (-t65 * (-t7 * pkin(3) + t8 * qJ(4)) + t61 * t8 + t62 * t7) * g(2) + (-t65 * (-t9 * pkin(3) + t10 * qJ(4)) + t62 * t9 + t61 * t10) * g(1), t65 * (-g(1) * t9 - g(2) * t7 - t27 * t53), (-g(1) * (t10 * t23 - t26 * t52) - g(2) * (t8 * t23 - t24 * t52) - g(3) * t36) * m(6)];
taug = t1(:);
