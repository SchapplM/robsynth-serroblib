% Calculate Gravitation load on the joints for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:00
% EndTime: 2019-12-31 21:11:01
% DurationCPUTime: 0.43s
% Computational Cost: add. (302->80), mult. (335->95), div. (0->0), fcn. (299->8), ass. (0->44)
t36 = sin(qJ(3));
t39 = cos(qJ(3));
t85 = (-mrSges(4,1) - mrSges(5,1)) * t39 + (mrSges(4,2) - mrSges(5,3)) * t36;
t84 = -mrSges(3,1) + t85;
t83 = -mrSges(5,2) + mrSges(6,3) - mrSges(4,3) + mrSges(3,2);
t35 = sin(qJ(5));
t38 = cos(qJ(5));
t48 = t39 * t35 - t36 * t38;
t34 = qJ(1) + qJ(2);
t29 = sin(t34);
t30 = cos(t34);
t80 = g(1) * t30 + g(2) * t29;
t32 = t39 * pkin(3);
t5 = t48 * t29;
t31 = t36 * qJ(4);
t58 = -pkin(2) - t31;
t47 = t36 * t35 + t39 * t38;
t6 = t47 * t29;
t77 = -pkin(3) - pkin(4);
t79 = t6 * mrSges(6,1) - t5 * mrSges(6,2) + t83 * t30 + (-m(5) * (t58 - t32) - m(6) * (t77 * t39 + t58) - t84) * t29;
t7 = t48 * t30;
t8 = t47 * t30;
t78 = -t8 * mrSges(6,1) + t7 * mrSges(6,2) + t83 * t29 + t84 * t30;
t37 = sin(qJ(1));
t74 = t37 * pkin(1);
t40 = cos(qJ(1));
t33 = t40 * pkin(1);
t73 = t30 * t39;
t66 = t30 * pkin(2) + t29 * pkin(7);
t65 = t32 + t31;
t64 = qJ(4) * t39;
t62 = t77 * t36;
t27 = t30 * pkin(7);
t61 = -t29 * pkin(2) + t27;
t60 = -t30 * pkin(8) + t27;
t57 = pkin(3) * t73 + t30 * t31 + t66;
t55 = -t5 * mrSges(6,1) - t6 * mrSges(6,2);
t54 = -t7 * mrSges(6,1) - t8 * mrSges(6,2);
t52 = -mrSges(6,1) * t47 + mrSges(6,2) * t48;
t44 = pkin(4) * t73 - t29 * pkin(8) + t57;
t43 = t39 * mrSges(5,3) + (-m(5) * pkin(3) - mrSges(5,1)) * t36;
t14 = t30 * t64;
t12 = t29 * t64;
t1 = [(-t40 * mrSges(2,1) + t37 * mrSges(2,2) - m(3) * t33 - m(4) * (t33 + t66) - m(5) * (t33 + t57) - m(6) * (t33 + t44) + t78) * g(2) + (t37 * mrSges(2,1) + t40 * mrSges(2,2) + m(3) * t74 - m(4) * (t61 - t74) - m(5) * (t27 - t74) - m(6) * (t60 - t74) + t79) * g(1), (-m(4) * t66 - m(5) * t57 - m(6) * t44 + t78) * g(2) + (-m(4) * t61 - m(5) * t27 - m(6) * t60 + t79) * g(1), t80 * (mrSges(4,1) * t36 + mrSges(4,2) * t39) + (-m(5) * t12 - t43 * t29 - m(6) * (t29 * t62 + t12) + t55) * g(2) + (-m(5) * t14 - t43 * t30 - m(6) * (t30 * t62 + t14) + t54) * g(1) + (-m(5) * t65 - m(6) * (t39 * pkin(4) + t65) + t52 + t85) * g(3), (g(3) * t39 - t36 * t80) * (m(5) + m(6)), -g(1) * t54 - g(2) * t55 - g(3) * t52];
taug = t1(:);
