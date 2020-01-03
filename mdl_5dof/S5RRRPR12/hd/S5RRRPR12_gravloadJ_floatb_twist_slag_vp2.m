% Calculate Gravitation load on the joints for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:37:05
% EndTime: 2019-12-31 21:37:08
% DurationCPUTime: 1.09s
% Computational Cost: add. (467->85), mult. (1030->120), div. (0->0), fcn. (1205->12), ass. (0->48)
t37 = sin(qJ(3));
t40 = cos(qJ(3));
t45 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) + m(6) * (-pkin(9) - qJ(4));
t32 = pkin(10) + qJ(5);
t29 = sin(t32);
t30 = cos(t32);
t33 = sin(pkin(10));
t35 = cos(pkin(10));
t43 = -m(6) * (pkin(4) * t35 + pkin(3)) - m(5) * pkin(3) - t35 * mrSges(5,1) + t33 * mrSges(5,2) - mrSges(4,1);
t96 = -t30 * mrSges(6,1) + t29 * mrSges(6,2) + t43;
t102 = t45 * t37 + t96 * t40 - mrSges(3,1);
t84 = m(5) + m(6);
t97 = m(4) + t84;
t101 = pkin(8) * t97;
t50 = t29 * mrSges(6,1) + t30 * mrSges(6,2);
t93 = -t35 * mrSges(5,2) - (m(6) * pkin(4) + mrSges(5,1)) * t33 + mrSges(3,2) - mrSges(4,3);
t100 = -t50 + t93;
t95 = pkin(2) * t97 - t102;
t85 = t100 - t101;
t80 = t93 - t101;
t78 = cos(qJ(1));
t77 = cos(qJ(2));
t34 = sin(pkin(5));
t38 = sin(qJ(2));
t73 = t34 * t38;
t39 = sin(qJ(1));
t72 = t34 * t39;
t71 = t34 * t40;
t68 = t78 * pkin(1) + pkin(7) * t72;
t67 = cos(pkin(5));
t63 = t34 * t78;
t62 = t34 * t77;
t58 = -pkin(1) * t39 + pkin(7) * t63;
t56 = t38 * t67;
t16 = t39 * t77 + t78 * t56;
t4 = t16 * t40 - t37 * t63;
t52 = t67 * t77;
t3 = t16 * t37 + t40 * t63;
t18 = -t39 * t56 + t78 * t77;
t17 = t78 * t38 + t39 * t52;
t15 = t38 * t39 - t78 * t52;
t14 = t67 * t37 + t38 * t71;
t13 = t37 * t73 - t67 * t40;
t8 = t18 * t40 + t37 * t72;
t7 = t18 * t37 - t39 * t71;
t2 = t17 * t29 + t30 * t8;
t1 = t17 * t30 - t29 * t8;
t5 = [(-t78 * mrSges(2,1) - m(3) * t68 - t18 * mrSges(3,1) - t2 * mrSges(6,1) - t1 * mrSges(6,2) + (-mrSges(3,3) * t34 + mrSges(2,2)) * t39 + t43 * t8 + t45 * t7 + t80 * t17 - t97 * (t18 * pkin(2) + t68)) * g(2) + (t39 * mrSges(2,1) + t78 * mrSges(2,2) - m(3) * t58 + t16 * mrSges(3,1) - mrSges(3,3) * t63 - t45 * t3 - t96 * t4 + (t50 - t80) * t15 + t97 * (t16 * pkin(2) - t58)) * g(1), (-t97 * (pkin(2) * t62 + pkin(8) * t73) + (t100 * t38 + t102 * t77) * t34) * g(3) + (t95 * t15 + t85 * t16) * g(2) + (t95 * t17 + t85 * t18) * g(1), (-t13 * t96 + t45 * t14) * g(3) + (-t3 * t96 + t45 * t4) * g(2) + (t45 * t8 - t7 * t96) * g(1), t84 * (-g(1) * t7 - g(2) * t3 - g(3) * t13), -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * ((t15 * t30 - t29 * t4) * mrSges(6,1) + (-t15 * t29 - t30 * t4) * mrSges(6,2)) - g(3) * ((-t14 * t29 - t30 * t62) * mrSges(6,1) + (-t14 * t30 + t29 * t62) * mrSges(6,2))];
taug = t5(:);
