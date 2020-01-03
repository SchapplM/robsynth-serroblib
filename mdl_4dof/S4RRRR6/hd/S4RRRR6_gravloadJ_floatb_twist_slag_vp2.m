% Calculate Gravitation load on the joints for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:22
% EndTime: 2019-12-31 17:29:24
% DurationCPUTime: 0.50s
% Computational Cost: add. (254->73), mult. (622->110), div. (0->0), fcn. (724->10), ass. (0->43)
t59 = m(4) + m(5);
t67 = t59 * pkin(7);
t30 = sin(qJ(4));
t34 = cos(qJ(4));
t64 = m(5) * pkin(3) + mrSges(5,1) * t34 - mrSges(5,2) * t30 + mrSges(4,1);
t44 = -m(5) * pkin(8) + mrSges(4,2) - mrSges(5,3);
t40 = t30 * mrSges(5,1) + t34 * mrSges(5,2);
t53 = mrSges(3,2) - mrSges(4,3);
t66 = -t40 + t53;
t31 = sin(qJ(3));
t35 = cos(qJ(3));
t63 = -t44 * t31 + t64 * t35 + mrSges(3,1);
t65 = pkin(2) * t59 + t63;
t61 = t53 - t67;
t60 = t66 - t67;
t58 = cos(qJ(1));
t29 = sin(pkin(4));
t32 = sin(qJ(2));
t57 = t29 * t32;
t33 = sin(qJ(1));
t56 = t29 * t33;
t55 = t29 * t35;
t36 = cos(qJ(2));
t54 = t29 * t36;
t51 = t58 * pkin(1) + pkin(6) * t56;
t50 = cos(pkin(4));
t45 = t33 * t50;
t18 = -t32 * t45 + t58 * t36;
t49 = t18 * pkin(2) + t51;
t47 = t29 * t58;
t46 = -pkin(1) * t33 + pkin(6) * t47;
t42 = t50 * t58;
t16 = t32 * t42 + t33 * t36;
t4 = t16 * t35 - t31 * t47;
t3 = -t16 * t31 - t35 * t47;
t17 = t58 * t32 + t36 * t45;
t15 = t33 * t32 - t36 * t42;
t14 = t50 * t31 + t32 * t55;
t8 = t18 * t35 + t31 * t56;
t7 = t18 * t31 - t33 * t55;
t2 = t17 * t30 + t34 * t8;
t1 = t17 * t34 - t30 * t8;
t5 = [(-t58 * mrSges(2,1) - m(3) * t51 - t18 * mrSges(3,1) - m(4) * t49 - t8 * mrSges(4,1) - m(5) * (pkin(3) * t8 + t49) - t2 * mrSges(5,1) - t1 * mrSges(5,2) + t44 * t7 + (-mrSges(3,3) * t29 + mrSges(2,2)) * t33 + t61 * t17) * g(2) + (t33 * mrSges(2,1) + t58 * mrSges(2,2) - m(3) * t46 + t16 * mrSges(3,1) - mrSges(3,3) * t47 + t44 * t3 + t64 * t4 + (t40 - t61) * t15 + t59 * (t16 * pkin(2) - t46)) * g(1), (-t59 * (pkin(2) * t54 + pkin(7) * t57) + (t66 * t32 - t63 * t36) * t29) * g(3) + (t65 * t15 + t60 * t16) * g(2) + (t65 * t17 + t60 * t18) * g(1), (t44 * t14 - t64 * (-t31 * t57 + t50 * t35)) * g(3) + (-t64 * t3 + t44 * t4) * g(2) + (t44 * t8 + t64 * t7) * g(1), -g(1) * (mrSges(5,1) * t1 - mrSges(5,2) * t2) - g(2) * ((t15 * t34 - t30 * t4) * mrSges(5,1) + (-t15 * t30 - t34 * t4) * mrSges(5,2)) - g(3) * ((-t14 * t30 - t34 * t54) * mrSges(5,1) + (-t14 * t34 + t30 * t54) * mrSges(5,2))];
taug = t5(:);
