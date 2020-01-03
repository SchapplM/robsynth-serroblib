% Calculate Gravitation load on the joints for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:21
% EndTime: 2020-01-03 11:57:21
% DurationCPUTime: 0.22s
% Computational Cost: add. (280->63), mult. (198->78), div. (0->0), fcn. (168->10), ass. (0->37)
t37 = sin(pkin(9));
t60 = t37 * mrSges(5,2) - mrSges(4,1);
t59 = mrSges(4,2) - mrSges(5,3);
t58 = m(5) + m(6);
t36 = qJ(1) + qJ(2);
t32 = sin(t36);
t29 = pkin(2) * t32;
t33 = cos(t36);
t30 = pkin(2) * t33;
t31 = pkin(8) + t36;
t27 = sin(t31);
t57 = t27 * t37;
t38 = cos(pkin(9));
t56 = t27 * t38;
t28 = cos(t31);
t55 = t28 * t37;
t54 = t28 * t38;
t39 = sin(qJ(5));
t52 = t38 * t39;
t41 = cos(qJ(5));
t51 = t38 * t41;
t50 = t28 * pkin(3) + t27 * qJ(4) + t30;
t49 = -m(3) * pkin(1) - mrSges(2,1);
t48 = t27 * pkin(3) - qJ(4) * t28 + t29;
t47 = pkin(4) * t54 + pkin(7) * t55 + t50;
t45 = pkin(4) * t56 + pkin(7) * t57 + t48;
t5 = -t27 * t52 - t28 * t41;
t6 = t27 * t51 - t28 * t39;
t44 = -t32 * mrSges(3,1) - mrSges(5,1) * t56 - t6 * mrSges(6,1) - t33 * mrSges(3,2) - t5 * mrSges(6,2) - mrSges(6,3) * t57 + t60 * t27 - t59 * t28;
t7 = -t27 * t41 + t28 * t52;
t8 = t27 * t39 + t28 * t51;
t43 = -t33 * mrSges(3,1) - mrSges(5,1) * t54 - t8 * mrSges(6,1) + t32 * mrSges(3,2) + t7 * mrSges(6,2) - mrSges(6,3) * t55 + t59 * t27 + t60 * t28;
t42 = cos(qJ(1));
t40 = sin(qJ(1));
t35 = t42 * pkin(1);
t34 = t40 * pkin(1);
t1 = [(t44 + t49 * t40 - mrSges(2,2) * t42 - m(4) * (t29 + t34) - m(5) * (t34 + t48) - m(6) * (t34 + t45)) * g(3) + (t43 + t49 * t42 - m(4) * (t30 + t35) - m(5) * (t35 + t50) - m(6) * (t35 + t47) + t40 * mrSges(2,2)) * g(2), (-m(4) * t29 - m(5) * t48 - m(6) * t45 + t44) * g(3) + (-m(4) * t30 - m(5) * t50 - m(6) * t47 + t43) * g(2), (-m(4) - t58) * g(1), t58 * (g(2) * t28 + g(3) * t27), -g(2) * (mrSges(6,1) * t5 - mrSges(6,2) * t6) - g(3) * (mrSges(6,1) * t7 + mrSges(6,2) * t8) - g(1) * (-mrSges(6,1) * t39 - mrSges(6,2) * t41) * t37];
taug = t1(:);
