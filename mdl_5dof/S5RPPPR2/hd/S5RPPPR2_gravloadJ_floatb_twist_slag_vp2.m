% Calculate Gravitation load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:04
% EndTime: 2022-01-23 08:59:06
% DurationCPUTime: 0.47s
% Computational Cost: add. (166->83), mult. (339->121), div. (0->0), fcn. (353->10), ass. (0->43)
t24 = sin(pkin(9));
t27 = cos(pkin(9));
t18 = t27 * pkin(4) + t24 * pkin(6) + pkin(3);
t25 = sin(pkin(8));
t28 = cos(pkin(8));
t38 = qJ(4) * t28 - qJ(2);
t62 = qJ(2) * (-m(3) - m(4)) + m(5) * (-t25 * pkin(3) + t38) + m(6) * (-t18 * t25 + t38) + mrSges(2,2) - mrSges(3,3);
t61 = m(5) + m(6);
t32 = cos(qJ(5));
t60 = mrSges(6,1) * t32;
t59 = -mrSges(4,2) + mrSges(5,3);
t58 = mrSges(5,2) - mrSges(6,3);
t30 = sin(qJ(5));
t50 = t25 * t30;
t10 = t27 * t50 + t32 * t28;
t29 = cos(pkin(7));
t49 = t25 * t32;
t26 = sin(pkin(7));
t48 = t27 * t28;
t9 = t26 * t24 + t29 * t48;
t2 = t29 * t49 - t9 * t30;
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t57 = -t31 * t10 + t2 * t33;
t51 = t25 * qJ(4) + pkin(2);
t52 = t26 * qJ(3) + pkin(1);
t54 = -m(4) * (pkin(2) * t29 + t52) - m(5) * ((t28 * pkin(3) + t51) * t29 + t52) - m(6) * ((t18 * t28 + t51) * t29 + pkin(1) + (t24 * pkin(4) - t27 * pkin(6) + qJ(3)) * t26) - mrSges(2,1) - m(3) * pkin(1) - t29 * mrSges(3,1) + t26 * mrSges(3,2);
t46 = t31 * t25;
t45 = t31 * t26;
t44 = t31 * t29;
t43 = t33 * t25;
t42 = t33 * t26;
t41 = t33 * t29;
t40 = m(4) + t61;
t37 = t29 * t50 + t9 * t32;
t15 = t28 * t41 + t46;
t14 = t25 * t41 - t31 * t28;
t13 = -t28 * t44 + t43;
t12 = t25 * t44 + t33 * t28;
t11 = t27 * t49 - t30 * t28;
t8 = -t29 * t24 + t26 * t48;
t7 = t33 * t11;
t1 = [(-t15 * mrSges(4,1) - mrSges(4,3) * t42 - (t15 * t27 + t24 * t42) * mrSges(5,1) - t27 * t46 * t60 - t57 * mrSges(6,2) + t58 * (t15 * t24 - t27 * t42) + (-t9 * t60 + t54) * t33 + (-t30 * mrSges(6,1) - t59) * t14 + t62 * t31) * g(2) + (-t13 * mrSges(4,1) + mrSges(4,3) * t45 - (t13 * t27 - t24 * t45) * mrSges(5,1) - t7 * mrSges(6,1) + t58 * (t13 * t24 + t27 * t45) + (t37 * mrSges(6,1) + t2 * mrSges(6,2) - t54) * t31 + t59 * t12 + (t10 * mrSges(6,2) + t62) * t33) * g(1), (-g(1) * t31 + g(2) * t33) * (m(3) + t40), (t29 * g(3) + (-g(1) * t33 - g(2) * t31) * t26) * t40, t61 * (-g(3) * t25 * t26 - g(1) * t14 - g(2) * t12), -g(1) * (t57 * mrSges(6,1) + (-t31 * t11 - t33 * t37) * mrSges(6,2)) - g(2) * ((-(-t27 * t43 + t9 * t31) * t30 + t32 * t12) * mrSges(6,1) + (-t37 * t31 + t7) * mrSges(6,2)) - g(3) * ((t26 * t49 - t30 * t8) * mrSges(6,1) + (-t26 * t50 - t8 * t32) * mrSges(6,2))];
taug = t1(:);
