% Calculate Gravitation load on the joints for
% S5RRPPR1
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:19
% EndTime: 2022-01-20 09:51:20
% DurationCPUTime: 0.21s
% Computational Cost: add. (242->50), mult. (156->50), div. (0->0), fcn. (114->10), ass. (0->26)
t25 = pkin(9) + qJ(5);
t19 = sin(t25);
t20 = cos(t25);
t46 = mrSges(6,1) * t20 - mrSges(6,2) * t19;
t45 = mrSges(4,2) - mrSges(6,3) - mrSges(5,3);
t28 = cos(pkin(9));
t44 = -t28 * mrSges(5,1) - mrSges(4,1) + sin(pkin(9)) * mrSges(5,2) - t46;
t43 = m(5) + m(6);
t26 = qJ(1) + qJ(2);
t23 = cos(t26);
t18 = pkin(2) * t23;
t38 = m(4) + t43;
t21 = pkin(8) + t26;
t15 = sin(t21);
t16 = cos(t21);
t37 = t16 * pkin(3) + t15 * qJ(4) + t18;
t17 = pkin(4) * t28 + pkin(3);
t29 = -pkin(7) - qJ(4);
t36 = -t15 * t29 + t16 * t17 + t18;
t22 = sin(t26);
t33 = -t23 * mrSges(3,1) + t22 * mrSges(3,2) + t45 * t15 + t44 * t16;
t32 = mrSges(3,2) * t23 + (t38 * pkin(2) + mrSges(3,1)) * t22 + (m(5) * pkin(3) + m(6) * t17 - t44) * t15 + (-m(5) * qJ(4) + m(6) * t29 + t45) * t16;
t31 = cos(qJ(1));
t30 = sin(qJ(1));
t24 = t31 * pkin(1);
t1 = [(t30 * mrSges(2,2) - m(4) * (t18 + t24) - m(5) * (t24 + t37) - m(6) * (t24 + t36) + (-m(3) * pkin(1) - mrSges(2,1)) * t31 + t33) * g(2) + (mrSges(2,2) * t31 + (mrSges(2,1) + (m(3) + t38) * pkin(1)) * t30 + t32) * g(1), (-m(4) * t18 - m(5) * t37 - m(6) * t36 + t33) * g(2) + t32 * g(1), -t38 * g(3), t43 * (-g(1) * t15 + g(2) * t16), -g(3) * t46 + (g(1) * t16 + g(2) * t15) * (mrSges(6,1) * t19 + mrSges(6,2) * t20)];
taug = t1(:);
