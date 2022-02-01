% Calculate Gravitation load on the joints for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:18:48
% EndTime: 2022-01-23 09:18:49
% DurationCPUTime: 0.20s
% Computational Cost: add. (223->47), mult. (143->49), div. (0->0), fcn. (104->10), ass. (0->26)
t24 = pkin(9) + qJ(5);
t18 = sin(t24);
t20 = cos(t24);
t47 = mrSges(6,1) * t20 - mrSges(6,2) * t18;
t46 = mrSges(4,2) - mrSges(6,3) - mrSges(5,3);
t27 = cos(pkin(9));
t45 = -t27 * mrSges(5,1) - mrSges(4,1) + sin(pkin(9)) * mrSges(5,2) - t47;
t44 = m(5) + m(6);
t25 = qJ(1) + pkin(8);
t22 = qJ(3) + t25;
t15 = sin(t22);
t16 = cos(t22);
t43 = t16 * pkin(3) + t15 * qJ(4);
t21 = cos(t25);
t30 = cos(qJ(1));
t38 = t30 * pkin(1) + pkin(2) * t21;
t37 = m(4) + t44;
t36 = m(3) + t37;
t17 = pkin(4) * t27 + pkin(3);
t28 = -pkin(7) - qJ(4);
t35 = -t15 * t28 + t16 * t17;
t32 = t46 * t15 + t45 * t16;
t31 = (m(5) * pkin(3) + m(6) * t17 - t45) * t15 + (-m(5) * qJ(4) + m(6) * t28 + t46) * t16;
t29 = sin(qJ(1));
t19 = sin(t25);
t1 = [(t29 * mrSges(2,2) - t21 * mrSges(3,1) + t19 * mrSges(3,2) - m(4) * t38 - m(5) * (t38 + t43) - m(6) * (t35 + t38) + (-m(3) * pkin(1) - mrSges(2,1)) * t30 + t32) * g(2) + (mrSges(2,2) * t30 + mrSges(3,2) * t21 + (t37 * pkin(2) + mrSges(3,1)) * t19 + (t36 * pkin(1) + mrSges(2,1)) * t29 + t31) * g(1), -t36 * g(3), (-m(5) * t43 - m(6) * t35 + t32) * g(2) + t31 * g(1), t44 * (-g(1) * t15 + g(2) * t16), -g(3) * t47 + (g(1) * t16 + g(2) * t15) * (mrSges(6,1) * t18 + mrSges(6,2) * t20)];
taug = t1(:);
