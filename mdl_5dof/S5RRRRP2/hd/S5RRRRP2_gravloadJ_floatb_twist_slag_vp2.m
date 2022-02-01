% Calculate Gravitation load on the joints for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:48:45
% EndTime: 2022-01-20 11:48:46
% DurationCPUTime: 0.31s
% Computational Cost: add. (289->60), mult. (244->61), div. (0->0), fcn. (186->8), ass. (0->33)
t33 = qJ(3) + qJ(4);
t26 = sin(t33);
t28 = cos(t33);
t52 = mrSges(5,2) + mrSges(6,2);
t53 = mrSges(5,1) + mrSges(6,1);
t68 = t52 * t26 - t53 * t28;
t35 = sin(qJ(3));
t67 = t35 * mrSges(4,2) + t68;
t37 = cos(qJ(3));
t66 = -t37 * mrSges(4,1) - mrSges(3,1) + t67;
t62 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t34 = qJ(1) + qJ(2);
t27 = sin(t34);
t29 = cos(t34);
t60 = g(1) * t29 + g(2) * t27;
t39 = -pkin(8) - pkin(7);
t24 = pkin(4) * t28;
t51 = t29 * pkin(2) + t27 * pkin(7);
t30 = t37 * pkin(3);
t50 = t24 + t30;
t49 = m(5) * pkin(3) + mrSges(4,1);
t2 = pkin(2) + t50;
t32 = qJ(5) - t39;
t48 = t29 * t2 + t27 * t32;
t25 = t30 + pkin(2);
t47 = t29 * t25 - t27 * t39;
t45 = t52 * t28;
t42 = t62 * t27 + t66 * t29;
t40 = (-m(4) * pkin(7) + m(5) * t39 - m(6) * t32 + t62) * t29 + (m(4) * pkin(2) + m(5) * t25 + m(6) * t2 - t66) * t27;
t38 = cos(qJ(1));
t36 = sin(qJ(1));
t31 = t38 * pkin(1);
t1 = [(t36 * mrSges(2,2) - m(4) * (t31 + t51) - m(5) * (t31 + t47) - m(6) * (t31 + t48) + (-m(3) * pkin(1) - mrSges(2,1)) * t38 + t42) * g(2) + (t38 * mrSges(2,2) + (mrSges(2,1) + (m(3) + m(4) + m(5) + m(6)) * pkin(1)) * t36 + t40) * g(1), (-m(4) * t51 - m(5) * t47 - m(6) * t48 + t42) * g(2) + t40 * g(1), (-m(6) * t50 - t49 * t37 + t67) * g(3) + t60 * (-m(6) * (-t35 * pkin(3) - pkin(4) * t26) + mrSges(4,2) * t37 + t53 * t26 + t49 * t35 + t45), (-m(6) * t24 + t68) * g(3) + t60 * (t45 + (m(6) * pkin(4) + t53) * t26), (-g(1) * t27 + g(2) * t29) * m(6)];
taug = t1(:);
