% Calculate Gravitation load on the joints for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:59:16
% EndTime: 2018-11-23 15:59:17
% DurationCPUTime: 0.68s
% Computational Cost: add. (433->89), mult. (458->95), div. (0->0), fcn. (414->10), ass. (0->49)
t77 = mrSges(6,1) + mrSges(7,1);
t76 = -mrSges(6,2) + mrSges(7,3);
t75 = m(6) + m(7);
t20 = pkin(9) + qJ(3);
t16 = sin(t20);
t27 = sin(qJ(1));
t28 = cos(qJ(1));
t69 = g(1) * t28 + g(2) * t27;
t74 = t69 * t16;
t73 = -mrSges(6,3) - mrSges(7,2);
t19 = pkin(10) + qJ(5);
t15 = sin(t19);
t17 = cos(t19);
t72 = t76 * t15 + t77 * t17;
t23 = cos(pkin(10));
t13 = pkin(4) * t23 + pkin(3);
t18 = cos(t20);
t25 = -pkin(8) - qJ(4);
t71 = t18 * t13 - t16 * t25;
t24 = cos(pkin(9));
t70 = -mrSges(2,1) - m(3) * pkin(1) - t24 * mrSges(3,1) + sin(pkin(9)) * mrSges(3,2);
t21 = sin(pkin(10));
t35 = m(5) * pkin(3) + t23 * mrSges(5,1) - t21 * mrSges(5,2);
t67 = t35 * t18;
t41 = t18 * mrSges(4,1) - t16 * mrSges(4,2);
t66 = t73 * t16 - t41;
t26 = -pkin(7) - qJ(2);
t65 = -m(3) * qJ(2) + m(5) * t26 - t21 * mrSges(5,1) - t23 * mrSges(5,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t63 = m(7) * pkin(5) + t77;
t62 = m(7) * qJ(6) + t76;
t61 = pkin(4) * t21;
t58 = g(3) * t16;
t54 = t16 * t28;
t53 = t18 * t28;
t52 = t27 * t15;
t51 = t27 * t17;
t50 = t28 * t15;
t49 = m(5) + t75;
t14 = pkin(2) * t24 + pkin(1);
t8 = t28 * t14;
t47 = -t27 * t26 + t8;
t45 = m(5) * qJ(4) + mrSges(5,3);
t37 = pkin(5) * t17 + qJ(6) * t15;
t32 = t45 * t16 + t67;
t4 = t17 * t53 + t52;
t3 = t18 * t50 - t51;
t2 = t18 * t51 - t50;
t1 = t17 * t28 + t18 * t52;
t5 = [(-m(4) * t47 - m(5) * t8 + t73 * t54 - t75 * (t13 * t53 - t25 * t54 + t27 * t61 + t47) - t63 * t4 - t62 * t3 + (-t41 - t32 + t70) * t28 + t65 * t27) * g(2) + (t62 * t1 + t63 * t2 + (m(4) * t14 - m(5) * (-qJ(4) * t16 - t14) + t16 * mrSges(5,3) + t67 - t66 - (-t14 - t71) * t75 - t70) * t27 + (-t75 * (-t26 + t61) + m(4) * t26 + t65) * t28) * g(1) (-g(1) * t27 + g(2) * t28) * (m(3) + m(4) + t49) (-t71 * t75 - t32 + t66) * g(3) + ((-m(7) * t37 - t72) * g(3) + t69 * (t75 * t25 + mrSges(4,2) - t45 + t73)) * t18 + (mrSges(4,1) + t35 + m(6) * t13 - m(7) * (-t13 - t37) + t72) * t74 (g(3) * t18 - t74) * t49 (t63 * t15 - t62 * t17) * t58 + (t63 * t1 - t62 * t2) * g(2) + (t63 * t3 - t62 * t4) * g(1) (-g(1) * t3 - g(2) * t1 - t15 * t58) * m(7)];
taug  = t5(:);
