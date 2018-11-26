% Calculate Gravitation load on the joints for
% S6RPRPRP3
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
% Datum: 2018-11-23 15:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:58:10
% EndTime: 2018-11-23 15:58:10
% DurationCPUTime: 0.66s
% Computational Cost: add. (443->81), mult. (435->92), div. (0->0), fcn. (396->10), ass. (0->46)
t77 = mrSges(6,1) + mrSges(7,1);
t76 = -mrSges(6,2) + mrSges(7,3);
t75 = m(6) + m(7);
t23 = sin(pkin(10));
t24 = cos(pkin(10));
t74 = m(5) * pkin(3) + t24 * mrSges(5,1) - t23 * mrSges(5,2) + mrSges(4,1);
t26 = sin(qJ(3));
t22 = qJ(1) + pkin(9);
t17 = sin(t22);
t19 = cos(t22);
t70 = g(1) * t19 + g(2) * t17;
t73 = t70 * t26;
t21 = pkin(10) + qJ(5);
t16 = sin(t21);
t18 = cos(t21);
t72 = t76 * t16 + t77 * t18;
t71 = mrSges(4,2) - mrSges(7,2) - mrSges(6,3);
t28 = cos(qJ(3));
t69 = -t71 * t26 + t74 * t28;
t15 = pkin(4) * t24 + pkin(3);
t11 = t28 * t15;
t44 = m(5) * qJ(4) + mrSges(5,3);
t25 = -pkin(8) - qJ(4);
t51 = t25 * t26;
t68 = -t44 * t26 - t69 - t75 * (t11 - t51);
t67 = -m(4) - m(5);
t64 = -mrSges(5,1) * t23 - mrSges(5,2) * t24 + mrSges(3,2) - mrSges(4,3);
t60 = m(7) * pkin(5) + t77;
t59 = m(7) * qJ(6) + t76;
t27 = sin(qJ(1));
t58 = pkin(1) * t27;
t57 = pkin(4) * t23;
t54 = g(3) * t26;
t29 = cos(qJ(1));
t20 = t29 * pkin(1);
t53 = t17 * t28;
t52 = t19 * t28;
t48 = m(5) + t75;
t47 = t19 * pkin(2) + t17 * pkin(7) + t20;
t46 = t19 * pkin(7) - t58;
t37 = pkin(5) * t18 + qJ(6) * t16;
t4 = t16 * t17 + t18 * t52;
t3 = t16 * t52 - t17 * t18;
t2 = -t19 * t16 + t18 * t53;
t1 = t16 * t53 + t18 * t19;
t5 = [(-m(3) * t20 - t29 * mrSges(2,1) + t27 * mrSges(2,2) + t67 * t47 - t75 * (t17 * t57 + t47) - t60 * t4 - t59 * t3 + t64 * t17 + (-mrSges(3,1) + t68) * t19) * g(2) + (m(3) * t58 + t27 * mrSges(2,1) + t29 * mrSges(2,2) + t67 * t46 - t75 * (t17 * t51 + t19 * t57 + t46) + t60 * t2 + t59 * t1 + t64 * t19 + (mrSges(3,1) + m(4) * pkin(2) - m(5) * (-qJ(4) * t26 - pkin(2)) + t26 * mrSges(5,3) - t75 * (-pkin(2) - t11) + t69) * t17) * g(1) (-m(3) - m(4) - t48) * g(3), t68 * g(3) + ((-m(7) * t37 - t72) * g(3) + t70 * (t75 * t25 - t44 + t71)) * t28 + (m(6) * t15 - m(7) * (-t15 - t37) + t72 + t74) * t73 (g(3) * t28 - t73) * t48 (t60 * t16 - t59 * t18) * t54 + (t60 * t1 - t59 * t2) * g(2) + (t60 * t3 - t59 * t4) * g(1) (-g(1) * t3 - g(2) * t1 - t16 * t54) * m(7)];
taug  = t5(:);
