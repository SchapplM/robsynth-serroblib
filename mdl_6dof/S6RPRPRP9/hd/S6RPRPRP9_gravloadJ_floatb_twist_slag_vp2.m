% Calculate Gravitation load on the joints for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2018-11-23 16:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:01:57
% EndTime: 2018-11-23 16:01:57
% DurationCPUTime: 0.70s
% Computational Cost: add. (298->86), mult. (452->100), div. (0->0), fcn. (408->8), ass. (0->47)
t23 = -pkin(8) - qJ(4);
t72 = m(6) + m(7);
t75 = t72 * t23;
t74 = mrSges(6,1) + mrSges(7,1);
t73 = -mrSges(6,2) + mrSges(7,3);
t71 = -m(4) - m(5);
t70 = mrSges(6,3) + mrSges(7,2);
t20 = pkin(9) + qJ(5);
t14 = sin(t20);
t15 = cos(t20);
t69 = t73 * t14 + t74 * t15;
t25 = sin(qJ(1));
t27 = cos(qJ(1));
t62 = -g(1) * t25 + g(2) * t27;
t21 = sin(pkin(9));
t22 = cos(pkin(9));
t68 = -t21 * mrSges(5,1) - t22 * mrSges(5,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t24 = sin(qJ(3));
t26 = cos(qJ(3));
t31 = m(5) * pkin(3) + mrSges(5,1) * t22 - mrSges(5,2) * t21;
t38 = t24 * mrSges(4,1) + t26 * mrSges(4,2);
t41 = m(5) * qJ(4) + mrSges(5,3);
t67 = -t31 * t24 + t41 * t26 + mrSges(2,2) - mrSges(3,3) - t38;
t61 = m(7) * pkin(5) + t74;
t60 = m(7) * qJ(6) + t73;
t13 = pkin(4) * t22 + pkin(3);
t34 = pkin(5) * t15 + qJ(6) * t14;
t65 = -m(7) * (-t13 - t34) + m(6) * t13 + t69;
t63 = -t70 + t75;
t58 = -pkin(1) - pkin(7);
t57 = pkin(4) * t21;
t54 = g(3) * t26;
t53 = t24 * t27;
t52 = t25 * t14;
t51 = t25 * t15;
t50 = t25 * t26;
t49 = t26 * t27;
t48 = t27 * t15;
t47 = t27 * pkin(1) + t25 * qJ(2);
t46 = m(5) + t72;
t44 = t27 * pkin(7) + t47;
t17 = t27 * qJ(2);
t4 = t24 * t48 - t52;
t3 = t14 * t53 + t51;
t2 = t14 * t27 + t24 * t51;
t1 = t24 * t52 - t48;
t5 = [(-m(3) * t47 + t70 * t50 + t71 * t44 - t72 * (t25 * t24 * t13 + t23 * t50 + t27 * t57 + t44) - t61 * t2 - t60 * t1 + t68 * t27 + t67 * t25) * g(2) + (-t72 * ((-t57 + t58) * t25 + t13 * t53 + t23 * t49 + t17) + t70 * t49 - t61 * t4 - t60 * t3 + (-m(3) + t71) * t17 + t67 * t27 + (m(3) * pkin(1) + t71 * t58 - t68) * t25) * g(1), t62 * (m(3) + m(4) + t46) (-t53 * t75 + (t70 * t24 + t65 * t26) * t27) * g(2) + (-t72 * t13 * t50 + ((-m(7) * t34 - t69) * t26 + t63 * t24) * t25) * g(1) + (t38 + (-t41 + t63) * t26 + (t31 + t65) * t24) * g(3) + ((mrSges(4,1) + t31) * t26 + (-mrSges(4,2) + t41) * t24) * t62 (-g(3) * t24 - t62 * t26) * t46 (t61 * t14 - t60 * t15) * t54 + (-t61 * t3 + t60 * t4) * g(2) + (t61 * t1 - t60 * t2) * g(1) (-g(1) * t1 + g(2) * t3 - t14 * t54) * m(7)];
taug  = t5(:);
