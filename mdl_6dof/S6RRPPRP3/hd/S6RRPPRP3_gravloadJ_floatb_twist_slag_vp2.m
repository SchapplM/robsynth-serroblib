% Calculate Gravitation load on the joints for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2018-11-23 16:46
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:46:03
% EndTime: 2018-11-23 16:46:04
% DurationCPUTime: 0.80s
% Computational Cost: add. (231->89), mult. (478->94), div. (0->0), fcn. (415->6), ass. (0->45)
t72 = -mrSges(6,1) - mrSges(7,1);
t71 = mrSges(6,2) + mrSges(7,2);
t82 = mrSges(3,1) + mrSges(4,1);
t81 = mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t80 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t23 = cos(qJ(5));
t12 = pkin(5) * t23 + pkin(4);
t20 = sin(qJ(5));
t79 = -m(6) * pkin(4) - m(7) * t12 + t71 * t20 + t72 * t23;
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t69 = g(1) * t25 + g(2) * t22;
t19 = -qJ(6) - pkin(8);
t62 = -pkin(2) - pkin(3);
t78 = -m(5) * t62 - m(6) * (-pkin(8) + t62) - m(7) * (t19 + t62) - t80;
t21 = sin(qJ(2));
t24 = cos(qJ(2));
t75 = t81 * t21 + t82 * t24;
t63 = m(7) * pkin(5);
t45 = m(5) + m(6) + m(7);
t73 = m(4) + t45;
t68 = t63 - t72;
t67 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t64 = t80 * t24 - t75;
t59 = g(3) * t24;
t16 = t24 * pkin(2);
t58 = t24 * pkin(8);
t57 = t19 * t24;
t56 = t20 * t25;
t54 = t22 * t20;
t53 = t22 * t23;
t50 = t24 * t25;
t49 = t25 * t23;
t13 = t21 * qJ(3);
t48 = t16 + t13;
t47 = t25 * pkin(1) + t22 * pkin(7);
t44 = t24 * pkin(3) + t48;
t42 = pkin(2) * t50 + t25 * t13 + t47;
t41 = -pkin(1) - t13;
t3 = -t21 * t56 - t53;
t1 = t21 * t54 - t49;
t17 = t25 * pkin(7);
t4 = t21 * t49 - t54;
t2 = -t21 * t53 - t56;
t5 = [(t54 * t63 - m(3) * t47 - m(4) * t42 + t72 * t4 - t45 * (pkin(3) * t50 - t22 * qJ(4) + t42) - t71 * t3 + t67 * t22 + (-mrSges(2,1) - m(6) * (pkin(4) * t21 + t58) - m(7) * (t12 * t21 - t57) + t64) * t25) * g(2) + (t56 * t63 - t45 * (-qJ(4) * t25 + t17) + t72 * t2 + (-m(3) - m(4)) * t17 - t71 * t1 + t67 * t25 + (mrSges(2,1) + m(3) * pkin(1) - m(4) * (t41 - t16) - m(5) * t41 - m(6) * (-pkin(1) + (-pkin(4) - qJ(3)) * t21) - m(7) * (-pkin(1) + (-qJ(3) - t12) * t21) + t78 * t24 + t75) * t22) * g(1) (-m(4) * t48 - m(5) * t44 - m(6) * (t44 + t58) - m(7) * (t44 - t57) + t79 * t21 + t64) * g(3) + ((m(4) * pkin(2) + t78 + t82) * t21 + (-qJ(3) * t73 + t79 - t81) * t24) * t69 (-t69 * t21 + t59) * t73 (g(1) * t22 - g(2) * t25) * t45 (-t68 * t20 - t71 * t23) * t59 + (t68 * t1 - t71 * t2) * g(2) + (-t68 * t3 + t71 * t4) * g(1) (-g(3) * t21 - t69 * t24) * m(7)];
taug  = t5(:);
