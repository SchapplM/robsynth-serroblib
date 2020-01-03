% Calculate Gravitation load on the joints for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:24
% EndTime: 2019-12-31 20:49:25
% DurationCPUTime: 0.37s
% Computational Cost: add. (283->64), mult. (241->67), div. (0->0), fcn. (187->8), ass. (0->34)
t27 = qJ(3) + pkin(8);
t21 = sin(t27);
t30 = sin(qJ(3));
t32 = cos(qJ(3));
t73 = mrSges(5,2) - mrSges(6,3);
t77 = -t32 * mrSges(4,1) + t30 * mrSges(4,2) + t73 * t21;
t74 = -mrSges(5,1) - mrSges(6,1);
t75 = m(5) + m(6);
t22 = cos(t27);
t72 = t74 * t22 + t77;
t71 = -mrSges(6,2) - mrSges(5,3) - mrSges(4,3) + mrSges(3,2);
t28 = qJ(1) + qJ(2);
t23 = sin(t28);
t24 = cos(t28);
t68 = g(1) * t24 + g(2) * t23;
t25 = t32 * pkin(3);
t20 = t25 + pkin(2);
t52 = t21 * qJ(5);
t42 = -t22 * pkin(4) - t52;
t66 = t71 * t24 + (-m(6) * (-t20 + t42) + mrSges(3,1) - t72) * t23;
t57 = t22 * t24;
t65 = t71 * t23 + t74 * t57 + (-mrSges(3,1) + t77) * t24;
t31 = sin(qJ(1));
t61 = t31 * pkin(1);
t33 = cos(qJ(1));
t26 = t33 * pkin(1);
t29 = -qJ(4) - pkin(7);
t56 = t24 * t29;
t53 = t24 * pkin(2) + t23 * pkin(7);
t50 = -t23 * pkin(2) + t24 * pkin(7);
t49 = t24 * t20 - t23 * t29;
t46 = pkin(4) * t57 + t24 * t52 + t49;
t41 = -t23 * t20 - t56;
t1 = [(-t33 * mrSges(2,1) + t31 * mrSges(2,2) - m(3) * t26 - m(4) * (t26 + t53) - m(5) * (t26 + t49) - m(6) * (t26 + t46) + t65) * g(2) + (t31 * mrSges(2,1) + t33 * mrSges(2,2) + m(3) * t61 - m(4) * (t50 - t61) - m(5) * (t41 - t61) - m(6) * (-t56 - t61) + t66) * g(1), (-m(4) * t53 - m(5) * t49 - m(6) * t46 + t65) * g(2) + (-m(4) * t50 - m(5) * t41 + m(6) * t56 + t66) * g(1), (-m(5) * t25 - m(6) * (t25 - t42) + t72) * g(3) + t68 * (mrSges(4,2) * t32 + (-m(6) * qJ(5) + t73) * t22 + (m(6) * pkin(4) - t74) * t21 + (t75 * pkin(3) + mrSges(4,1)) * t30), t75 * (-g(1) * t23 + g(2) * t24), (g(3) * t22 - t68 * t21) * m(6)];
taug = t1(:);
