% Calculate Gravitation load on the joints for
% S5RRRPP4
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
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:52
% EndTime: 2019-12-31 20:54:54
% DurationCPUTime: 0.44s
% Computational Cost: add. (279->69), mult. (269->71), div. (0->0), fcn. (208->8), ass. (0->36)
t24 = sin(qJ(2));
t26 = cos(qJ(2));
t23 = qJ(2) + qJ(3);
t18 = pkin(8) + t23;
t14 = sin(t18);
t15 = cos(t18);
t19 = sin(t23);
t20 = cos(t23);
t66 = -t20 * mrSges(4,1) + t19 * mrSges(4,2) + (-mrSges(5,1) - mrSges(6,1)) * t15 + (mrSges(5,2) - mrSges(6,3)) * t14;
t73 = t26 * mrSges(3,1) - t24 * mrSges(3,2) - t66;
t72 = mrSges(4,1) * t19 + mrSges(5,1) * t14 + mrSges(4,2) * t20 + mrSges(5,2) * t15;
t71 = (m(6) * qJ(5) + mrSges(6,3)) * t15;
t25 = sin(qJ(1));
t27 = cos(qJ(1));
t67 = g(1) * t27 + g(2) * t25;
t70 = m(5) + m(6);
t36 = t15 * pkin(4) + t14 * qJ(5);
t21 = t26 * pkin(2);
t64 = mrSges(2,1) + m(4) * (t21 + pkin(1)) + m(3) * pkin(1) + t73;
t28 = -pkin(7) - pkin(6);
t63 = -m(3) * pkin(6) + m(4) * t28 + mrSges(2,2) - mrSges(6,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t62 = pkin(3) * t19;
t16 = pkin(3) * t20;
t61 = pkin(4) * t14;
t58 = t24 * pkin(2);
t56 = t14 * mrSges(6,1);
t53 = t16 + t21;
t51 = t16 + t36;
t50 = t71 * t25;
t49 = t71 * t27;
t3 = -t58 - t62;
t32 = m(6) * (t3 - t61) - t56;
t29 = m(6) * (-t61 - t62) - t56;
t22 = -qJ(4) + t28;
t2 = pkin(1) + t53;
t1 = [(-t70 * (t27 * t2 - t25 * t22) + (-m(6) * t36 - t64) * t27 + t63 * t25) * g(2) + ((t70 * t22 + t63) * t27 + (m(5) * t2 - m(6) * (-t2 - t36) + t64) * t25) * g(1), -g(1) * (t32 * t27 + t49) - g(2) * (t32 * t25 + t50) + (-m(4) * t21 - m(5) * t53 - m(6) * (t21 + t51) - t73) * g(3) + t67 * (m(4) * t58 - m(5) * t3 + mrSges(3,1) * t24 + mrSges(3,2) * t26 + t72), -g(1) * (t29 * t27 + t49) - g(2) * (t29 * t25 + t50) + (-m(5) * t16 - m(6) * t51 + t66) * g(3) + (m(5) * t62 + t72) * t67, t70 * (-g(1) * t25 + g(2) * t27), (g(3) * t15 - t67 * t14) * m(6)];
taug = t1(:);
