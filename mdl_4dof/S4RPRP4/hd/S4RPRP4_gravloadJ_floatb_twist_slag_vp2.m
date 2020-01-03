% Calculate Gravitation load on the joints for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:38
% EndTime: 2019-12-31 16:43:38
% DurationCPUTime: 0.20s
% Computational Cost: add. (100->34), mult. (114->38), div. (0->0), fcn. (85->6), ass. (0->18)
t29 = mrSges(4,1) + mrSges(5,1);
t28 = -mrSges(4,2) + mrSges(5,3);
t10 = cos(qJ(3));
t8 = sin(qJ(3));
t27 = t29 * t10 + t28 * t8;
t26 = -m(4) - m(5);
t14 = pkin(3) * t10 + qJ(4) * t8;
t25 = -m(5) * t14 - t27;
t7 = qJ(1) + pkin(6);
t4 = sin(t7);
t5 = cos(t7);
t24 = g(1) * t5 + g(2) * t4;
t23 = mrSges(3,2) - mrSges(4,3) - mrSges(5,2);
t9 = sin(qJ(1));
t22 = pkin(1) * t9;
t11 = cos(qJ(1));
t6 = t11 * pkin(1);
t1 = [(-m(3) * t6 - mrSges(2,1) * t11 + t9 * mrSges(2,2) + t26 * (t5 * pkin(2) + t4 * pkin(5) + t6) + (-mrSges(3,1) + t25) * t5 + t23 * t4) * g(2) + (m(3) * t22 + t9 * mrSges(2,1) + mrSges(2,2) * t11 + t26 * (t5 * pkin(5) - t22) + t23 * t5 + (mrSges(3,1) + m(4) * pkin(2) - m(5) * (-pkin(2) - t14) + t27) * t4) * g(1), (-m(3) + t26) * g(3), t25 * g(3) + ((m(5) * pkin(3) + t29) * t8 + (-m(5) * qJ(4) - t28) * t10) * t24, (g(3) * t10 - t24 * t8) * m(5)];
taug = t1(:);
