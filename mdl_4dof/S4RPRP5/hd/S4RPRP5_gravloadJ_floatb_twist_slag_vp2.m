% Calculate Gravitation load on the joints for
% S4RPRP5
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:46
% EndTime: 2019-12-31 16:44:46
% DurationCPUTime: 0.18s
% Computational Cost: add. (102->38), mult. (131->38), div. (0->0), fcn. (99->6), ass. (0->17)
t26 = m(4) + m(5);
t10 = cos(qJ(1));
t9 = sin(qJ(1));
t25 = g(1) * t10 + g(2) * t9;
t5 = pkin(6) + qJ(3);
t3 = sin(t5);
t4 = cos(t5);
t17 = mrSges(4,1) * t4 - mrSges(4,2) * t3;
t7 = cos(pkin(6));
t24 = -mrSges(2,1) - m(3) * pkin(1) - t7 * mrSges(3,1) + sin(pkin(6)) * mrSges(3,2) - t17;
t23 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(5,2) - mrSges(3,3) - mrSges(4,3);
t15 = t4 * mrSges(5,1) + t3 * mrSges(5,3);
t14 = pkin(3) * t4 + qJ(4) * t3;
t12 = m(5) * t14 + t15;
t8 = -pkin(5) - qJ(2);
t2 = pkin(2) * t7 + pkin(1);
t1 = [(-t26 * (t10 * t2 - t9 * t8) + t23 * t9 + (-t12 + t24) * t10) * g(2) + ((m(4) * t2 - m(5) * (-t14 - t2) + t15 - t24) * t9 + (t26 * t8 + t23) * t10) * g(1), (-g(1) * t9 + g(2) * t10) * (m(3) + t26), (-t12 - t17) * g(3) + ((-m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3)) * t4 + (m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1)) * t3) * t25, (g(3) * t4 - t25 * t3) * m(5)];
taug = t1(:);
