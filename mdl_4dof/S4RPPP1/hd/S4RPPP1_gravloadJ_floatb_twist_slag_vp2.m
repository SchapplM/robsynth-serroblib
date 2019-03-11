% Calculate Gravitation load on the joints for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:20
% EndTime: 2019-03-08 18:26:20
% DurationCPUTime: 0.15s
% Computational Cost: add. (88->37), mult. (203->48), div. (0->0), fcn. (199->6), ass. (0->22)
t32 = m(4) + m(5);
t33 = t32 * qJ(3) - mrSges(3,2) + mrSges(5,2) + mrSges(4,3);
t13 = sin(pkin(4));
t31 = g(3) * t13;
t16 = sin(qJ(1));
t17 = cos(qJ(1));
t26 = qJ(2) * t13;
t30 = t17 * pkin(1) + t16 * t26;
t12 = sin(pkin(6));
t29 = t12 * t16;
t15 = cos(pkin(4));
t28 = t15 * t17;
t14 = cos(pkin(6));
t27 = t16 * t14;
t22 = -t16 * pkin(1) + t17 * t26;
t20 = m(5) * qJ(4) + mrSges(3,1) - mrSges(4,2) + mrSges(5,3);
t18 = mrSges(2,2) + (-m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1) - mrSges(3,3)) * t13;
t6 = t14 * t17 - t15 * t29;
t5 = t12 * t17 + t15 * t27;
t4 = t12 * t28 + t27;
t3 = -t14 * t28 + t29;
t1 = [(-mrSges(2,1) * t17 - m(3) * t30 - t20 * t6 - t33 * t5 + t18 * t16 - t32 * (t6 * pkin(2) + t30)) * g(2) + (t16 * mrSges(2,1) - m(3) * t22 + t20 * t4 + t33 * t3 + t18 * t17 + t32 * (t4 * pkin(2) - t22)) * g(1) (-g(3) * t15 + (-g(1) * t16 + g(2) * t17) * t13) * (m(3) + t32) t32 * (-g(1) * t5 - g(2) * t3 + t14 * t31) (-g(1) * t6 - g(2) * t4 - t12 * t31) * m(5)];
taug  = t1(:);
