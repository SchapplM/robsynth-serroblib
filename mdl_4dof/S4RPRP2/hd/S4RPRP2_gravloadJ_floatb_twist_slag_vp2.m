% Calculate Gravitation load on the joints for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:42
% EndTime: 2019-03-08 18:30:42
% DurationCPUTime: 0.16s
% Computational Cost: add. (67->31), mult. (119->34), div. (0->0), fcn. (108->4), ass. (0->17)
t10 = sin(qJ(3));
t11 = sin(qJ(1));
t12 = cos(qJ(1));
t18 = cos(qJ(3));
t1 = -t11 * t10 - t12 * t18;
t15 = mrSges(4,2) + mrSges(5,2);
t23 = t15 * t1;
t2 = t12 * t10 - t11 * t18;
t22 = t15 * t2;
t16 = mrSges(4,1) + mrSges(5,1);
t14 = m(3) + m(4) + m(5);
t20 = m(5) * pkin(3);
t21 = t20 + t16;
t17 = mrSges(2,1) + mrSges(3,1);
t13 = -t10 * t20 + mrSges(2,2) - mrSges(3,3);
t6 = t18 * pkin(3) + pkin(2);
t3 = [(-t16 * t2 + t23 + (m(3) * pkin(1) - m(4) * (-pkin(1) - pkin(2)) - m(5) * (-pkin(1) - t6) + t17) * t11 + (-t14 * qJ(2) + t13) * t12) * g(1) + (t22 + (-m(4) * pkin(2) - m(5) * t6 - t17) * t12 + t13 * t11 + t16 * t1 - t14 * (t12 * pkin(1) + t11 * qJ(2))) * g(2) (-g(1) * t11 + g(2) * t12) * t14 (-t21 * t1 - t22) * g(2) + (t21 * t2 - t23) * g(1), g(3) * m(5)];
taug  = t3(:);
