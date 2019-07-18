% Calculate Gravitation load on the joints for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(2,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:22
% EndTime: 2019-07-18 13:27:22
% DurationCPUTime: 0.11s
% Computational Cost: add. (68->20), mult. (54->21), div. (0->0), fcn. (32->6), ass. (0->16)
t17 = m(4) + m(5);
t8 = qJ(2) + qJ(3);
t16 = m(5) * pkin(2) + mrSges(4,1);
t7 = qJ(4) + t8;
t3 = sin(t7);
t4 = cos(t7);
t15 = -t4 * mrSges(5,1) + t3 * mrSges(5,2);
t14 = t17 * pkin(1) + mrSges(3,1);
t13 = -t3 * mrSges(5,1) - t4 * mrSges(5,2);
t5 = sin(t8);
t6 = cos(t8);
t12 = -t5 * mrSges(4,2) + t16 * t6 - t15;
t11 = t6 * mrSges(4,2) + t16 * t5 - t13;
t10 = cos(qJ(2));
t9 = sin(qJ(2));
t1 = [(m(2) + m(3) + t17) * g(2), (t10 * mrSges(3,2) + t14 * t9 + t11) * g(3) + (-t9 * mrSges(3,2) + t14 * t10 + t12) * g(1), t12 * g(1) + t11 * g(3), -g(1) * t15 - g(3) * t13];
taug  = t1(:);
