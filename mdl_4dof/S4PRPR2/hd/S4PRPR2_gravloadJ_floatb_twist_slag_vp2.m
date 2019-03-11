% Calculate Gravitation load on the joints for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:48
% EndTime: 2019-03-08 18:21:48
% DurationCPUTime: 0.09s
% Computational Cost: add. (46->20), mult. (40->20), div. (0->0), fcn. (22->6), ass. (0->14)
t13 = m(4) + m(5);
t14 = t13 * pkin(2) + mrSges(3,1);
t7 = qJ(2) + pkin(6);
t12 = m(5) * pkin(3) + mrSges(4,1);
t6 = qJ(4) + t7;
t2 = sin(t6);
t3 = cos(t6);
t11 = t3 * mrSges(5,1) - t2 * mrSges(5,2);
t10 = -t2 * mrSges(5,1) - t3 * mrSges(5,2);
t9 = cos(qJ(2));
t8 = sin(qJ(2));
t5 = cos(t7);
t4 = sin(t7);
t1 = [(-m(2) - m(3) - t13) * g(2) (t8 * mrSges(3,2) + t4 * mrSges(4,2) - t12 * t5 - t14 * t9 - t11) * g(2) + (t9 * mrSges(3,2) + t5 * mrSges(4,2) + t12 * t4 + t14 * t8 - t10) * g(1), -t13 * g(3), -g(1) * t10 - g(2) * t11];
taug  = t1(:);
