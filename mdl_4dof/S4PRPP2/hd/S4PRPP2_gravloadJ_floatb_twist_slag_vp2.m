% Calculate Gravitation load on the joints for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
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
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP2_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:18:49
% EndTime: 2019-03-08 18:18:49
% DurationCPUTime: 0.08s
% Computational Cost: add. (38->18), mult. (40->18), div. (0->0), fcn. (22->4), ass. (0->10)
t8 = m(4) + m(5);
t9 = t8 * pkin(2) + mrSges(3,1);
t7 = m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1);
t6 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3);
t5 = cos(qJ(2));
t4 = sin(qJ(2));
t3 = qJ(2) + pkin(5);
t2 = cos(t3);
t1 = sin(t3);
t10 = [(-m(2) - m(3) - t8) * g(2) (t4 * mrSges(3,2) + t6 * t1 - t7 * t2 - t9 * t5) * g(2) + (t5 * mrSges(3,2) + t7 * t1 + t6 * t2 + t9 * t4) * g(1), -t8 * g(3) (-g(1) * t1 + g(2) * t2) * m(5)];
taug  = t10(:);
