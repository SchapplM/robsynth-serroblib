% Calculate Gravitation load on the joints for
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
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
% Datum: 2019-03-08 18:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:54
% EndTime: 2019-03-08 18:17:55
% DurationCPUTime: 0.10s
% Computational Cost: add. (56->18), mult. (48->20), div. (0->0), fcn. (28->2), ass. (0->7)
t11 = m(4) + m(5);
t10 = mrSges(3,1) + mrSges(5,3) - mrSges(4,2);
t9 = mrSges(3,2) - mrSges(4,3) - mrSges(5,2);
t6 = pkin(5) + qJ(2);
t5 = cos(t6);
t4 = sin(t6);
t1 = [(-m(2) - m(3) - t11) * g(3) (-t11 * (t5 * pkin(2) + t4 * qJ(3)) + (-m(5) * qJ(4) - t10) * t5 + t9 * t4) * g(2) + ((m(4) * pkin(2) - m(5) * (-pkin(2) - qJ(4)) + t10) * t4 + (-t11 * qJ(3) + t9) * t5) * g(1), t11 * (-g(1) * t4 + g(2) * t5) (-g(1) * t5 - g(2) * t4) * m(5)];
taug  = t1(:);
