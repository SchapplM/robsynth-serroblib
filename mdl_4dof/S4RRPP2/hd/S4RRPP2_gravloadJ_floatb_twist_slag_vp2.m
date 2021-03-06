% Calculate Gravitation load on the joints for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:59
% EndTime: 2019-03-08 18:34:00
% DurationCPUTime: 0.14s
% Computational Cost: add. (103->32), mult. (87->34), div. (0->0), fcn. (58->4), ass. (0->17)
t32 = -mrSges(3,1) - mrSges(4,1) - mrSges(5,1);
t31 = -mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t15 = qJ(1) + qJ(2);
t12 = sin(t15);
t13 = cos(t15);
t30 = t31 * t13 + (-m(5) * (-pkin(2) - pkin(3)) - t32) * t12;
t29 = t31 * t12 + t32 * t13;
t16 = sin(qJ(1));
t26 = pkin(1) * t16;
t17 = cos(qJ(1));
t14 = t17 * pkin(1);
t25 = t13 * pkin(2) + t12 * qJ(3);
t24 = t14 + t25;
t2 = t13 * qJ(3);
t23 = -pkin(2) * t12 + t2;
t10 = t13 * pkin(3);
t1 = [(-mrSges(2,1) * t17 + t16 * mrSges(2,2) - m(3) * t14 - m(4) * t24 - m(5) * (t10 + t24) + t29) * g(2) + (t16 * mrSges(2,1) + mrSges(2,2) * t17 + m(3) * t26 - m(4) * (t23 - t26) - m(5) * (t2 - t26) + t30) * g(1) (-m(4) * t25 - m(5) * (t10 + t25) + t29) * g(2) + (-m(4) * t23 - m(5) * t2 + t30) * g(1) (m(4) + m(5)) * (-g(1) * t12 + g(2) * t13) g(3) * m(5)];
taug  = t1(:);
