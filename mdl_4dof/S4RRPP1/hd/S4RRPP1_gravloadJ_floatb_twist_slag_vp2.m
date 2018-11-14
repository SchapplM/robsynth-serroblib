% Calculate Gravitation load on the joints for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4RRPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:29
% EndTime: 2018-11-14 13:51:29
% DurationCPUTime: 0.14s
% Computational Cost: add. (118->33), mult. (80->33), div. (0->0), fcn. (52->6), ass. (0->17)
t24 = -mrSges(4,1) - mrSges(5,1);
t23 = mrSges(4,2) - mrSges(5,3);
t22 = m(4) + m(5);
t16 = qJ(1) + qJ(2);
t14 = cos(t16);
t11 = pkin(2) * t14;
t12 = pkin(6) + t16;
t10 = cos(t12);
t9 = sin(t12);
t21 = t10 * pkin(3) + t9 * qJ(4) + t11;
t13 = sin(t16);
t20 = -t14 * mrSges(3,1) + t13 * mrSges(3,2) + t24 * t10 + t23 * t9;
t19 = t14 * mrSges(3,2) + (m(5) * pkin(3) - t24) * t9 + (t22 * pkin(2) + mrSges(3,1)) * t13 + (-m(5) * qJ(4) + t23) * t10;
t18 = cos(qJ(1));
t17 = sin(qJ(1));
t15 = t18 * pkin(1);
t1 = [(t17 * mrSges(2,2) - m(4) * (t11 + t15) - m(5) * (t15 + t21) + (-m(3) * pkin(1) - mrSges(2,1)) * t18 + t20) * g(2) + (mrSges(2,2) * t18 + (mrSges(2,1) + (m(3) + t22) * pkin(1)) * t17 + t19) * g(1) (-m(4) * t11 - m(5) * t21 + t20) * g(2) + t19 * g(1), -t22 * g(3) (-g(1) * t9 + g(2) * t10) * m(5)];
taug  = t1(:);
