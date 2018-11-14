% Calculate Gravitation load on the joints for
% S4RPRP1
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4RPRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:27
% EndTime: 2018-11-14 13:48:28
% DurationCPUTime: 0.13s
% Computational Cost: add. (103->31), mult. (69->32), div. (0->0), fcn. (44->6), ass. (0->17)
t25 = -mrSges(4,1) - mrSges(5,1);
t24 = mrSges(4,2) - mrSges(5,3);
t23 = m(4) + m(5);
t15 = qJ(1) + pkin(6);
t13 = qJ(3) + t15;
t10 = cos(t13);
t9 = sin(t13);
t22 = t10 * pkin(3) + t9 * qJ(4);
t12 = cos(t15);
t17 = cos(qJ(1));
t21 = t17 * pkin(1) + pkin(2) * t12;
t20 = m(3) + t23;
t19 = t25 * t10 + t24 * t9;
t18 = (m(5) * pkin(3) - t25) * t9 + (-m(5) * qJ(4) + t24) * t10;
t16 = sin(qJ(1));
t11 = sin(t15);
t1 = [(t16 * mrSges(2,2) - t12 * mrSges(3,1) + t11 * mrSges(3,2) - m(4) * t21 - m(5) * (t21 + t22) + (-m(3) * pkin(1) - mrSges(2,1)) * t17 + t19) * g(2) + (mrSges(2,2) * t17 + mrSges(3,2) * t12 + (t23 * pkin(2) + mrSges(3,1)) * t11 + (t20 * pkin(1) + mrSges(2,1)) * t16 + t18) * g(1), -t20 * g(3) (-m(5) * t22 + t19) * g(2) + t18 * g(1) (-g(1) * t9 + g(2) * t10) * m(5)];
taug  = t1(:);
