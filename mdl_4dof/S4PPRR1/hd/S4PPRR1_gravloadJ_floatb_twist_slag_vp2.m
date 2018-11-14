% Calculate Gravitation load on the joints for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4PPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:18
% EndTime: 2018-11-14 13:40:18
% DurationCPUTime: 0.10s
% Computational Cost: add. (50->17), mult. (68->25), div. (0->0), fcn. (62->6), ass. (0->16)
t21 = -m(5) * pkin(3) - mrSges(4,1);
t11 = qJ(3) + qJ(4);
t10 = cos(t11);
t12 = sin(pkin(6));
t13 = cos(pkin(6));
t9 = sin(t11);
t5 = -t10 * t13 - t12 * t9;
t6 = -t10 * t12 + t13 * t9;
t19 = t5 * mrSges(5,1) + t6 * mrSges(5,2);
t18 = -t6 * mrSges(5,1) + t5 * mrSges(5,2);
t17 = m(3) + m(4) + m(5);
t14 = sin(qJ(3));
t15 = cos(qJ(3));
t16 = t12 * t15 - t13 * t14;
t7 = -t12 * t14 - t13 * t15;
t1 = [(-m(2) - t17) * g(3) (-g(1) * t12 + g(2) * t13) * t17 (t16 * mrSges(4,2) + t21 * t7 - t19) * g(2) + (-t7 * mrSges(4,2) + t21 * t16 - t18) * g(1), -g(1) * t18 - g(2) * t19];
taug  = t1(:);
