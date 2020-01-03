% Calculate Gravitation load on the joints for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:25
% EndTime: 2019-12-31 17:33:25
% DurationCPUTime: 0.14s
% Computational Cost: add. (74->20), mult. (144->26), div. (0->0), fcn. (146->6), ass. (0->15)
t24 = m(5) + m(6);
t25 = m(6) * pkin(6) + pkin(3) * t24 + mrSges(4,1) - mrSges(5,2) + mrSges(6,3);
t14 = sin(pkin(7));
t15 = cos(pkin(7));
t16 = sin(qJ(3));
t17 = cos(qJ(3));
t3 = -t14 * t16 - t15 * t17;
t4 = -t14 * t17 + t15 * t16;
t23 = -g(1) * t4 + g(2) * t3;
t22 = m(3) + m(4) + t24;
t8 = sin(qJ(5));
t9 = cos(qJ(5));
t11 = t8 * mrSges(6,1) + t9 * mrSges(6,2);
t20 = t24 * qJ(4) - mrSges(4,2) + mrSges(5,3) + t11;
t1 = [(-m(2) - t22) * g(3), (-t14 * g(1) + t15 * g(2)) * t22, (t20 * t4 - t25 * t3) * g(2) + (t20 * t3 + t25 * t4) * g(1), t24 * t23, -g(3) * t11 + t23 * (mrSges(6,1) * t9 - mrSges(6,2) * t8)];
taug = t1(:);
