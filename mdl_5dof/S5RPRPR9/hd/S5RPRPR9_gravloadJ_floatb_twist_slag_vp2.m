% Calculate Gravitation load on the joints for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:43
% EndTime: 2019-12-31 18:23:44
% DurationCPUTime: 0.41s
% Computational Cost: add. (199->67), mult. (232->77), div. (0->0), fcn. (194->8), ass. (0->36)
t58 = -mrSges(4,1) + mrSges(5,2);
t57 = mrSges(4,2) - mrSges(5,3);
t17 = qJ(1) + pkin(8);
t12 = sin(t17);
t13 = cos(t17);
t50 = g(1) * t13 + g(2) * t12;
t46 = m(5) + m(6);
t54 = -m(4) - t46;
t19 = sin(qJ(3));
t22 = cos(qJ(3));
t51 = t57 * t19 + t58 * t22;
t49 = -t22 * mrSges(6,3) + t51;
t47 = -m(6) * pkin(4) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t43 = g(3) * t22;
t20 = sin(qJ(1));
t42 = t20 * pkin(1);
t15 = t22 * pkin(3);
t23 = cos(qJ(1));
t16 = t23 * pkin(1);
t41 = t13 * t22;
t18 = sin(qJ(5));
t40 = t19 * t18;
t21 = cos(qJ(5));
t39 = t19 * t21;
t14 = t19 * qJ(4);
t37 = t15 + t14;
t35 = t13 * pkin(2) + t12 * pkin(6) + t16;
t33 = -pkin(2) - t14;
t32 = pkin(3) * t41 + t13 * t14 + t35;
t31 = m(6) * (-pkin(3) - pkin(7)) - mrSges(6,3);
t27 = t18 * mrSges(6,1) + t21 * mrSges(6,2);
t4 = -t12 * t40 + t13 * t21;
t3 = t12 * t39 + t13 * t18;
t2 = t12 * t21 + t13 * t40;
t1 = -t12 * t18 + t13 * t39;
t5 = [(-t23 * mrSges(2,1) + t20 * mrSges(2,2) - m(3) * t16 - m(4) * t35 - m(5) * t32 - m(6) * (pkin(7) * t41 + t32) - t2 * mrSges(6,1) - t1 * mrSges(6,2) + t47 * t12 + (-mrSges(3,1) + t49) * t13) * g(2) + (m(3) * t42 + t20 * mrSges(2,1) - t4 * mrSges(6,1) + t23 * mrSges(2,2) + t3 * mrSges(6,2) + t54 * (t13 * pkin(6) - t42) + t47 * t13 + (mrSges(3,1) + m(4) * pkin(2) - m(5) * (t33 - t15) - m(6) * t33 - t31 * t22 - t51) * t12) * g(1), (-m(3) + t54) * g(3), (-m(5) * t37 - m(6) * (t22 * pkin(7) + t37) - t27 * t19 + t49) * g(3) + ((m(5) * pkin(3) - t31 - t58) * t19 + (-qJ(4) * t46 - t27 + t57) * t22) * t50, (-t19 * t50 + t43) * t46, -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * (mrSges(6,1) * t3 + t4 * mrSges(6,2)) - (-mrSges(6,1) * t21 + mrSges(6,2) * t18) * t43];
taug = t5(:);
