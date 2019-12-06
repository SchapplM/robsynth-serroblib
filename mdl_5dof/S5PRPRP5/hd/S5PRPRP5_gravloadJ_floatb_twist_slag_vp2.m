% Calculate Gravitation load on the joints for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:34
% EndTime: 2019-12-05 15:37:36
% DurationCPUTime: 0.35s
% Computational Cost: add. (166->43), mult. (257->52), div. (0->0), fcn. (224->8), ass. (0->25)
t42 = m(5) + m(6);
t45 = t42 * (-pkin(6) - qJ(3)) - m(4) * qJ(3) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t44 = mrSges(5,1) + mrSges(6,1);
t43 = -mrSges(5,2) + mrSges(6,3);
t11 = pkin(8) + qJ(4);
t10 = cos(t11);
t14 = cos(pkin(8));
t9 = sin(t11);
t41 = m(4) * pkin(2) + t14 * mrSges(4,1) - sin(pkin(8)) * mrSges(4,2) + mrSges(3,1) + t43 * t9 + t44 * t10;
t13 = sin(pkin(7));
t15 = cos(pkin(7));
t39 = g(1) * t15 + g(2) * t13;
t36 = m(6) * pkin(4) + t44;
t35 = -m(6) * qJ(5) - t43;
t17 = sin(qJ(2));
t32 = g(3) * t17;
t18 = cos(qJ(2));
t31 = t18 * t9;
t30 = t18 * t10;
t29 = m(4) + t42;
t24 = pkin(4) * t10 + qJ(5) * t9;
t8 = pkin(3) * t14 + pkin(2);
t3 = -t13 * t10 + t15 * t31;
t1 = t15 * t10 + t13 * t31;
t2 = [(-m(2) - m(3) - t29) * g(3), ((-m(6) * t24 - t42 * t8 - t41) * g(3) + t39 * t45) * t18 + (t45 * g(3) + t39 * (m(5) * t8 - m(6) * (-t24 - t8) + t41)) * t17, (g(3) * t18 - t39 * t17) * t29, (t35 * t10 + t36 * t9) * t32 + (t35 * (t13 * t30 - t15 * t9) + t36 * t1) * g(2) + (t35 * (t13 * t9 + t15 * t30) + t36 * t3) * g(1), (-g(1) * t3 - g(2) * t1 - t9 * t32) * m(6)];
taug = t2(:);
