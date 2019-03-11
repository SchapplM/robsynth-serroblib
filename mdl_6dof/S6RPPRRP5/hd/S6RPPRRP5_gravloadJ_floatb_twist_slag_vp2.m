% Calculate Gravitation load on the joints for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:41
% EndTime: 2019-03-09 02:07:42
% DurationCPUTime: 0.50s
% Computational Cost: add. (168->61), mult. (335->71), div. (0->0), fcn. (284->6), ass. (0->31)
t59 = mrSges(6,1) + mrSges(7,1);
t52 = mrSges(6,2) + mrSges(7,2);
t16 = cos(qJ(5));
t58 = m(7) * (t16 * pkin(5) + pkin(4)) + m(6) * pkin(4);
t13 = sin(qJ(5));
t55 = -t52 * t13 + t59 * t16 + t58;
t54 = m(7) * (-qJ(6) - pkin(8)) - mrSges(7,3) - m(6) * pkin(8) - mrSges(6,3);
t15 = sin(qJ(1));
t18 = cos(qJ(1));
t51 = g(1) * t18 + g(2) * t15;
t45 = m(7) * pkin(5);
t50 = -m(5) - m(6) - m(7);
t49 = t45 + t59;
t32 = -m(4) + t50;
t47 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t14 = sin(qJ(4));
t17 = cos(qJ(4));
t23 = t14 * mrSges(5,1) + t17 * mrSges(5,2);
t46 = t58 * t14 + t54 * t17 + mrSges(2,1) - mrSges(3,2) + mrSges(4,3) + t23;
t41 = t18 * pkin(1) + t15 * qJ(2);
t40 = t15 * t13;
t39 = t15 * t16;
t36 = t18 * t13;
t35 = t18 * t16;
t33 = t18 * qJ(3) + t41;
t3 = -t14 * t36 - t39;
t1 = t14 * t40 - t35;
t10 = t18 * qJ(2);
t4 = t14 * t35 - t40;
t2 = -t14 * t39 - t36;
t5 = [(t40 * t45 - m(3) * t41 - m(4) * t33 - t59 * t4 - t52 * t3 + t50 * (-t15 * pkin(7) + t33) + t47 * t15 - t46 * t18) * g(2) + (t36 * t45 + t50 * (-t18 * pkin(7) + t10) - t59 * t2 + (-m(3) - m(4)) * t10 - t52 * t1 + t47 * t18 + (m(3) * pkin(1) + t32 * (-pkin(1) - qJ(3)) + t46) * t15) * g(1) (-g(1) * t15 + g(2) * t18) * (m(3) - t32) t51 * t32, g(3) * t23 + (t54 * g(3) + t51 * (-mrSges(5,1) - t55)) * t17 + (t55 * g(3) + t51 * (mrSges(5,2) + t54)) * t14 (t49 * t13 + t52 * t16) * g(3) * t17 + (t49 * t1 - t52 * t2) * g(2) + (-t49 * t3 + t52 * t4) * g(1) (-g(3) * t14 + t51 * t17) * m(7)];
taug  = t5(:);
