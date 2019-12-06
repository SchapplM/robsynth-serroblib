% Calculate Gravitation load on the joints for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:18
% EndTime: 2019-12-05 15:26:19
% DurationCPUTime: 0.32s
% Computational Cost: add. (121->38), mult. (168->48), div. (0->0), fcn. (131->8), ass. (0->21)
t30 = m(5) + m(6);
t22 = m(4) + t30;
t42 = mrSges(4,1) + mrSges(6,3) - mrSges(5,2);
t11 = sin(qJ(5));
t13 = cos(qJ(5));
t41 = -mrSges(6,1) * t11 - mrSges(6,2) * t13 + mrSges(4,2) - mrSges(5,3);
t10 = cos(pkin(7));
t9 = sin(pkin(7));
t34 = g(1) * t10 + g(2) * t9;
t8 = qJ(2) + pkin(8);
t6 = cos(t8);
t31 = g(3) * t6;
t27 = t11 * t9;
t26 = t13 * t9;
t14 = cos(qJ(2));
t7 = t14 * pkin(2);
t24 = t10 * t11;
t23 = t10 * t13;
t12 = sin(qJ(2));
t5 = sin(t8);
t1 = [(-m(2) - m(3) - t22) * g(3), (-m(4) * t7 - mrSges(3,1) * t14 + t12 * mrSges(3,2) - t30 * (pkin(3) * t6 + qJ(4) * t5 + t7) + (-m(6) * pkin(6) - t42) * t6 + t41 * t5) * g(3) + (mrSges(3,2) * t14 + (-m(6) * (-pkin(3) - pkin(6)) + m(5) * pkin(3) + t42) * t5 + (-qJ(4) * t30 + t41) * t6 + (pkin(2) * t22 + mrSges(3,1)) * t12) * t34, (-g(1) * t9 + g(2) * t10) * t22, (-t34 * t5 + t31) * t30, -g(1) * ((t23 * t5 - t27) * mrSges(6,1) + (-t24 * t5 - t26) * mrSges(6,2)) - g(2) * ((t26 * t5 + t24) * mrSges(6,1) + (-t27 * t5 + t23) * mrSges(6,2)) - (-mrSges(6,1) * t13 + mrSges(6,2) * t11) * t31];
taug = t1(:);
