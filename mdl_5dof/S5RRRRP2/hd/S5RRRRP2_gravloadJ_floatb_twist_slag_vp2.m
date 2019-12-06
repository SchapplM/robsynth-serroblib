% Calculate Gravitation load on the joints for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:47:33
% EndTime: 2019-12-05 18:47:34
% DurationCPUTime: 0.25s
% Computational Cost: add. (289->56), mult. (244->54), div. (0->0), fcn. (186->8), ass. (0->33)
t30 = qJ(3) + qJ(4);
t24 = sin(t30);
t47 = mrSges(5,2) + mrSges(6,2);
t62 = t47 * t24;
t26 = cos(t30);
t43 = t47 * t26;
t48 = mrSges(5,1) + mrSges(6,1);
t61 = t48 * t26;
t34 = cos(qJ(3));
t28 = t34 * pkin(3);
t22 = pkin(4) * t26;
t46 = t22 + t28;
t32 = sin(qJ(3));
t49 = t32 * mrSges(4,2);
t60 = -t49 + m(4) * pkin(2) + m(5) * (t28 + pkin(2)) + m(6) * (pkin(2) + t46) + t34 * mrSges(4,1) + t61 + mrSges(3,1);
t36 = -pkin(8) - pkin(7);
t59 = -m(4) * pkin(7) - mrSges(6,3) - mrSges(5,3) - mrSges(4,3) + mrSges(3,2) + m(6) * (-qJ(5) + t36) + m(5) * t36;
t44 = m(5) * pkin(3) + mrSges(4,1);
t55 = -t44 * t32 + m(6) * (-t32 * pkin(3) - pkin(4) * t24) - mrSges(4,2) * t34;
t54 = m(6) * pkin(4);
t31 = qJ(1) + qJ(2);
t27 = cos(t31);
t53 = g(3) * t27;
t25 = sin(t31);
t52 = t24 * t25;
t45 = -t25 * t43 - t48 * t52;
t41 = -t61 + t62;
t40 = mrSges(2,1) + (m(3) + m(4) + m(5) + m(6)) * pkin(1);
t38 = -t59 * t25 + (t60 - t62) * t27;
t37 = t60 * t25 + t59 * t27 - t47 * t52;
t35 = cos(qJ(1));
t33 = sin(qJ(1));
t1 = [(t35 * mrSges(2,2) + t40 * t33 + t37) * g(3) + (-t33 * mrSges(2,2) + t40 * t35 + t38) * g(2), t38 * g(2) + t37 * g(3), (t55 * t25 + t45) * g(2) + (-m(6) * t46 - t44 * t34 + t41 + t49) * g(1) + (t48 * t24 + t43 - t55) * t53, (-t52 * t54 + t45) * g(2) + (-m(6) * t22 + t41) * g(1) + (t43 + (t48 + t54) * t24) * t53, (-g(2) * t27 - g(3) * t25) * m(6)];
taug = t1(:);
