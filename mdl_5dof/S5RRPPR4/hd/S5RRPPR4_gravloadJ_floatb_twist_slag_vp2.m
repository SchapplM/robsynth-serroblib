% Calculate Gravitation load on the joints for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:40
% EndTime: 2019-12-31 19:27:41
% DurationCPUTime: 0.18s
% Computational Cost: add. (247->48), mult. (219->51), div. (0->0), fcn. (208->8), ass. (0->27)
t50 = -mrSges(3,1) - mrSges(4,1);
t49 = mrSges(3,2) - mrSges(4,3);
t48 = mrSges(5,2) - mrSges(6,3);
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t47 = -t30 * mrSges(6,1) + t28 * mrSges(6,2);
t46 = mrSges(5,1) - t47;
t45 = m(5) + m(6);
t27 = qJ(1) + qJ(2);
t24 = sin(t27);
t25 = cos(t27);
t41 = t25 * pkin(2) + t24 * qJ(3);
t40 = cos(pkin(8));
t39 = sin(pkin(8));
t38 = m(4) + t45;
t22 = t25 * pkin(3);
t37 = t22 + t41;
t31 = cos(qJ(1));
t26 = t31 * pkin(1);
t36 = t26 + t41;
t11 = -t24 * t39 - t25 * t40;
t12 = -t24 * t40 + t25 * t39;
t34 = -t11 * pkin(4) + pkin(7) * t12 + t37;
t33 = t46 * t11 + t48 * t12 + t49 * t24 + t50 * t25;
t32 = (m(4) * pkin(2) + t45 * (pkin(2) + pkin(3)) - t50) * t24 + (-t38 * qJ(3) + t49) * t25 + (-m(6) * pkin(7) + t48) * t11 + (-m(6) * pkin(4) - t46) * t12;
t29 = sin(qJ(1));
t1 = [(t29 * mrSges(2,2) - m(4) * t36 - m(5) * (t22 + t36) - m(6) * (t26 + t34) + (-m(3) * pkin(1) - mrSges(2,1)) * t31 + t33) * g(2) + (mrSges(2,2) * t31 + (mrSges(2,1) + (m(3) + t38) * pkin(1)) * t29 + t32) * g(1), (-m(4) * t41 - m(5) * t37 - m(6) * t34 + t33) * g(2) + t32 * g(1), (-g(1) * t24 + g(2) * t25) * t38, t45 * g(3), -g(3) * t47 + (-g(1) * t11 - g(2) * t12) * (mrSges(6,1) * t28 + mrSges(6,2) * t30)];
taug = t1(:);
