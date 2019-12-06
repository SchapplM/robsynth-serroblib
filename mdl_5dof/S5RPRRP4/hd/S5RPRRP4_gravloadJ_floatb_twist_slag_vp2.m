% Calculate Gravitation load on the joints for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:05:30
% EndTime: 2019-12-05 18:05:33
% DurationCPUTime: 0.52s
% Computational Cost: add. (232->78), mult. (337->91), div. (0->0), fcn. (316->8), ass. (0->45)
t67 = mrSges(5,1) + mrSges(6,1);
t71 = mrSges(5,2) + mrSges(6,2);
t62 = m(3) + m(4) + m(5) + m(6);
t70 = t62 * qJ(2) - mrSges(2,2) + mrSges(3,3);
t59 = m(5) * pkin(3);
t29 = qJ(3) + qJ(4);
t24 = sin(t29);
t32 = sin(qJ(3));
t53 = pkin(3) * t32;
t18 = pkin(4) * t24 + t53;
t69 = m(6) * t18;
t68 = mrSges(4,1) + t59;
t25 = cos(t29);
t65 = mrSges(5,1) * t24 + t71 * t25;
t31 = cos(pkin(8));
t35 = cos(qJ(1));
t43 = t35 * t24;
t33 = sin(qJ(1));
t46 = t33 * t25;
t11 = t31 * t43 - t46;
t49 = t31 * t35;
t12 = -t24 * t33 - t25 * t49;
t64 = t67 * t11 - t71 * t12;
t10 = t31 * t46 - t43;
t50 = t31 * t33;
t9 = t24 * t50 + t25 * t35;
t63 = -t71 * t10 - t67 * t9;
t34 = cos(qJ(3));
t27 = t34 * pkin(3);
t19 = pkin(4) * t25 + t27;
t30 = sin(pkin(8));
t36 = -pkin(7) - pkin(6);
t61 = mrSges(2,1) + m(3) * pkin(1) + t31 * mrSges(3,1) - m(4) * (-pkin(2) * t31 - pkin(1)) - m(5) * (-(t27 + pkin(2)) * t31 - pkin(1)) - m(6) * (-(pkin(2) + t19) * t31 - pkin(1)) + (-mrSges(3,2) + m(4) * pkin(6) + mrSges(4,3) - m(5) * t36 + mrSges(5,3) - m(6) * (-qJ(5) + t36) + mrSges(6,3)) * t30;
t60 = m(5) * t53 + t69;
t58 = m(6) * pkin(4);
t52 = g(1) * t30;
t48 = t32 * t33;
t47 = t32 * t35;
t45 = t33 * t34;
t44 = t34 * t35;
t15 = t31 * t47 - t45;
t13 = t31 * t48 + t44;
t16 = -t31 * t44 - t48;
t14 = t31 * t45 - t47;
t1 = [(t14 * mrSges(4,1) - t13 * mrSges(4,2) + t67 * t10 + t61 * t33 - t47 * t59 - t71 * t9 + (-t69 - t70) * t35) * g(3) + (-t16 * mrSges(4,1) - t15 * mrSges(4,2) - t67 * t12 - t71 * t11 + t61 * t35 + (t60 + t70) * t33) * g(2), -(g(2) * t35 + g(3) * t33) * t62, (mrSges(4,1) * t32 + mrSges(6,1) * t24 + mrSges(4,2) * t34 + t60 + t65) * t52 + (-mrSges(4,2) * t16 - m(6) * (-t18 * t49 + t19 * t33) + t68 * t15 + t64) * g(3) + (-mrSges(4,2) * t14 - m(6) * (t18 * t50 + t19 * t35) - t68 * t13 + t63) * g(2), (-(-mrSges(6,1) - t58) * t24 + t65) * t52 + (t11 * t58 + t64) * g(3) + (-t58 * t9 + t63) * g(2), (g(1) * t31 + (g(2) * t33 - g(3) * t35) * t30) * m(6)];
taug = t1(:);
