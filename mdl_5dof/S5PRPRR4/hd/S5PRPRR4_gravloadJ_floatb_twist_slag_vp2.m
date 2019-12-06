% Calculate Gravitation load on the joints for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:40
% EndTime: 2019-12-05 15:49:44
% DurationCPUTime: 0.69s
% Computational Cost: add. (331->69), mult. (850->113), div. (0->0), fcn. (1035->12), ass. (0->41)
t68 = m(5) + m(6);
t33 = sin(qJ(5));
t36 = cos(qJ(5));
t65 = -m(6) * pkin(4) - t36 * mrSges(6,1) + t33 * mrSges(6,2) - mrSges(5,1);
t64 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t29 = sin(pkin(9));
t31 = cos(pkin(9));
t35 = sin(qJ(2));
t32 = cos(pkin(5));
t38 = cos(qJ(2));
t57 = t32 * t38;
t67 = -t29 * t35 + t31 * t57;
t28 = sin(pkin(10));
t54 = cos(pkin(10));
t45 = -t35 * t28 + t38 * t54;
t34 = sin(qJ(4));
t37 = cos(qJ(4));
t63 = t64 * t34 + t65 * t37 - mrSges(4,1);
t62 = t33 * mrSges(6,1) + t36 * mrSges(6,2) + t68 * pkin(7) - mrSges(4,2) + mrSges(5,3);
t30 = sin(pkin(5));
t60 = t30 * t34;
t59 = t30 * t37;
t58 = t32 * t35;
t53 = m(4) + t68;
t26 = t30 * t38 * pkin(2);
t49 = t67 * pkin(2);
t22 = -t38 * t28 - t35 * t54;
t20 = t22 * t32;
t9 = -t31 * t20 + t29 * t45;
t10 = -t29 * t20 - t31 * t45;
t46 = -t29 * t57 - t31 * t35;
t42 = t46 * pkin(2);
t41 = t45 * t32;
t19 = t22 * t30;
t18 = t45 * t30;
t14 = -t19 * t37 + t32 * t34;
t11 = t31 * t22 - t29 * t41;
t8 = t29 * t22 + t31 * t41;
t4 = -t10 * t37 + t29 * t60;
t2 = -t31 * t60 + t9 * t37;
t1 = [(-m(2) - m(3) - t53) * g(3), (-(mrSges(3,1) * t38 - mrSges(3,2) * t35) * t30 - m(4) * t26 - t68 * (t18 * pkin(3) + t26) + t62 * t19 + t63 * t18) * g(3) + (-t67 * mrSges(3,1) - (-t29 * t38 - t31 * t58) * mrSges(3,2) - m(4) * t49 - t68 * (t8 * pkin(3) + t49) + t63 * t8 - t62 * t9) * g(2) + (-t46 * mrSges(3,1) - (t29 * t58 - t31 * t38) * mrSges(3,2) - m(4) * t42 - t68 * (t11 * pkin(3) + t42) + t63 * t11 + t62 * t10) * g(1), (-t32 * g(3) + (-g(1) * t29 + g(2) * t31) * t30) * t53, (t64 * t14 + t65 * (t19 * t34 + t32 * t37)) * g(3) + (t64 * t2 + t65 * (-t31 * t59 - t9 * t34)) * g(2) + (t64 * t4 + t65 * (t10 * t34 + t29 * t59)) * g(1), -g(1) * ((-t11 * t36 - t4 * t33) * mrSges(6,1) + (t11 * t33 - t4 * t36) * mrSges(6,2)) - g(2) * ((-t2 * t33 - t8 * t36) * mrSges(6,1) + (-t2 * t36 + t8 * t33) * mrSges(6,2)) - g(3) * ((-t14 * t33 - t18 * t36) * mrSges(6,1) + (-t14 * t36 + t18 * t33) * mrSges(6,2))];
taug = t1(:);
