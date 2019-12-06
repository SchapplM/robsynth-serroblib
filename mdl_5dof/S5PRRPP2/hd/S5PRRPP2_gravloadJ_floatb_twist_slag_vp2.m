% Calculate Gravitation load on the joints for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:43
% EndTime: 2019-12-05 16:08:46
% DurationCPUTime: 0.52s
% Computational Cost: add. (183->54), mult. (299->69), div. (0->0), fcn. (267->8), ass. (0->28)
t56 = mrSges(5,1) + mrSges(6,1);
t55 = -mrSges(5,2) + mrSges(6,3);
t44 = m(5) + m(6);
t54 = t44 * (-qJ(4) - pkin(6)) - m(4) * pkin(6) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t48 = pkin(3) * t44 + mrSges(4,1);
t12 = qJ(3) + pkin(8);
t10 = sin(t12);
t11 = cos(t12);
t16 = sin(qJ(3));
t18 = cos(qJ(3));
t53 = m(4) * pkin(2) + t18 * mrSges(4,1) - t16 * mrSges(4,2) + t55 * t10 + t56 * t11 + mrSges(3,1);
t13 = sin(pkin(7));
t14 = cos(pkin(7));
t49 = g(1) * t14 + g(2) * t13;
t46 = m(6) * pkin(4) + t56;
t45 = -m(6) * qJ(5) - t55;
t17 = sin(qJ(2));
t40 = g(3) * t17;
t19 = cos(qJ(2));
t38 = t16 * t19;
t37 = t18 * t19;
t36 = t19 * t10;
t35 = t19 * t11;
t27 = pkin(4) * t11 + qJ(5) * t10;
t9 = t18 * pkin(3) + pkin(2);
t3 = -t13 * t11 + t14 * t36;
t1 = t14 * t11 + t13 * t36;
t2 = [(-m(2) - m(3) - m(4) - t44) * g(3), ((-m(6) * t27 - t44 * t9 - t53) * g(3) + t49 * t54) * t19 + (t54 * g(3) + t49 * (m(5) * t9 - m(6) * (-t27 - t9) + t53)) * t17, (mrSges(4,2) * t18 + t46 * t10 + t45 * t11 + t48 * t16) * t40 + (-(-t13 * t37 + t14 * t16) * mrSges(4,2) + t45 * (-t14 * t10 + t13 * t35) + t46 * t1 - t48 * (-t13 * t38 - t14 * t18)) * g(2) + (-(-t13 * t16 - t14 * t37) * mrSges(4,2) + t45 * (t13 * t10 + t14 * t35) + t46 * t3 + t48 * (-t13 * t18 + t14 * t38)) * g(1), (t19 * g(3) - t17 * t49) * t44, (-g(1) * t3 - g(2) * t1 - t10 * t40) * m(6)];
taug = t2(:);
