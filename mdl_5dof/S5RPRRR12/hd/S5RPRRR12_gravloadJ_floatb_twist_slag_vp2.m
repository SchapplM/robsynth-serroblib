% Calculate Gravitation load on the joints for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:29
% EndTime: 2019-12-31 19:12:30
% DurationCPUTime: 0.39s
% Computational Cost: add. (193->75), mult. (269->94), div. (0->0), fcn. (224->8), ass. (0->47)
t23 = qJ(3) + qJ(4);
t18 = sin(t23);
t19 = cos(t23);
t74 = t18 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t19;
t73 = -m(3) - m(4);
t72 = -m(5) - m(6);
t41 = -t18 * pkin(4) + t19 * pkin(8);
t71 = m(6) * t41;
t29 = cos(qJ(1));
t24 = sin(qJ(5));
t56 = mrSges(6,2) * t24;
t45 = t19 * t56;
t57 = mrSges(5,2) * t18;
t70 = (-t45 - t57) * t29;
t26 = sin(qJ(1));
t27 = cos(qJ(5));
t51 = t27 * mrSges(6,1);
t44 = t19 * t51;
t54 = t26 * t19;
t55 = t18 * t26;
t69 = -mrSges(5,1) * t54 - mrSges(6,3) * t55 - (t44 - t45) * t26;
t68 = -(-t51 + t56) * t18 + t74;
t67 = mrSges(5,1) * t19 + t18 * mrSges(6,3) + t44;
t66 = -g(1) * t26 + g(2) * t29;
t65 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t25 = sin(qJ(3));
t28 = cos(qJ(3));
t38 = t25 * mrSges(4,1) + t28 * mrSges(4,2);
t64 = mrSges(2,2) - mrSges(3,3) - t38 + t71 - t74;
t63 = pkin(3) * t28;
t59 = t25 * pkin(3);
t53 = t26 * t24;
t52 = t26 * t27;
t50 = t29 * t24;
t49 = t29 * t27;
t48 = pkin(4) * t54 + pkin(8) * t55;
t47 = t29 * pkin(1) + t26 * qJ(2);
t46 = m(5) * t63;
t21 = t29 * qJ(2);
t42 = -t26 * pkin(1) + t21;
t40 = -pkin(4) * t19 - pkin(8) * t18;
t30 = -pkin(7) - pkin(6);
t4 = t18 * t49 - t53;
t3 = t18 * t50 + t52;
t2 = t18 * t52 + t50;
t1 = -t18 * t53 + t49;
t5 = [(-t2 * mrSges(6,1) - t1 * mrSges(6,2) + t73 * t47 + t72 * (t26 * t59 - t29 * t30 + t47) + (-m(4) * pkin(6) - t65) * t29 + t64 * t26) * g(2) + (-m(3) * t42 - m(4) * t21 - t4 * mrSges(6,1) + t3 * mrSges(6,2) + t72 * (t26 * t30 + t29 * t59 + t42) + (-m(4) * (-pkin(1) - pkin(6)) + t65) * t26 + t64 * t29) * g(1), t66 * (-t72 - t73), t66 * (mrSges(4,1) * t28 - mrSges(4,2) * t25) + ((t46 - m(6) * (t40 - t63) + t67) * t29 + t70) * g(2) + (-(t46 - t57) * t26 - m(6) * (t26 * t63 + t48) + t69) * g(1) + (t38 + m(5) * t59 - m(6) * (t41 - t59) + t68) * g(3), (t68 - t71) * g(3) + ((-m(6) * t40 + t67) * t29 + t70) * g(2) + (-m(6) * t48 + mrSges(5,2) * t55 + t69) * g(1), -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * (mrSges(6,1) * t3 + mrSges(6,2) * t4) - g(3) * (-mrSges(6,1) * t24 - mrSges(6,2) * t27) * t19];
taug = t5(:);
