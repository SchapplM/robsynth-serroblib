% Calculate Gravitation load on the joints for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:01:57
% EndTime: 2018-11-23 15:01:58
% DurationCPUTime: 0.70s
% Computational Cost: add. (906->90), mult. (1083->125), div. (0->0), fcn. (1029->14), ass. (0->59)
t88 = -m(5) - m(6) - m(7);
t71 = m(4) - t88;
t93 = mrSges(6,1) + mrSges(7,1);
t89 = -mrSges(6,2) - mrSges(7,2);
t47 = cos(qJ(5));
t92 = m(6) * pkin(4) + m(7) * (pkin(5) * t47 + pkin(4)) + mrSges(5,1);
t91 = m(6) * pkin(9) - m(7) * (-qJ(6) - pkin(9)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t82 = m(7) * pkin(5);
t87 = mrSges(3,1) - mrSges(4,2) + mrSges(5,3);
t44 = sin(qJ(5));
t86 = t89 * t44 + t93 * t47 + t92;
t85 = -t82 - t93;
t45 = sin(qJ(4));
t48 = cos(qJ(4));
t83 = -t71 * qJ(3) - t92 * t45 + t91 * t48 + mrSges(3,2) - mrSges(4,3);
t40 = sin(pkin(10));
t46 = sin(qJ(2));
t73 = pkin(6) - qJ(2);
t60 = cos(t73) / 0.2e1;
t72 = pkin(6) + qJ(2);
t63 = cos(t72);
t57 = t60 + t63 / 0.2e1;
t74 = cos(pkin(10));
t21 = t40 * t46 - t74 * t57;
t81 = t21 * t44;
t24 = t40 * t57 + t74 * t46;
t80 = t24 * t44;
t61 = sin(t72);
t58 = t61 / 0.2e1;
t62 = sin(t73);
t59 = t62 / 0.2e1;
t33 = t58 + t59;
t79 = t33 * t44;
t41 = sin(pkin(6));
t78 = t40 * t41;
t49 = cos(qJ(2));
t77 = t40 * t49;
t76 = t44 * t45;
t75 = t45 * t47;
t65 = t41 * t74;
t64 = t74 * t49;
t56 = t59 - t61 / 0.2e1;
t55 = t58 - t62 / 0.2e1;
t42 = cos(pkin(6));
t34 = t60 - t63 / 0.2e1;
t32 = t33 * pkin(2);
t28 = -t33 * t45 + t42 * t48;
t27 = t33 * t48 + t42 * t45;
t26 = t40 * t56 + t64;
t25 = -t40 * t55 + t64;
t23 = -t74 * t56 + t77;
t22 = t74 * t55 + t77;
t20 = t24 * pkin(2);
t19 = t21 * pkin(2);
t14 = -t21 * t45 + t48 * t65;
t13 = t21 * t48 + t45 * t65;
t12 = t24 * t45 + t48 * t78;
t11 = -t24 * t48 + t45 * t78;
t1 = [(-m(2) - m(3) - t71) * g(3) (-t79 * t82 - m(4) * t32 + t88 * (pkin(8) * t33 + t32) - t93 * (t34 * t75 + t79) + t89 * (t33 * t47 - t34 * t76) - t87 * t33 + t83 * t34) * g(3) + (t81 * t82 + m(4) * t19 + t88 * (-pkin(8) * t21 - t19) - t93 * (t23 * t75 - t81) + t89 * (-t21 * t47 - t23 * t76) + t87 * t21 + t83 * t23) * g(2) + (t80 * t82 + m(4) * t20 - t93 * (t26 * t75 - t80) + t89 * (-t24 * t47 - t26 * t76) + t88 * (-pkin(8) * t24 - t20) + t87 * t24 + t83 * t26) * g(1) (-g(1) * t24 - g(2) * t21 + g(3) * t33) * t71 (t86 * t27 - t28 * t91) * g(3) + (-t86 * t13 + t14 * t91) * g(2) + (t86 * t11 - t12 * t91) * g(1) (t89 * (-t28 * t47 - t34 * t44) + t85 * (-t28 * t44 + t34 * t47)) * g(3) + (t89 * (t14 * t47 - t22 * t44) + t85 * (t14 * t44 + t22 * t47)) * g(2) + (t89 * (-t12 * t47 - t25 * t44) + t85 * (-t12 * t44 + t25 * t47)) * g(1) (-g(1) * t11 + g(2) * t13 - g(3) * t27) * m(7)];
taug  = t1(:);
