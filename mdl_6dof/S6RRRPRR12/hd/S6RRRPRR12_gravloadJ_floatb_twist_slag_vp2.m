% Calculate Gravitation load on the joints for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 18:00
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:59:52
% EndTime: 2018-11-23 17:59:54
% DurationCPUTime: 1.40s
% Computational Cost: add. (1685->127), mult. (1879->158), div. (0->0), fcn. (1881->18), ass. (0->68)
t53 = pkin(12) + qJ(5);
t49 = qJ(6) + t53;
t44 = sin(t49);
t45 = cos(t49);
t47 = sin(t53);
t48 = cos(t53);
t55 = cos(pkin(12));
t46 = t55 * pkin(4) + pkin(3);
t54 = sin(pkin(12));
t65 = -mrSges(4,1) - m(5) * pkin(3) - t55 * mrSges(5,1) + t54 * mrSges(5,2) - m(7) * (pkin(5) * t48 + t46) - m(6) * t46;
t143 = -t48 * mrSges(6,1) - t45 * mrSges(7,1) + t47 * mrSges(6,2) + t44 * mrSges(7,2) + t65;
t56 = -pkin(10) - qJ(4);
t66 = mrSges(4,2) - m(5) * qJ(4) - mrSges(5,3) + m(7) * (-pkin(11) + t56) - mrSges(7,3) + m(6) * t56 - mrSges(6,3);
t57 = sin(qJ(3));
t60 = cos(qJ(3));
t129 = m(5) + m(6) + m(7);
t99 = m(4) + t129;
t136 = pkin(2) * t99 - t143 * t60 - t66 * t57 + mrSges(3,1);
t140 = t47 * mrSges(6,1) + t44 * mrSges(7,1) + t48 * mrSges(6,2) + t45 * mrSges(7,2);
t135 = -t55 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t115 = pkin(4) * t54;
t35 = pkin(5) * t47 + t115;
t121 = -t54 * mrSges(5,1) + t135 + (-m(4) - m(5)) * pkin(9) - m(6) * (pkin(9) + t115) - m(7) * (pkin(9) + t35) - t140;
t132 = -m(7) * pkin(5) - mrSges(6,1);
t122 = -m(7) * t35 - t99 * pkin(9) - (m(6) * pkin(4) + mrSges(5,1)) * t54 + t135;
t59 = sin(qJ(1));
t61 = cos(qJ(2));
t106 = t59 * t61;
t112 = cos(qJ(1));
t101 = pkin(6) + qJ(2);
t88 = sin(t101);
t83 = t88 / 0.2e1;
t102 = pkin(6) - qJ(2);
t89 = sin(t102);
t76 = t83 - t89 / 0.2e1;
t24 = t112 * t76 + t106;
t103 = sin(pkin(6));
t86 = t112 * t103;
t12 = t24 * t60 - t57 * t86;
t58 = sin(qJ(2));
t85 = cos(t101) / 0.2e1;
t90 = cos(t102);
t68 = t90 / 0.2e1 + t85;
t23 = -t112 * t68 + t58 * t59;
t119 = (-t12 * t44 + t23 * t45) * mrSges(7,1) + (-t12 * t45 - t23 * t44) * mrSges(7,2);
t96 = t112 * t61;
t27 = -t59 * t76 + t96;
t91 = t59 * t103;
t16 = t27 * t60 + t57 * t91;
t26 = t112 * t58 + t59 * t68;
t5 = -t16 * t44 + t26 * t45;
t6 = t16 * t45 + t26 * t44;
t118 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t104 = cos(pkin(6));
t34 = t85 - t90 / 0.2e1;
t22 = t104 * t57 - t34 * t60;
t84 = t89 / 0.2e1;
t33 = t83 + t84;
t113 = (-t22 * t44 - t33 * t45) * mrSges(7,1) + (-t22 * t45 + t33 * t44) * mrSges(7,2);
t105 = t112 * pkin(1) + pkin(8) * t91;
t93 = -pkin(1) * t59 + pkin(8) * t86;
t7 = -t16 * t47 + t26 * t48;
t77 = t84 - t88 / 0.2e1;
t11 = t24 * t57 + t60 * t86;
t21 = -t104 * t60 - t34 * t57;
t15 = t27 * t57 - t60 * t91;
t8 = t16 * t48 + t26 * t47;
t1 = [(-t112 * mrSges(2,1) - m(3) * t105 - t27 * mrSges(3,1) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + (-t103 * mrSges(3,3) + mrSges(2,2)) * t59 + t65 * t16 + t122 * t26 + t66 * t15 - t99 * (t27 * pkin(2) + t105)) * g(2) + (t59 * mrSges(2,1) + t112 * mrSges(2,2) - m(3) * t93 + t24 * mrSges(3,1) - mrSges(3,3) * t86 - t143 * t12 + (-t122 + t140) * t23 - t66 * t11 + t99 * (t24 * pkin(2) - t93)) * g(1) (-t121 * t34 - t136 * t33) * g(3) + (t121 * (-t112 * t77 + t106) + t136 * t23) * g(2) + (t121 * (t59 * t77 + t96) + t136 * t26) * g(1) (-t143 * t21 + t22 * t66) * g(3) + (-t11 * t143 + t12 * t66) * g(2) + (-t143 * t15 + t16 * t66) * g(1), t129 * (-g(1) * t15 - g(2) * t11 - g(3) * t21) (-(-t22 * t48 + t33 * t47) * mrSges(6,2) - t113 + t132 * (-t22 * t47 - t33 * t48)) * g(3) + (-(-t12 * t48 - t23 * t47) * mrSges(6,2) - t119 + t132 * (-t12 * t47 + t23 * t48)) * g(2) + (t8 * mrSges(6,2) + t132 * t7 - t118) * g(1), -g(1) * t118 - g(2) * t119 - g(3) * t113];
taug  = t1(:);
