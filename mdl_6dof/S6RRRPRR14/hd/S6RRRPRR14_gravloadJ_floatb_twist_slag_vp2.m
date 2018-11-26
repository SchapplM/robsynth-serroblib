% Calculate Gravitation load on the joints for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2018-11-23 18:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR14_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:02:02
% EndTime: 2018-11-23 18:02:04
% DurationCPUTime: 1.34s
% Computational Cost: add. (1600->150), mult. (1920->183), div. (0->0), fcn. (1919->16), ass. (0->88)
t155 = m(7) * pkin(5);
t59 = sin(qJ(5));
t108 = t59 * t155;
t146 = mrSges(4,2) - mrSges(5,3);
t57 = qJ(5) + qJ(6);
t54 = sin(t57);
t55 = cos(t57);
t63 = cos(qJ(5));
t148 = -t54 * mrSges(7,1) - t63 * mrSges(6,2) - t55 * mrSges(7,2);
t156 = t59 * mrSges(6,1) + t108 - t146 - t148;
t154 = mrSges(4,1) - mrSges(5,2) + mrSges(6,3) - m(7) * (-pkin(11) - pkin(10)) + mrSges(7,3);
t152 = t63 * mrSges(6,1) + t55 * mrSges(7,1) - t59 * mrSges(6,2) - t54 * mrSges(7,2);
t141 = m(5) + m(6) + m(7);
t151 = qJ(4) * t141;
t140 = mrSges(3,2) - mrSges(5,1) - mrSges(4,3);
t53 = pkin(5) * t63 + pkin(4);
t130 = -m(6) * (pkin(4) + pkin(9)) + t140 - m(7) * (pkin(9) + t53) - t152;
t60 = sin(qJ(3));
t64 = cos(qJ(3));
t132 = t154 * t64 + t156 * t60 + mrSges(3,1);
t71 = -m(6) * pkin(10) - t154;
t147 = pkin(3) * t141 - t71;
t111 = qJ(4) * t60;
t122 = cos(qJ(1));
t61 = sin(qJ(2));
t62 = sin(qJ(1));
t106 = pkin(6) + qJ(2);
t88 = cos(t106) / 0.2e1;
t107 = pkin(6) - qJ(2);
t92 = cos(t107);
t69 = t92 / 0.2e1 + t88;
t32 = -t122 * t69 + t61 * t62;
t121 = t32 * t64;
t145 = -pkin(3) * t121 - t32 * t111;
t35 = t122 * t61 + t62 * t69;
t120 = t35 * t64;
t144 = -pkin(3) * t120 - t35 * t111;
t90 = sin(t106);
t86 = t90 / 0.2e1;
t91 = sin(t107);
t87 = t91 / 0.2e1;
t43 = t86 + t87;
t119 = t43 * t64;
t143 = pkin(3) * t119 + t43 * t111;
t142 = -mrSges(6,1) - t155;
t139 = -m(6) * pkin(4) - m(7) * t53 - pkin(9) * (m(4) + t141) + t140;
t135 = -t151 - t156;
t58 = sin(pkin(6));
t100 = t58 * t122;
t65 = cos(qJ(2));
t115 = t62 * t65;
t77 = t86 - t91 / 0.2e1;
t33 = t122 * t77 + t115;
t15 = t100 * t64 + t33 * t60;
t127 = (t15 * t55 - t32 * t54) * mrSges(7,1) + (-t15 * t54 - t32 * t55) * mrSges(7,2);
t118 = t58 * t62;
t99 = t122 * t65;
t36 = -t62 * t77 + t99;
t19 = -t118 * t64 + t36 * t60;
t5 = t19 * t55 - t35 * t54;
t6 = t19 * t54 + t35 * t55;
t126 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t110 = cos(pkin(6));
t44 = t88 - t92 / 0.2e1;
t30 = -t110 * t64 - t44 * t60;
t123 = (t30 * t55 + t43 * t54) * mrSges(7,1) + (-t30 * t54 + t43 * t55) * mrSges(7,2);
t112 = t122 * pkin(1) + pkin(8) * t118;
t26 = t32 * pkin(2);
t105 = -t26 + t145;
t28 = t35 * pkin(2);
t104 = -t28 + t144;
t103 = t36 * pkin(2) + t112;
t42 = t43 * pkin(2);
t102 = t42 + t143;
t97 = -t62 * pkin(1) + pkin(8) * t100;
t78 = t87 - t90 / 0.2e1;
t34 = -t122 * t78 + t115;
t96 = t34 * pkin(9) - t26;
t37 = t62 * t78 + t99;
t95 = t37 * pkin(9) - t28;
t94 = -t44 * pkin(9) + t42;
t16 = -t60 * t100 + t33 * t64;
t89 = -t33 * pkin(2) + t97;
t7 = t19 * t63 - t35 * t59;
t76 = t146 - t151;
t20 = t118 * t60 + t36 * t64;
t8 = t19 * t59 + t35 * t63;
t1 = [(-t122 * mrSges(2,1) - m(3) * t112 - t36 * mrSges(3,1) - m(4) * t103 - t8 * mrSges(6,1) - t7 * mrSges(6,2) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + (-mrSges(3,3) * t58 + mrSges(2,2)) * t62 + (t76 - t108) * t19 + t71 * t20 + t139 * t35 - t141 * (t20 * pkin(3) + t103)) * g(2) + (t62 * mrSges(2,1) + t122 * mrSges(2,2) - m(3) * t97 + t33 * mrSges(3,1) - mrSges(3,3) * t100 - m(4) * t89 - t71 * t16 - (t142 * t59 + t148 + t76) * t15 + (-t139 + t152) * t32 + t141 * (pkin(3) * t16 - t89)) * g(1) (-m(4) * t94 - m(5) * (t94 + t143) - m(6) * (pkin(10) * t119 + t102) - m(7) * t102 - t130 * t44 - t132 * t43) * g(3) + (-m(4) * t96 - m(5) * (t96 + t145) - m(6) * (-pkin(10) * t121 + t105) - m(7) * t105 + t130 * t34 + t132 * t32) * g(2) + (-m(4) * t95 - m(5) * (t95 + t144) - m(6) * (-pkin(10) * t120 + t104) - m(7) * t104 + t130 * t37 + t132 * t35) * g(1) (t135 * (t110 * t60 - t44 * t64) + t147 * t30) * g(3) + (t135 * t16 + t147 * t15) * g(2) + (t135 * t20 + t147 * t19) * g(1), t141 * (-g(1) * t19 - g(2) * t15 - g(3) * t30) (-(-t30 * t59 + t43 * t63) * mrSges(6,2) - t123 + t142 * (t30 * t63 + t43 * t59)) * g(3) + (-(-t15 * t59 - t32 * t63) * mrSges(6,2) - t127 + t142 * (t15 * t63 - t32 * t59)) * g(2) + (t8 * mrSges(6,2) + t142 * t7 - t126) * g(1), -g(1) * t126 - g(2) * t127 - g(3) * t123];
taug  = t1(:);
