% Calculate Gravitation load on the joints for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:31
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:30:50
% EndTime: 2018-11-23 18:30:52
% DurationCPUTime: 1.70s
% Computational Cost: add. (1779->171), mult. (1894->220), div. (0->0), fcn. (1882->16), ass. (0->83)
t184 = mrSges(6,2) + mrSges(7,2);
t192 = mrSges(6,1) + mrSges(7,1);
t94 = sin(qJ(5));
t98 = cos(qJ(5));
t198 = -t184 * t94 + t98 * t192 + mrSges(5,1);
t190 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t93 = -qJ(6) - pkin(11);
t107 = -m(6) * pkin(11) + m(7) * t93 + t190;
t100 = cos(qJ(2));
t173 = cos(qJ(1));
t134 = pkin(6) + qJ(2);
t115 = sin(t134) / 0.2e1;
t135 = pkin(6) - qJ(2);
t121 = sin(t135);
t71 = t115 - t121 / 0.2e1;
t97 = sin(qJ(1));
t111 = t97 * t100 + t173 * t71;
t92 = sin(pkin(6));
t130 = t92 * t173;
t91 = qJ(3) + qJ(4);
t88 = sin(t91);
t89 = cos(t91);
t31 = t111 * t88 + t130 * t89;
t32 = t111 * t89 - t88 * t130;
t196 = t190 * t32 + t198 * t31;
t110 = t100 * t173 - t97 * t71;
t155 = t92 * t97;
t35 = t110 * t88 - t155 * t89;
t36 = t110 * t89 + t155 * t88;
t195 = t190 * t36 + t198 * t35;
t136 = cos(pkin(6));
t116 = cos(t134) / 0.2e1;
t122 = cos(t135);
t72 = t116 - t122 / 0.2e1;
t54 = -t136 * t89 - t72 * t88;
t55 = t136 * t88 - t72 * t89;
t194 = t190 * t55 + t198 * t54;
t95 = sin(qJ(3));
t99 = cos(qJ(3));
t188 = t136 * t99 + t72 * t95;
t39 = -t110 * t95 + t99 * t155;
t86 = pkin(5) * t98 + pkin(4);
t187 = -m(4) * pkin(2) - t99 * mrSges(4,1) + t95 * mrSges(4,2) - mrSges(3,1) + (-m(6) * pkin(4) - m(7) * t86 - mrSges(5,1)) * t89 + t107 * t88;
t102 = t122 / 0.2e1 + t116;
t96 = sin(qJ(2));
t60 = -t102 * t173 + t96 * t97;
t186 = t32 * t94 - t60 * t98;
t185 = -t32 * t98 - t60 * t94;
t174 = m(7) * pkin(5);
t183 = -m(5) - m(6) - m(7);
t182 = -t174 - t192;
t181 = -m(4) * pkin(9) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t164 = t111 * t94;
t162 = t110 * t94;
t161 = t72 * t94;
t157 = t89 * t94;
t156 = t89 * t98;
t152 = -t31 * t86 - t32 * t93;
t151 = -t35 * t86 - t36 * t93;
t144 = -t54 * t86 - t55 * t93;
t137 = t173 * pkin(1) + pkin(8) * t155;
t133 = t95 * t155;
t129 = -pkin(1) * t97 + pkin(8) * t130;
t128 = -t31 * pkin(4) + t32 * pkin(11);
t127 = -t35 * pkin(4) + pkin(11) * t36;
t126 = -t54 * pkin(4) + pkin(11) * t55;
t78 = t95 * t130;
t124 = -t111 * t99 + t78;
t120 = t39 * pkin(3);
t119 = t188 * pkin(3);
t101 = -pkin(10) - pkin(9);
t63 = t102 * t97 + t173 * t96;
t87 = pkin(3) * t99 + pkin(2);
t118 = pkin(3) * t133 - t63 * t101 + t110 * t87 + t137;
t5 = -t36 * t94 + t63 * t98;
t112 = pkin(3) * t78 + t60 * t101 - t111 * t87 + t129;
t108 = t111 * t95 + t130 * t99;
t106 = t174 * t94 - t181;
t105 = t108 * pkin(3);
t70 = t115 + t121 / 0.2e1;
t40 = t110 * t99 + t133;
t6 = t36 * t98 + t63 * t94;
t1 = [(-t173 * mrSges(2,1) - m(3) * t137 - t110 * mrSges(3,1) - m(4) * (pkin(2) * t110 + t137) - t40 * mrSges(4,1) - t39 * mrSges(4,2) - m(5) * t118 - t36 * mrSges(5,1) - m(6) * (pkin(4) * t36 + t118) - m(7) * (t36 * t86 + t118) + (-mrSges(3,3) * t92 + mrSges(2,2)) * t97 - t192 * t6 - t184 * t5 - t106 * t63 + t107 * t35) * g(2) + (t97 * mrSges(2,1) + t173 * mrSges(2,2) - m(3) * t129 + t111 * mrSges(3,1) - mrSges(3,3) * t130 - m(4) * (-pkin(2) * t111 + t129) - t124 * mrSges(4,1) - t108 * mrSges(4,2) - m(5) * t112 + t32 * mrSges(5,1) - m(6) * (-pkin(4) * t32 + t112) - m(7) * (-t32 * t86 + t112) - t192 * t185 - t184 * t186 + t106 * t60 - t107 * t31) * g(1) (t161 * t174 - t192 * (t156 * t70 - t161) - t184 * (-t157 * t70 - t72 * t98) + t183 * (t72 * t101 + t70 * t87) - t181 * t72 + t187 * t70) * g(3) + (-t164 * t174 - t192 * (-t156 * t60 + t164) - t184 * (t111 * t98 + t157 * t60) + t183 * (-t101 * t111 - t60 * t87) + t181 * t111 - t187 * t60) * g(2) + (-t162 * t174 - t184 * (t110 * t98 + t157 * t63) + t183 * (-t101 * t110 - t63 * t87) - t192 * (-t156 * t63 + t162) + t181 * t110 - t187 * t63) * g(1) (-t188 * mrSges(4,1) - (-t136 * t95 + t72 * t99) * mrSges(4,2) - m(5) * t119 - m(6) * (t119 + t126) - m(7) * (t119 + t144) + t194) * g(3) + (mrSges(4,1) * t108 - mrSges(4,2) * t124 + m(5) * t105 - m(6) * (-t105 + t128) - m(7) * (-t105 + t152) + t196) * g(2) + (-mrSges(4,1) * t39 + mrSges(4,2) * t40 - m(5) * t120 - m(6) * (t120 + t127) - m(7) * (t120 + t151) + t195) * g(1) (-m(6) * t126 - m(7) * t144 + t194) * g(3) + (-m(6) * t128 - m(7) * t152 + t196) * g(2) + (-m(6) * t127 - m(7) * t151 + t195) * g(1) (-t184 * (-t55 * t98 + t70 * t94) + t182 * (-t55 * t94 - t70 * t98)) * g(3) + (-t182 * t186 - t184 * t185) * g(2) + (t182 * t5 + t184 * t6) * g(1) (-g(1) * t35 - g(2) * t31 - g(3) * t54) * m(7)];
taug  = t1(:);
