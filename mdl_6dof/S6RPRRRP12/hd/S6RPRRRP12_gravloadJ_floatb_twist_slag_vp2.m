% Calculate Gravitation load on the joints for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:31
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:31:23
% EndTime: 2018-11-23 16:31:25
% DurationCPUTime: 1.59s
% Computational Cost: add. (4459->160), mult. (4771->222), div. (0->0), fcn. (4786->22), ass. (0->100)
t155 = pkin(6) - pkin(12);
t126 = cos(t155) / 0.2e1;
t154 = pkin(6) + pkin(12);
t141 = cos(t154);
t114 = t126 + t141 / 0.2e1;
t159 = sin(pkin(12));
t172 = sin(qJ(1));
t95 = cos(qJ(1));
t108 = -t114 * t95 + t172 * t159;
t160 = cos(pkin(7));
t87 = sin(pkin(6));
t146 = t87 * t160;
t86 = sin(pkin(7));
t186 = t108 * t86 - t95 * t146;
t157 = pkin(7) - qJ(3);
t145 = cos(t157);
t136 = t145 / 0.2e1;
t156 = pkin(7) + qJ(3);
t144 = cos(t156);
t123 = t136 - t144 / 0.2e1;
t119 = t87 * t123;
t134 = sin(t156) / 0.2e1;
t143 = sin(t157);
t178 = t134 - t143 / 0.2e1;
t125 = sin(t154) / 0.2e1;
t140 = sin(t155);
t76 = t125 - t140 / 0.2e1;
t88 = cos(pkin(12));
t70 = t172 * t88 + t95 * t76;
t94 = cos(qJ(3));
t184 = -t108 * t178 + t70 * t94;
t45 = t119 * t95 - t184;
t90 = sin(qJ(4));
t93 = cos(qJ(4));
t20 = t186 * t90 - t45 * t93;
t115 = t134 + t143 / 0.2e1;
t111 = t87 * t115;
t135 = t144 / 0.2e1;
t118 = t136 + t135;
t91 = sin(qJ(3));
t41 = t108 * t118 + t111 * t95 + t70 * t91;
t89 = sin(qJ(5));
t92 = cos(qJ(5));
t1 = t20 * t89 - t41 * t92;
t190 = t20 * t92 + t41 * t89;
t142 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t139 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t19 = t186 * t93 + t45 * t90;
t174 = -m(6) - m(7);
t158 = -m(5) + t174;
t105 = t114 * t172 + t159 * t95;
t187 = t105 * t86 + t172 * t146;
t71 = -t172 * t76 + t95 * t88;
t185 = -t105 * t178 + t71 * t94;
t113 = t125 + t140 / 0.2e1;
t77 = t126 - t141 / 0.2e1;
t183 = t113 * t178 + t77 * t94;
t181 = mrSges(4,2) - mrSges(5,3);
t180 = mrSges(6,3) + mrSges(7,2);
t179 = t93 * mrSges(5,1) - t90 * mrSges(5,2) + mrSges(4,1);
t177 = mrSges(5,2) - t180;
t176 = t139 * t89 - t142 * t92 - mrSges(5,1);
t173 = pkin(4) * t93;
t171 = t41 * t90;
t46 = t105 * t118 - t111 * t172 + t71 * t91;
t170 = t46 * t90;
t161 = cos(pkin(6));
t54 = -t113 * t118 - t115 * t161 + t77 * t91;
t169 = t54 * t90;
t165 = t87 * t95;
t164 = t89 * t93;
t163 = t92 * t93;
t153 = t87 * t172;
t162 = t95 * pkin(1) + qJ(2) * t153;
t117 = t135 - t145 / 0.2e1;
t112 = t87 * t117;
t43 = t112 * t95 + t184;
t152 = -t41 * pkin(3) + pkin(10) * t43;
t48 = -t112 * t172 + t185;
t151 = -t46 * pkin(3) + pkin(10) * t48;
t56 = -t117 * t161 + t183;
t150 = -t54 * pkin(3) + pkin(10) * t56;
t137 = -pkin(1) * t172 + qJ(2) * t165;
t121 = pkin(10) * t158 + t181;
t120 = pkin(11) * t174 + t177;
t102 = -t70 * pkin(2) - t186 * pkin(9) + t137;
t101 = t45 * pkin(3) + t102;
t100 = t71 * pkin(2) + t187 * pkin(9) + t162;
t47 = t119 * t172 + t185;
t98 = t47 * pkin(3) + t100;
t69 = -t113 * t86 + t160 * t161;
t55 = t123 * t161 + t183;
t27 = t55 * t93 + t69 * t90;
t26 = -t55 * t90 + t69 * t93;
t24 = t187 * t90 + t47 * t93;
t23 = -t187 * t93 + t47 * t90;
t11 = t27 * t89 - t54 * t92;
t6 = t24 * t92 + t46 * t89;
t5 = t24 * t89 - t46 * t92;
t2 = [(-m(3) * t162 - m(4) * t100 - m(5) * t98 - t95 * mrSges(2,1) - t71 * mrSges(3,1) - t47 * mrSges(4,1) - t24 * mrSges(5,1) + t172 * mrSges(2,2) + t105 * mrSges(3,2) - mrSges(3,3) * t153 - t187 * mrSges(4,3) + t120 * t23 + t121 * t46 + t139 * t5 - t142 * t6 + t174 * (t24 * pkin(4) + t98)) * g(2) + (t172 * mrSges(2,1) + t95 * mrSges(2,2) - m(3) * t137 + t70 * mrSges(3,1) - t108 * mrSges(3,2) - mrSges(3,3) * t165 - m(4) * t102 - t45 * mrSges(4,1) + t186 * mrSges(4,3) - m(5) * t101 + t20 * mrSges(5,1) - t121 * t41 + t142 * t190 - t139 * t1 + t120 * t19 + t174 * (-pkin(4) * t20 + t101)) * g(1) ((t172 * g(1) - t95 * g(2)) * t87 + t161 * g(3)) * (t158 - m(3) - m(4)) (-m(5) * t150 + t181 * t56 + t179 * t54 + t174 * (-pkin(11) * t169 - t54 * t173 + t150) - t142 * (-t163 * t54 + t56 * t89) + t139 * (-t164 * t54 - t56 * t92) + t180 * t169) * g(3) + (-m(5) * t152 + t174 * (-pkin(11) * t171 - t41 * t173 + t152) - t142 * (-t163 * t41 + t43 * t89) + t139 * (-t164 * t41 - t43 * t92) + t181 * t43 + t179 * t41 + t180 * t171) * g(2) + (-m(5) * t151 + t174 * (-pkin(11) * t170 - t46 * t173 + t151) + t139 * (-t164 * t46 - t48 * t92) + t181 * t48 + t179 * t46 + t180 * t170 - t142 * (-t163 * t46 + t48 * t89)) * g(1) (t174 * (t26 * pkin(4) + pkin(11) * t27) + t177 * t27 + t176 * t26) * g(3) + (t174 * (t19 * pkin(4) + pkin(11) * t20) + t177 * t20 + t176 * t19) * g(2) + (t174 * (-t23 * pkin(4) + pkin(11) * t24) + t177 * t24 - t176 * t23) * g(1) (t139 * (t27 * t92 + t54 * t89) + t142 * t11) * g(3) + (t142 * t1 + t139 * t190) * g(2) + (t139 * t6 + t142 * t5) * g(1) (-g(1) * t5 - g(2) * t1 - g(3) * t11) * m(7)];
taug  = t2(:);
