% Calculate Gravitation load on the joints for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2018-11-23 16:41
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:40:47
% EndTime: 2018-11-23 16:40:49
% DurationCPUTime: 2.35s
% Computational Cost: add. (8984->190), mult. (9130->262), div. (0->0), fcn. (9078->30), ass. (0->111)
t100 = sin(qJ(5));
t103 = cos(qJ(5));
t104 = cos(qJ(4));
t105 = cos(qJ(3));
t106 = cos(qJ(1));
t179 = sin(pkin(6));
t168 = t106 * t179;
t172 = pkin(6) - pkin(14);
t150 = cos(t172) / 0.2e1;
t171 = pkin(6) + pkin(14);
t161 = cos(t171);
t133 = t150 + t161 / 0.2e1;
t177 = sin(pkin(14));
t193 = sin(qJ(1));
t82 = -t106 * t133 + t193 * t177;
t149 = sin(t171) / 0.2e1;
t160 = sin(t172);
t132 = t149 - t160 / 0.2e1;
t180 = cos(pkin(14));
t83 = t106 * t132 + t193 * t180;
t175 = pkin(7) + qJ(3);
t154 = sin(t175) / 0.2e1;
t176 = pkin(7) - qJ(3);
t164 = sin(t176);
t85 = t154 - t164 / 0.2e1;
t157 = cos(t175) / 0.2e1;
t167 = cos(t176);
t86 = t157 - t167 / 0.2e1;
t58 = t83 * t105 + t86 * t168 - t82 * t85;
t136 = t154 + t164 / 0.2e1;
t126 = t136 * t179;
t139 = t167 / 0.2e1 + t157;
t192 = sin(qJ(3));
t60 = t106 * t126 + t82 * t139 + t83 * t192;
t173 = pkin(8) + qJ(4);
t153 = sin(t173) / 0.2e1;
t174 = pkin(8) - qJ(4);
t163 = sin(t174);
t84 = t153 - t163 / 0.2e1;
t115 = t58 * t104 - t60 * t84;
t181 = cos(pkin(7));
t151 = t181 * t179;
t178 = sin(pkin(7));
t129 = -t106 * t151 + t82 * t178;
t166 = cos(t174);
t156 = t166 / 0.2e1;
t165 = cos(t173);
t146 = t156 - t165 / 0.2e1;
t20 = t129 * t146 + t115;
t97 = sin(pkin(8));
t98 = cos(pkin(8));
t44 = t129 * t98 + t60 * t97;
t4 = t100 * t44 + t103 * t20;
t208 = -t100 * t20 + t103 * t44;
t102 = cos(qJ(6));
t201 = m(6) + m(7);
t147 = -pkin(12) * t201 + mrSges(5,2) - mrSges(6,3);
t99 = sin(qJ(6));
t207 = -t99 * mrSges(7,1) - t102 * mrSges(7,2) + t147;
t206 = m(5) + t201;
t199 = m(7) * pkin(5) + mrSges(7,1) * t102 - mrSges(7,2) * t99 + mrSges(6,1);
t162 = -m(7) * pkin(13) + mrSges(6,2) - mrSges(7,3);
t121 = t106 * t177 + t193 * t133;
t122 = t106 * t180 - t193 * t132;
t112 = t121 * t139 + t122 * t192 - t193 * t126;
t203 = t121 * t178 + t193 * t151;
t204 = t112 * t97 + t203 * t98;
t202 = pkin(4) * t201 - t100 * t162 + t199 * t103 + mrSges(5,1);
t101 = sin(qJ(4));
t135 = t153 + t163 / 0.2e1;
t155 = t165 / 0.2e1;
t138 = t156 + t155;
t19 = t101 * t58 - t129 * t135 + t60 * t138;
t158 = t179 * t193;
t185 = t106 * pkin(1) + qJ(2) * t158;
t184 = t100 * t97;
t183 = t103 * t97;
t182 = cos(pkin(6));
t159 = -t193 * pkin(1) + qJ(2) * t168;
t140 = -mrSges(4,2) + (t206 * pkin(11) + mrSges(5,3)) * t97;
t137 = t155 - t166 / 0.2e1;
t134 = t150 - t161 / 0.2e1;
t131 = t149 + t160 / 0.2e1;
t127 = -t83 * pkin(2) - pkin(10) * t129 + t159;
t124 = -t58 * pkin(3) - pkin(11) * t44 + t127;
t120 = t131 * t178 - t182 * t181;
t118 = t122 * pkin(2) + pkin(10) * t203 + t185;
t68 = t134 * t105 + t131 * t85 - t182 * t86;
t62 = t122 * t105 - t121 * t85 - t86 * t158;
t114 = -t131 * t139 + t134 * t192 - t182 * t136;
t113 = t68 * t104 - t114 * t84;
t110 = t62 * t104 - t112 * t84;
t109 = t62 * pkin(3) + pkin(11) * t204 + t118;
t25 = t146 * t203 + t110;
t108 = t25 * pkin(4) + t109;
t67 = t114 * pkin(3);
t56 = t112 * pkin(3);
t54 = t60 * pkin(3);
t53 = t114 * t97 - t120 * t98;
t41 = -t114 * t104 - t68 * t84;
t37 = -t120 * t146 + t113;
t36 = t101 * t68 + t114 * t138 + t120 * t135;
t34 = -t112 * t104 - t62 * t84;
t32 = -t60 * t104 - t58 * t84;
t24 = t101 * t62 + t112 * t138 - t135 * t203;
t14 = t100 * t53 + t103 * t37;
t8 = t100 * t204 + t25 * t103;
t7 = t100 * t25 - t103 * t204;
t2 = t102 * t8 + t24 * t99;
t1 = t102 * t24 - t8 * t99;
t3 = [(-t106 * mrSges(2,1) + t193 * mrSges(2,2) - m(3) * t185 - t122 * mrSges(3,1) + t121 * mrSges(3,2) - mrSges(3,3) * t158 - m(4) * t118 - t62 * mrSges(4,1) + t112 * mrSges(4,2) - t203 * mrSges(4,3) - m(5) * t109 - t25 * mrSges(5,1) - t204 * mrSges(5,3) - m(6) * t108 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t108) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t162 * t7 + t147 * t24) * g(2) + (t193 * mrSges(2,1) - m(3) * t159 + t83 * mrSges(3,1) - t82 * mrSges(3,2) - m(4) * t127 + t58 * mrSges(4,1) - t60 * mrSges(4,2) + t129 * mrSges(4,3) - m(5) * t124 + t20 * mrSges(5,1) + t44 * mrSges(5,3) + t162 * t208 + t199 * t4 - t207 * t19 + (-t179 * mrSges(3,3) + mrSges(2,2)) * t106 + t201 * (pkin(4) * t20 - t124)) * g(1) (-g(1) * t158 + g(2) * t168 - g(3) * t182) * (m(3) + m(4) + t206) (t114 * mrSges(4,1) + m(5) * t67 - t41 * mrSges(5,1) - t140 * t68 - t199 * (t103 * t41 + t184 * t68) + t207 * (-t114 * t101 + t138 * t68) + t162 * (t100 * t41 - t183 * t68) - t201 * (t41 * pkin(4) - t67)) * g(3) + (t60 * mrSges(4,1) + m(5) * t54 - t32 * mrSges(5,1) + t162 * (t100 * t32 - t183 * t58) - t140 * t58 - t199 * (t103 * t32 + t184 * t58) + t207 * (-t60 * t101 + t138 * t58) - t201 * (t32 * pkin(4) - t54)) * g(2) + (t112 * mrSges(4,1) + m(5) * t56 - t34 * mrSges(5,1) - t140 * t62 - t199 * (t103 * t34 + t184 * t62) + t207 * (-t112 * t101 + t138 * t62) + t162 * (t100 * t34 - t183 * t62) - t201 * (t34 * pkin(4) - t56)) * g(1) (t207 * (t120 * t137 + t113) + t202 * t36) * g(3) + (t207 * (-t129 * t137 + t115) + t202 * t19) * g(2) + (t207 * (-t137 * t203 + t110) + t202 * t24) * g(1) (t162 * t14 - t199 * (-t100 * t37 + t103 * t53)) * g(3) + (t162 * t4 - t199 * t208) * g(2) + (t162 * t8 + t199 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t102 * t19 - t4 * t99) * mrSges(7,1) + (-t102 * t4 - t19 * t99) * mrSges(7,2)) - g(3) * ((t102 * t36 - t14 * t99) * mrSges(7,1) + (-t102 * t14 - t36 * t99) * mrSges(7,2))];
taug  = t3(:);
