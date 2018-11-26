% Calculate Gravitation load on the joints for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:46
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:45:02
% EndTime: 2018-11-23 18:45:04
% DurationCPUTime: 2.08s
% Computational Cost: add. (5407->193), mult. (5650->259), div. (0->0), fcn. (5605->24), ass. (0->106)
t161 = pkin(7) - qJ(3);
t155 = cos(t161);
t147 = t155 / 0.2e1;
t160 = pkin(7) + qJ(3);
t154 = cos(t160);
t138 = t147 - t154 / 0.2e1;
t165 = sin(pkin(6));
t130 = t138 * t165;
t180 = cos(qJ(1));
t162 = pkin(6) + qJ(2);
t148 = cos(t162) / 0.2e1;
t163 = pkin(6) - qJ(2);
t156 = cos(t163);
t128 = t156 / 0.2e1 + t148;
t178 = sin(qJ(2));
t179 = sin(qJ(1));
t119 = -t128 * t180 + t179 * t178;
t144 = sin(t160) / 0.2e1;
t152 = sin(t161);
t193 = t144 - t152 / 0.2e1;
t100 = cos(qJ(2));
t145 = sin(t162) / 0.2e1;
t153 = sin(t163);
t77 = t145 - t153 / 0.2e1;
t68 = t100 * t179 + t180 * t77;
t99 = cos(qJ(3));
t200 = -t119 * t193 + t68 * t99;
t31 = t130 * t180 - t200;
t166 = cos(pkin(7));
t139 = t166 * t165;
t93 = sin(pkin(7));
t52 = t119 * t93 - t180 * t139;
t95 = sin(qJ(4));
t98 = cos(qJ(4));
t12 = -t31 * t98 + t52 * t95;
t207 = t31 * t95 + t52 * t98;
t198 = -m(5) - m(6);
t97 = cos(qJ(5));
t88 = pkin(5) * t97 + pkin(4);
t92 = qJ(5) + qJ(6);
t89 = sin(t92);
t90 = cos(t92);
t94 = sin(qJ(5));
t191 = m(6) * pkin(4) + m(7) * t88 + t97 * mrSges(6,1) + t90 * mrSges(7,1) - t94 * mrSges(6,2) - t89 * mrSges(7,2) + mrSges(5,1);
t190 = mrSges(5,2) + m(7) * (-pkin(13) - pkin(12)) - mrSges(7,3) - m(6) * pkin(12) - mrSges(6,3);
t112 = t128 * t179 + t178 * t180;
t102 = t112 * t93 + t179 * t139;
t204 = -t89 * mrSges(7,1) - t97 * mrSges(6,2) - t90 * mrSges(7,2);
t192 = m(7) - t198;
t203 = pkin(3) * t192 - t190 * t95 + t191 * t98 + mrSges(4,1);
t195 = mrSges(4,2) - mrSges(5,3);
t202 = t195 + t204;
t70 = t100 * t180 - t179 * t77;
t201 = -t112 * t193 + t70 * t99;
t125 = t145 + t153 / 0.2e1;
t78 = t148 - t156 / 0.2e1;
t199 = t125 * t193 - t78 * t99;
t197 = m(7) * (pkin(5) * t94 + pkin(11));
t196 = t94 * mrSges(6,1);
t194 = -m(7) * pkin(5) - mrSges(6,1);
t123 = t144 + t152 / 0.2e1;
t120 = t123 * t165;
t146 = t154 / 0.2e1;
t127 = t147 + t146;
t96 = sin(qJ(3));
t27 = t119 * t127 + t120 * t180 + t68 * t96;
t187 = t195 - t197;
t186 = t198 * pkin(11) - t196 - t197 + t202;
t183 = (-t12 * t89 + t27 * t90) * mrSges(7,1) + (-t12 * t90 - t27 * t89) * mrSges(7,2);
t33 = t130 * t179 + t201;
t16 = t102 * t95 + t33 * t98;
t32 = t112 * t127 - t120 * t179 + t70 * t96;
t5 = -t16 * t89 + t32 * t90;
t6 = t16 * t90 + t32 * t89;
t182 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t167 = cos(pkin(6));
t45 = t138 * t167 + t199;
t67 = -t125 * t93 + t166 * t167;
t22 = t45 * t98 + t67 * t95;
t44 = -t123 * t167 - t125 * t127 - t78 * t96;
t181 = (-t22 * t89 + t44 * t90) * mrSges(7,1) + (-t22 * t90 - t44 * t89) * mrSges(7,2);
t174 = t93 * t95;
t173 = t93 * t98;
t149 = t165 * t179;
t168 = t180 * pkin(1) + pkin(9) * t149;
t150 = t180 * t165;
t151 = -pkin(1) * t179 + pkin(9) * t150;
t7 = -t16 * t94 + t32 * t97;
t129 = -mrSges(3,2) + (mrSges(4,3) + (m(4) + t192) * pkin(10)) * t93;
t126 = t146 - t155 / 0.2e1;
t121 = t126 * t165;
t111 = -pkin(11) * t192 + t194 * t94 + t202;
t109 = -t68 * pkin(2) - t52 * pkin(10) + t151;
t107 = t31 * pkin(3) + t109;
t106 = t70 * pkin(2) + t102 * pkin(10) + t168;
t105 = t33 * pkin(3) + t106;
t103 = t32 * pkin(11) + t105;
t76 = t125 * pkin(2);
t65 = t112 * pkin(2);
t63 = t119 * pkin(2);
t49 = t125 * t99 + t193 * t78;
t42 = -t112 * t99 - t193 * t70;
t40 = -t119 * t99 - t193 * t68;
t15 = -t102 * t98 + t33 * t95;
t8 = t16 * t97 + t32 * t94;
t1 = [(-t180 * mrSges(2,1) + t179 * mrSges(2,2) - m(3) * t168 - t70 * mrSges(3,1) + t112 * mrSges(3,2) - mrSges(3,3) * t149 - m(4) * t106 - t33 * mrSges(4,1) - t102 * mrSges(4,3) - m(5) * t103 - t16 * mrSges(5,1) - m(6) * (t16 * pkin(4) + t103) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t16 * t88 + t105) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + t187 * t32 + t190 * t15) * g(2) + (-m(3) * t151 - m(4) * t109 - m(7) * t107 + t179 * mrSges(2,1) + t68 * mrSges(3,1) - t31 * mrSges(4,1) + t180 * mrSges(2,2) - t119 * mrSges(3,2) - mrSges(3,3) * t150 + t52 * mrSges(4,3) + t198 * (-pkin(11) * t27 + t107) + t191 * t12 - (t187 - t196 + t204) * t27 + t190 * t207) * g(1) (-t125 * mrSges(3,1) - m(4) * t76 - t49 * mrSges(4,1) + t129 * t78 - t191 * (-t174 * t78 + t49 * t98) + t111 * (t125 * t96 - t127 * t78) + t190 * (t173 * t78 + t49 * t95)) * g(3) + (t119 * mrSges(3,1) + m(4) * t63 - t40 * mrSges(4,1) - t129 * t68 - t191 * (t174 * t68 + t40 * t98) + t111 * (-t119 * t96 + t127 * t68) + t190 * (-t173 * t68 + t40 * t95)) * g(2) + (t112 * mrSges(3,1) + m(4) * t65 - t42 * mrSges(4,1) - t129 * t70 - t191 * (t174 * t70 + t42 * t98) + t111 * (-t112 * t96 + t127 * t70) + t190 * (-t173 * t70 + t42 * t95)) * g(1) + (-g(2) * (t40 * pkin(3) - t63) - g(1) * (t42 * pkin(3) - t65) - g(3) * (t49 * pkin(3) + t76)) * t192 (t186 * (-t126 * t167 + t199) + t203 * t44) * g(3) + (t186 * (t121 * t180 + t200) + t203 * t27) * g(2) + (t186 * (-t121 * t179 + t201) + t203 * t32) * g(1) (t190 * t22 - t191 * (-t45 * t95 + t67 * t98)) * g(3) + (t190 * t12 - t191 * t207) * g(2) + (t15 * t191 + t16 * t190) * g(1) (-(-t22 * t97 - t44 * t94) * mrSges(6,2) - t181 + t194 * (-t22 * t94 + t44 * t97)) * g(3) + (-(-t12 * t97 - t27 * t94) * mrSges(6,2) - t183 + t194 * (-t12 * t94 + t27 * t97)) * g(2) + (mrSges(6,2) * t8 + t194 * t7 - t182) * g(1), -g(1) * t182 - g(2) * t183 - g(3) * t181];
taug  = t1(:);
