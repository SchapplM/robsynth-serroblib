% Calculate Gravitation load on the joints for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:59:41
% EndTime: 2019-03-10 02:59:47
% DurationCPUTime: 2.07s
% Computational Cost: add. (1772->199), mult. (4850->296), div. (0->0), fcn. (6195->14), ass. (0->98)
t112 = sin(qJ(5));
t115 = cos(qJ(5));
t113 = sin(qJ(4));
t116 = cos(qJ(4));
t111 = sin(pkin(7));
t114 = sin(qJ(3));
t117 = cos(qJ(1));
t177 = sin(pkin(6));
t158 = t117 * t177;
t178 = cos(pkin(7));
t159 = t114 * t178;
t190 = cos(qJ(3));
t179 = cos(pkin(6));
t191 = cos(qJ(2));
t152 = t179 * t191;
t188 = sin(qJ(2));
t189 = sin(qJ(1));
t94 = -t117 * t152 + t189 * t188;
t151 = t179 * t188;
t95 = t117 * t151 + t189 * t191;
t54 = -t111 * t114 * t158 - t94 * t159 + t95 * t190;
t142 = t178 * t177;
t134 = t117 * t142;
t79 = t94 * t111 - t134;
t22 = t113 * t79 + t116 * t54;
t148 = t190 * t177;
t150 = t178 * t190;
t53 = t117 * t111 * t148 + t95 * t114 + t94 * t150;
t1 = t112 * t22 - t53 * t115;
t206 = t112 * t53 + t115 * t22;
t155 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t154 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t125 = t117 * t188 + t189 * t152;
t203 = t125 * t111 + t189 * t142;
t201 = m(6) + m(7);
t170 = -m(5) - t201;
t202 = ((m(4) - t170) * pkin(10) + mrSges(4,3)) * t111 - mrSges(3,2);
t21 = -t113 * t54 + t116 * t79;
t200 = mrSges(4,2) - mrSges(5,3);
t199 = mrSges(6,3) + mrSges(7,2);
t198 = t116 * mrSges(5,1) - mrSges(5,2) * t113 + mrSges(4,1);
t147 = t177 * t189;
t197 = -t111 * t147 + t125 * t178;
t196 = t179 * t111 + t191 * t142;
t195 = mrSges(5,2) - t199;
t194 = t154 * t112 - t155 * t115 - mrSges(5,1);
t187 = pkin(4) * t116;
t64 = -t95 * t159 - t94 * t190;
t89 = t94 * pkin(2);
t186 = t64 * pkin(3) - t89;
t96 = t117 * t191 - t189 * t151;
t66 = -t125 * t190 - t96 * t159;
t91 = t125 * pkin(2);
t185 = t66 * pkin(3) - t91;
t184 = t113 * t53;
t57 = t114 * t96 + t197 * t190;
t183 = t113 * t57;
t146 = t177 * t188;
t77 = t114 * t146 - t196 * t190;
t182 = t113 * t77;
t176 = t111 * t113;
t175 = t111 * t116;
t174 = t112 * t116;
t173 = t115 * t116;
t137 = t111 * t146;
t149 = t191 * t177;
t172 = pkin(2) * t149 + pkin(10) * t137;
t171 = t117 * pkin(1) + pkin(9) * t147;
t131 = t188 * t142;
t88 = -t114 * t131 + t191 * t148;
t167 = t88 * pkin(3) + t172;
t165 = -t53 * pkin(3) + pkin(11) * t54;
t58 = -t197 * t114 + t96 * t190;
t164 = -t57 * pkin(3) + pkin(11) * t58;
t78 = t196 * t114 + t190 * t146;
t163 = -t77 * pkin(3) + pkin(11) * t78;
t153 = -t189 * pkin(1) + pkin(9) * t158;
t135 = t170 * pkin(11) + t200;
t133 = -t201 * pkin(12) + t195;
t132 = -t95 * pkin(2) + pkin(10) * t134 + t153;
t128 = -pkin(3) * t54 + t132;
t121 = t96 * pkin(2) + t203 * pkin(10) + t171;
t120 = t58 * pkin(3) + t121;
t93 = -t111 * t149 + t179 * t178;
t87 = t114 * t149 + t190 * t131;
t69 = t113 * t137 + t88 * t116;
t65 = -t125 * t114 + t96 * t150;
t63 = -t114 * t94 + t95 * t150;
t52 = t113 * t93 + t116 * t78;
t51 = -t113 * t78 + t116 * t93;
t34 = t116 * t66 + t96 * t176;
t32 = t116 * t64 + t95 * t176;
t26 = t113 * t203 + t58 * t116;
t25 = t113 * t58 - t116 * t203;
t15 = t112 * t52 - t77 * t115;
t6 = t112 * t57 + t115 * t26;
t5 = t112 * t26 - t57 * t115;
t2 = [(-m(3) * t171 - m(4) * t121 - m(5) * t120 - t117 * mrSges(2,1) - t96 * mrSges(3,1) - t58 * mrSges(4,1) - t26 * mrSges(5,1) + t189 * mrSges(2,2) + t125 * mrSges(3,2) - mrSges(3,3) * t147 - t203 * mrSges(4,3) + t133 * t25 + t135 * t57 + t154 * t5 - t155 * t6 - t201 * (t26 * pkin(4) + t120)) * g(2) + (t189 * mrSges(2,1) - m(3) * t153 + t95 * mrSges(3,1) - m(4) * t132 + t54 * mrSges(4,1) - m(5) * t128 + t22 * mrSges(5,1) + t202 * t94 - t135 * t53 + t155 * t206 - t154 * t1 + (-t177 * mrSges(3,3) - mrSges(4,3) * t142 + mrSges(2,2)) * t117 + t133 * t21 + t201 * (pkin(4) * t22 - t128)) * g(1) (-mrSges(3,1) * t149 + mrSges(3,2) * t146 - m(4) * t172 - t88 * mrSges(4,1) - mrSges(4,3) * t137 - m(5) * t167 - t69 * mrSges(5,1) + t135 * t87 - t155 * (t112 * t87 + t115 * t69) + t154 * (t112 * t69 - t87 * t115) + t133 * (t113 * t88 - t116 * t137) - t201 * (t69 * pkin(4) + t167)) * g(3) + (mrSges(3,1) * t94 + m(4) * t89 - t64 * mrSges(4,1) - m(5) * t186 - t32 * mrSges(5,1) - t202 * t95 + t135 * t63 - t155 * (t112 * t63 + t115 * t32) + t154 * (t112 * t32 - t63 * t115) + t133 * (t113 * t64 - t95 * t175) - t201 * (t32 * pkin(4) + t186)) * g(2) + (t125 * mrSges(3,1) + m(4) * t91 - t66 * mrSges(4,1) - m(5) * t185 - t34 * mrSges(5,1) - t202 * t96 + t135 * t65 - t155 * (t112 * t65 + t115 * t34) + t154 * (t112 * t34 - t65 * t115) + t133 * (t113 * t66 - t96 * t175) - t201 * (t34 * pkin(4) + t185)) * g(1) (-m(5) * t163 + t200 * t78 + t198 * t77 - t201 * (-pkin(12) * t182 - t77 * t187 + t163) - t155 * (t112 * t78 - t77 * t173) + t154 * (-t78 * t115 - t77 * t174) + t199 * t182) * g(3) + (-m(5) * t165 - t201 * (-pkin(12) * t184 - t53 * t187 + t165) - t155 * (t112 * t54 - t53 * t173) + t154 * (-t54 * t115 - t53 * t174) + t200 * t54 + t198 * t53 + t199 * t184) * g(2) + (-m(5) * t164 - t201 * (-pkin(12) * t183 - t57 * t187 + t164) + t154 * (-t58 * t115 - t57 * t174) + t200 * t58 + t198 * t57 + t199 * t183 - t155 * (t112 * t58 - t57 * t173)) * g(1) (-t201 * (t51 * pkin(4) + pkin(12) * t52) + t195 * t52 + t194 * t51) * g(3) + (-t201 * (t21 * pkin(4) + pkin(12) * t22) + t195 * t22 + t194 * t21) * g(2) + (-t201 * (-t25 * pkin(4) + pkin(12) * t26) + t195 * t26 - t194 * t25) * g(1) (t154 * (t112 * t77 + t115 * t52) + t155 * t15) * g(3) + (t155 * t1 + t154 * t206) * g(2) + (t154 * t6 + t155 * t5) * g(1) (-g(1) * t5 - g(2) * t1 - g(3) * t15) * m(7)];
taug  = t2(:);
