% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:36
% EndTime: 2019-12-31 20:57:44
% DurationCPUTime: 3.58s
% Computational Cost: add. (2221->314), mult. (5868->387), div. (0->0), fcn. (3602->4), ass. (0->144)
t211 = Ifges(6,4) + Ifges(5,5);
t216 = Ifges(5,4) + Ifges(4,5);
t221 = -Ifges(6,5) + t216;
t210 = Ifges(6,2) + Ifges(5,3);
t207 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t220 = -mrSges(5,1) - mrSges(4,1);
t219 = mrSges(5,2) + mrSges(4,3);
t138 = sin(qJ(3));
t139 = sin(qJ(2));
t140 = cos(qJ(2));
t195 = cos(qJ(3));
t114 = t138 * t140 + t139 * t195;
t101 = t114 * qJD(1);
t218 = t211 * t101;
t156 = t195 * t140;
t169 = qJD(1) * t139;
t100 = -qJD(1) * t156 + t138 * t169;
t217 = t211 * t100;
t215 = -Ifges(4,6) + Ifges(5,6);
t137 = qJD(2) + qJD(3);
t214 = t218 + (Ifges(5,6) - Ifges(6,6)) * t137 + t210 * t100;
t159 = pkin(2) * t140 + pkin(1);
t124 = t159 * qJD(1);
t213 = -qJ(4) * t101 - t124;
t96 = Ifges(4,4) * t100;
t212 = t207 * t101 + t221 * t137 + t217 - t96;
t184 = mrSges(6,3) * t100;
t80 = mrSges(6,2) * t137 + t184;
t85 = -mrSges(5,2) * t100 + mrSges(5,3) * t137;
t189 = t80 + t85;
t185 = mrSges(4,3) * t100;
t81 = -mrSges(4,2) * t137 - t185;
t188 = t81 + t85;
t178 = t101 * mrSges(4,3);
t186 = mrSges(5,2) * t101;
t209 = -t220 * t137 - t178 - t186;
t202 = -pkin(7) - pkin(6);
t125 = t202 * t139;
t118 = qJD(1) * t125;
t110 = qJD(2) * pkin(2) + t118;
t126 = t202 * t140;
t119 = qJD(1) * t126;
t157 = t195 * t119;
t73 = t138 * t110 - t157;
t170 = t138 * t119;
t72 = t195 * t110 + t170;
t160 = qJD(4) - t72;
t79 = t138 * t125 - t195 * t126;
t168 = qJD(1) * t140;
t173 = Ifges(3,6) * qJD(2);
t183 = Ifges(3,4) * t139;
t208 = t173 / 0.2e1 + (t140 * Ifges(3,2) + t183) * qJD(1) / 0.2e1 + pkin(6) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t168);
t164 = -Ifges(4,4) + t211;
t206 = -0.2e1 * pkin(1);
t141 = -pkin(3) - pkin(4);
t201 = -t100 / 0.2e1;
t200 = t100 / 0.2e1;
t198 = t101 / 0.2e1;
t197 = -t137 / 0.2e1;
t196 = t137 / 0.2e1;
t194 = pkin(2) * t138;
t193 = pkin(4) * t101;
t192 = pkin(6) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t169);
t158 = qJD(2) * t202;
t149 = qJD(1) * t158;
t111 = t139 * t149;
t147 = t140 * t149;
t15 = t73 * qJD(3) + t138 * t111 - t195 * t147;
t78 = -t195 * t125 - t126 * t138;
t190 = t15 * t78;
t64 = t101 * pkin(3) + t100 * qJ(4);
t77 = t137 * t114;
t63 = t77 * qJD(1);
t182 = qJ(4) * t63;
t43 = -pkin(3) * t137 + t160;
t179 = t100 * t43;
t177 = t101 * Ifges(4,4);
t136 = t137 * qJ(4);
t171 = qJ(5) * t100;
t37 = t171 + t73;
t27 = t136 + t37;
t176 = t101 * t27;
t130 = qJ(4) + t194;
t175 = t130 * t63;
t174 = Ifges(3,5) * qJD(2);
t75 = t195 * t118 + t170;
t167 = qJD(2) * t139;
t166 = qJD(3) * t138;
t165 = qJD(1) * qJD(2);
t163 = t195 * pkin(2);
t162 = pkin(2) * t169;
t145 = -t138 * t139 + t156;
t76 = t137 * t145;
t62 = t76 * qJD(1);
t155 = -t63 * mrSges(6,1) + t62 * mrSges(6,2);
t154 = qJD(3) * t195;
t153 = t139 * t165;
t152 = t140 * t165;
t74 = t118 * t138 - t157;
t150 = pkin(2) * t154;
t132 = -t163 - pkin(3);
t44 = t162 + t64;
t146 = qJ(4) * t114 + t159;
t14 = t110 * t154 + t195 * t111 + t119 * t166 + t138 * t147;
t120 = t139 * t158;
t121 = t140 * t158;
t18 = t195 * t120 + t138 * t121 + t125 * t154 + t126 * t166;
t12 = t137 * qJD(4) + t14;
t144 = -pkin(2) * t167 + qJ(4) * t76 + qJD(4) * t114;
t143 = -pkin(2) * t153 + qJ(4) * t62 + qJD(4) * t101;
t19 = t79 * qJD(3) + t138 * t120 - t195 * t121;
t93 = t101 * qJ(5);
t16 = t137 * t141 + t160 - t93;
t17 = t100 * t141 + qJD(5) - t213;
t3 = qJ(5) * t63 + qJD(5) * t100 + t12;
t4 = -t62 * qJ(5) - t101 * qJD(5) + t15;
t42 = pkin(3) * t100 + t213;
t50 = -t100 * Ifges(4,2) + t137 * Ifges(4,6) + t177;
t56 = t136 + t73;
t142 = -t42 * (mrSges(5,1) * t101 + mrSges(5,3) * t100) - t17 * (-mrSges(6,1) * t101 - mrSges(6,2) * t100) + t124 * (mrSges(4,1) * t101 - mrSges(4,2) * t100) + (-Ifges(6,5) * t100 + Ifges(6,6) * t101) * t196 + t50 * t198 - t16 * t184 - t72 * t185 + t56 * t186 + t12 * mrSges(5,3) - t14 * mrSges(4,2) + t3 * mrSges(6,2) - t4 * mrSges(6,1) + (t101 * t210 - t217) * t201 + (-t100 * t216 + t215 * t101) * t197 + t220 * t15 + (-Ifges(6,6) + t215) * t63 + t221 * t62 + (-Ifges(4,2) * t101 + t212 - t96) * t200 - (-t207 * t100 - t177 + t214 + t218) * t101 / 0.2e1;
t134 = Ifges(3,4) * t168;
t129 = -pkin(4) + t132;
t127 = t150 + qJD(4);
t99 = Ifges(3,1) * t169 + t134 + t174;
t82 = -mrSges(6,1) * t137 - t101 * mrSges(6,3);
t71 = -pkin(3) * t145 - t146;
t70 = mrSges(4,1) * t100 + mrSges(4,2) * t101;
t69 = -mrSges(6,1) * t100 + mrSges(6,2) * t101;
t68 = mrSges(5,1) * t100 - mrSges(5,3) * t101;
t55 = -qJ(5) * t145 + t79;
t54 = -qJ(5) * t114 + t78;
t41 = t93 + t75;
t40 = t74 + t171;
t39 = -t141 * t145 + t146;
t36 = t93 + t72;
t35 = -t64 - t193;
t26 = -t44 - t193;
t9 = pkin(3) * t77 - t144;
t8 = pkin(3) * t63 - t143;
t7 = -t76 * qJ(5) - t114 * qJD(5) + t19;
t6 = qJ(5) * t77 - qJD(5) * t145 + t18;
t5 = t141 * t77 + t144;
t1 = t141 * t63 + t143;
t2 = [m(4) * (t14 * t79 + t18 * t73 - t19 * t72 + t190) + m(5) * (t12 * t79 + t18 * t56 + t19 * t43 + t42 * t9 + t71 * t8 + t190) + (mrSges(6,2) * t1 - mrSges(5,3) * t8 - mrSges(6,3) * t4 + t219 * t15 + t164 * t63 + t207 * t62) * t114 - t209 * t19 + t188 * t18 + t39 * t155 + (t99 / 0.2e1 - t192 + t174 / 0.2e1 + (mrSges(3,2) * t206 + 0.3e1 / 0.2e1 * Ifges(3,4) * t140) * qJD(1)) * t140 * qJD(2) - (mrSges(5,1) * t8 - mrSges(6,1) * t1 - mrSges(5,2) * t12 - mrSges(4,3) * t14 + mrSges(6,3) * t3 + (Ifges(4,2) + t210) * t63 + t164 * t62) * t145 - t159 * (mrSges(4,1) * t63 + mrSges(4,2) * t62) + t6 * t80 + t7 * t82 + (-t54 * t62 + t55 * t63) * mrSges(6,3) + t9 * t68 + t5 * t69 + t71 * (mrSges(5,1) * t63 - mrSges(5,3) * t62) + (-t173 / 0.2e1 + (mrSges(3,1) * t206 - 0.3e1 / 0.2e1 * t183 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t140) * qJD(1) + (-0.2e1 * m(4) * t124 + t70 + qJD(1) * (-mrSges(4,1) * t145 + mrSges(4,2) * t114)) * pkin(2) - t208) * t167 + m(6) * (t1 * t39 + t16 * t7 + t17 * t5 + t27 * t6 + t3 * t55 + t4 * t54) + t219 * (t62 * t78 - t63 * t79) + (t207 * t198 + t211 * t200 - t72 * mrSges(4,3) + t43 * mrSges(5,2) - t16 * mrSges(6,3) + Ifges(6,5) * t197 + Ifges(4,4) * t201 - t124 * mrSges(4,2) - t42 * mrSges(5,3) + t17 * mrSges(6,2) + t216 * t196 + t212 / 0.2e1) * t76 + (t164 * t198 + t210 * t200 - t73 * mrSges(4,3) - t56 * mrSges(5,2) + t27 * mrSges(6,3) + Ifges(6,6) * t197 - Ifges(4,2) * t201 - t124 * mrSges(4,1) - t50 / 0.2e1 + t42 * mrSges(5,1) - t17 * mrSges(6,1) + t214 / 0.2e1 + t215 * t196) * t77; (m(5) * t43 + m(6) * t16 - t209 + t82) * pkin(2) * t166 + t208 * t169 + t209 * t74 - t188 * t75 + t189 * t127 - t70 * t162 - (Ifges(3,5) * t140 - Ifges(3,6) * t139) * t165 / 0.2e1 - Ifges(3,6) * t153 + (t12 * t130 + t132 * t15 - t42 * t44 - t43 * t74 + (t127 - t75) * t56) * m(5) + t142 + (t129 * t4 + t130 * t3 - t16 * t40 - t17 * t26 + (t127 - t41) * t27) * m(6) - (-Ifges(3,2) * t169 + t134 + t99) * t168 / 0.2e1 + (-t139 * (Ifges(3,1) * t140 - t183) / 0.2e1 + pkin(1) * (mrSges(3,1) * t139 + mrSges(3,2) * t140)) * qJD(1) ^ 2 + (t124 * t162 + t72 * t74 - t73 * t75 + (-t195 * t15 + t138 * t14 + (-t138 * t72 + t195 * t73) * qJD(3)) * pkin(2)) * m(4) + t168 * t192 + t73 * t178 - t41 * t80 - t40 * t82 - t44 * t68 - t26 * t69 + t81 * t150 + Ifges(3,5) * t152 + (t132 * t62 - t175 + t179) * mrSges(5,2) + (-t62 * t163 - t63 * t194) * mrSges(4,3) + (-t129 * t62 + t175 - t176) * mrSges(6,3) + (-mrSges(3,1) * t152 + mrSges(3,2) * t153) * pkin(6); t142 + (t178 + t209) * t73 - t188 * t72 + t189 * qJD(4) + (-t141 * t62 - t176 + t182) * mrSges(6,3) + (-pkin(3) * t62 + t179 - t182) * mrSges(5,2) - t36 * t80 - t37 * t82 - t64 * t68 - t35 * t69 + (t3 * qJ(4) + t4 * t141 - t16 * t37 - t17 * t35 + (-t36 + qJD(4)) * t27) * m(6) + (-t15 * pkin(3) + t12 * qJ(4) + t160 * t56 - t42 * t64 - t43 * t73) * m(5); (mrSges(5,2) - mrSges(6,3)) * t62 - t189 * t137 + (t68 - t69) * t101 + (-t17 * t101 - t27 * t137 + t4) * m(6) + (t101 * t42 - t137 * t56 + t15) * m(5); -t100 * t80 + t101 * t82 + 0.2e1 * (t1 / 0.2e1 + t27 * t201 + t16 * t198) * m(6) + t155;];
tauc = t2(:);
