% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:29
% EndTime: 2019-12-31 19:12:36
% DurationCPUTime: 2.88s
% Computational Cost: add. (3801->337), mult. (8198->484), div. (0->0), fcn. (5070->6), ass. (0->175)
t138 = sin(qJ(5));
t140 = cos(qJ(5));
t162 = mrSges(6,1) * t138 + mrSges(6,2) * t140;
t135 = qJD(3) + qJD(4);
t233 = sin(qJ(4));
t139 = sin(qJ(3));
t142 = -pkin(1) - pkin(6);
t118 = qJD(1) * t142 + qJD(2);
t171 = pkin(7) * qJD(1) - t118;
t95 = t171 * t139;
t185 = t233 * t95;
t234 = cos(qJ(4));
t141 = cos(qJ(3));
t109 = t141 * t118;
t195 = qJD(1) * t141;
t96 = -pkin(7) * t195 + t109;
t93 = qJD(3) * pkin(3) + t96;
t58 = t234 * t93 + t185;
t53 = -t135 * pkin(4) - t58;
t152 = t53 * t162;
t157 = Ifges(6,5) * t140 - Ifges(6,6) * t138;
t217 = Ifges(6,4) * t140;
t159 = -Ifges(6,2) * t138 + t217;
t218 = Ifges(6,4) * t138;
t161 = Ifges(6,1) * t140 - t218;
t235 = t140 / 0.2e1;
t107 = t139 * t234 + t141 * t233;
t102 = t107 * qJD(1);
t98 = qJD(5) + t102;
t240 = t98 / 0.2e1;
t177 = t233 * t139;
t168 = qJD(1) * t177;
t178 = t234 * t141;
t103 = qJD(1) * t178 - t168;
t88 = t103 * t140 + t135 * t138;
t242 = t88 / 0.2e1;
t87 = -t103 * t138 + t135 * t140;
t244 = t87 / 0.2e1;
t232 = Ifges(6,4) * t88;
t35 = t87 * Ifges(6,2) + t98 * Ifges(6,6) + t232;
t256 = -t35 / 0.2e1;
t86 = Ifges(6,4) * t87;
t36 = t88 * Ifges(6,1) + t98 * Ifges(6,5) + t86;
t257 = t138 * t256 + t157 * t240 + t159 * t244 + t161 * t242 + t36 * t235 + t152;
t214 = t103 * mrSges(5,3);
t222 = -mrSges(5,1) * t135 - mrSges(6,1) * t87 + mrSges(6,2) * t88 + t214;
t255 = (qJ(2) * (m(3) + m(4)));
t190 = qJD(1) * qJD(3);
t172 = t141 * t190;
t173 = t234 * qJD(4);
t254 = t234 * qJD(3) + t173;
t175 = qJD(4) * t233;
t253 = -qJD(3) * t233 - t175;
t191 = qJD(5) * t140;
t192 = qJD(5) * t138;
t55 = -mrSges(6,2) * t98 + mrSges(6,3) * t87;
t56 = mrSges(6,1) * t98 - mrSges(6,3) * t88;
t252 = -t56 * t191 - t55 * t192;
t186 = t234 * t95;
t59 = t233 * t93 - t186;
t54 = t135 * pkin(8) + t59;
t196 = qJD(1) * t139;
t115 = pkin(3) * t196 + qJD(1) * qJ(2);
t60 = pkin(4) * t102 - pkin(8) * t103 + t115;
t14 = -t138 * t54 + t140 * t60;
t15 = t138 * t60 + t140 * t54;
t251 = -t138 * t14 + t140 * t15;
t75 = t135 * t102;
t42 = qJD(5) * t87 - t140 * t75;
t147 = t254 * t141;
t76 = qJD(1) * t147 - t135 * t168;
t11 = mrSges(6,1) * t76 - mrSges(6,3) * t42;
t43 = -qJD(5) * t88 + t138 * t75;
t12 = -mrSges(6,2) * t76 + mrSges(6,3) * t43;
t250 = -t138 * t11 + t140 * t12;
t194 = qJD(3) * t139;
t153 = t171 * t194;
t193 = qJD(3) * t141;
t94 = t171 * t193;
t17 = t58 * qJD(4) + t153 * t233 - t234 * t94;
t110 = pkin(3) * t172 + qJD(1) * qJD(2);
t26 = pkin(4) * t76 + pkin(8) * t75 + t110;
t2 = qJD(5) * t14 + t138 * t26 + t140 * t17;
t3 = -qJD(5) * t15 - t138 * t17 + t140 * t26;
t249 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t42 + Ifges(6,6) * t43;
t221 = mrSges(5,3) * t102;
t90 = -mrSges(5,2) * t135 - t221;
t248 = -t138 * t56 + t140 * t55 + t90;
t231 = t138 * t3;
t247 = -t14 * t191 - t15 * t192 - t231;
t246 = qJD(1) ^ 2;
t245 = -t87 / 0.2e1;
t243 = -t88 / 0.2e1;
t241 = -t98 / 0.2e1;
t236 = t138 / 0.2e1;
t230 = t140 * t2;
t18 = t59 * qJD(4) - t153 * t234 - t233 * t94;
t223 = pkin(7) - t142;
t111 = t223 * t139;
t112 = t223 * t141;
t84 = -t111 * t233 + t112 * t234;
t229 = t18 * t84;
t228 = t75 * mrSges(5,3);
t227 = t76 * mrSges(5,3);
t226 = t87 * Ifges(6,6);
t225 = t88 * Ifges(6,5);
t224 = t98 * Ifges(6,3);
t220 = Ifges(4,4) * t139;
t219 = Ifges(4,4) * t141;
t97 = Ifges(5,4) * t102;
t215 = t102 * Ifges(5,2);
t213 = t103 * Ifges(5,1);
t212 = t103 * Ifges(5,4);
t106 = t177 - t178;
t211 = t106 * t18;
t210 = t135 * Ifges(5,5);
t209 = t135 * Ifges(5,6);
t202 = Ifges(4,5) * qJD(3);
t201 = Ifges(4,6) * qJD(3);
t200 = t102 * t138;
t199 = t102 * t140;
t114 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t195;
t198 = t139 * t114;
t128 = t139 * pkin(3) + qJ(2);
t119 = pkin(3) * t193 + qJD(2);
t189 = t234 * pkin(3);
t188 = t233 * pkin(3);
t187 = pkin(3) * t195;
t166 = t223 * t194;
t79 = pkin(4) * t103 + pkin(8) * t102;
t82 = -t254 * t139 + t253 * t141;
t83 = t253 * t139 + t147;
t165 = t58 * t82 + t59 * t83;
t164 = mrSges(4,1) * t139 + mrSges(4,2) * t141;
t163 = -mrSges(6,1) * t140 + mrSges(6,2) * t138;
t160 = Ifges(6,1) * t138 + t217;
t158 = Ifges(6,2) * t140 + t218;
t156 = Ifges(6,5) * t138 + Ifges(6,6) * t140;
t155 = t138 * t15 + t14 * t140;
t154 = -t138 * t55 - t140 * t56;
t80 = pkin(4) * t107 + pkin(8) * t106 + t128;
t85 = -t111 * t234 - t112 * t233;
t38 = -t138 * t85 + t140 * t80;
t39 = t138 * t80 + t140 * t85;
t151 = qJ(2) * (mrSges(4,1) * t141 - mrSges(4,2) * t139);
t146 = -qJD(5) * t155 - t231;
t145 = m(6) * (t230 + t247) + t250;
t34 = t224 + t225 + t226;
t67 = t209 + t212 - t215;
t68 = t210 - t97 + t213;
t8 = t42 * Ifges(6,4) + t43 * Ifges(6,2) + t76 * Ifges(6,6);
t9 = t42 * Ifges(6,1) + t43 * Ifges(6,4) + t76 * Ifges(6,5);
t144 = t103 * t67 / 0.2e1 - Ifges(5,5) * t75 - t17 * mrSges(5,2) + t42 * t160 / 0.2e1 + t43 * t158 / 0.2e1 + t102 * t152 + t200 * t256 - t135 * (-Ifges(5,5) * t102 - Ifges(5,6) * t103) / 0.2e1 - t115 * (mrSges(5,1) * t103 - mrSges(5,2) * t102) + (Ifges(6,3) * t103 - t102 * t157) * t241 + (Ifges(6,5) * t103 - t102 * t161) * t243 + (Ifges(6,6) * t103 - t102 * t159) * t245 - t14 * (mrSges(6,1) * t103 + mrSges(6,3) * t199) - t15 * (-mrSges(6,2) * t103 + mrSges(6,3) * t200) - t58 * t221 + mrSges(6,3) * t230 + t8 * t235 + t9 * t236 + t36 * t199 / 0.2e1 + (-Ifges(5,6) + t156 / 0.2e1) * t76 + (-Ifges(5,2) * t103 + t68 - t97) * t102 / 0.2e1 - (-Ifges(5,1) * t102 - t212 + t34) * t103 / 0.2e1 + (-mrSges(5,1) + t163) * t18 + t257 * qJD(5);
t131 = -t189 - pkin(4);
t113 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t196;
t108 = qJD(1) * t164;
t105 = qJD(3) * t112;
t101 = t202 + (Ifges(4,1) * t141 - t220) * qJD(1);
t100 = t201 + (-Ifges(4,2) * t139 + t219) * qJD(1);
t78 = mrSges(5,1) * t102 + mrSges(5,2) * t103;
t72 = Ifges(6,3) * t76;
t66 = t79 + t187;
t62 = t234 * t96 + t185;
t61 = t233 * t96 - t186;
t45 = qJD(4) * t85 - t105 * t233 - t166 * t234;
t44 = -qJD(4) * t84 - t105 * t234 + t166 * t233;
t37 = pkin(4) * t83 - pkin(8) * t82 + t119;
t25 = t138 * t79 + t140 * t58;
t24 = -t138 * t58 + t140 * t79;
t23 = t138 * t66 + t140 * t62;
t22 = -t138 * t62 + t140 * t66;
t10 = -mrSges(6,1) * t43 + mrSges(6,2) * t42;
t5 = -qJD(5) * t39 - t138 * t44 + t140 * t37;
t4 = qJD(5) * t38 + t138 * t37 + t140 * t44;
t1 = [t128 * (mrSges(5,1) * t76 - mrSges(5,2) * t75) + qJD(2) * t108 + t119 * t78 + t44 * t90 + t84 * t10 + t4 * t55 + t5 * t56 + t38 * t11 + t39 * t12 + (-t209 / 0.2e1 + t115 * mrSges(5,1) + t215 / 0.2e1 - t212 / 0.2e1 + t34 / 0.2e1 - t67 / 0.2e1 + t224 / 0.2e1 + t226 / 0.2e1 + t225 / 0.2e1 + t14 * mrSges(6,1) - t15 * mrSges(6,2)) * t83 + (t210 / 0.2e1 + t115 * mrSges(5,2) - t97 / 0.2e1 + t213 / 0.2e1 + t68 / 0.2e1 - t155 * mrSges(6,3) + t257) * t82 + t222 * t45 + m(5) * (t110 * t128 + t115 * t119 + t17 * t85 + t44 * t59 - t45 * t58 + t229) + m(6) * (t14 * t5 + t15 * t4 + t2 * t39 + t3 * t38 + t45 * t53 + t229) + (-t75 * t84 - t76 * t85 - t165) * mrSges(5,3) + ((t142 * t113 - t100 / 0.2e1 - t201 / 0.2e1) * t141 + (-t142 * t114 - t101 / 0.2e1 - t202 / 0.2e1) * t139) * qJD(3) + (-t17 * mrSges(5,3) + Ifges(5,4) * t75 + t110 * mrSges(5,1) + t72 / 0.2e1 + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t76 + t249) * t107 + (-t42 * t161 / 0.2e1 - t43 * t159 / 0.2e1 - t140 * t9 / 0.2e1 + t8 * t236 - t110 * mrSges(5,2) + Ifges(5,1) * t75 + (-mrSges(5,3) - t162) * t18 + (t138 * t2 + t140 * t3) * mrSges(6,3) + (t251 * mrSges(6,3) + t156 * t240 + t158 * t244 + t160 * t242 + t53 * t163 + t35 * t235 + t36 * t236) * qJD(5) + (-t157 / 0.2e1 + Ifges(5,4)) * t76) * t106 + (((2 * mrSges(3,3)) + t164 + (2 * t255)) * qJD(2) + (0.2e1 * t151 + (0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,1)) * t139 * t141 + (0.3e1 / 0.2e1 * t139 ^ 2 - 0.3e1 / 0.2e1 * t141 ^ 2) * Ifges(4,4)) * qJD(3)) * qJD(1); -t222 * t82 + (t10 - t228) * t106 + (t141 * t113 - t198) * qJD(3) + m(6) * (-t53 * t82 + t211) + m(5) * (t165 + t211) + (-mrSges(3,3) - t255) * t246 + (m(5) * t17 + qJD(5) * t154 + t145 - t227) * t107 + (-m(5) * t115 - m(6) * t155 - t108 + t154 - t78) * qJD(1) + (m(6) * t251 + t248) * t83; t131 * t10 - t62 * t90 - t23 * t55 - t22 * t56 - t188 * t227 + t189 * t228 - t190 * Ifges(4,5) * t139 / 0.2e1 + t101 * t196 / 0.2e1 - t113 * t109 + t100 * t195 / 0.2e1 - t78 * t187 + t144 + t59 * t214 - Ifges(4,6) * t172 / 0.2e1 - t222 * t61 + (t18 * t131 - t14 * t22 - t15 * t23 - t53 * t61) * m(6) + (-t115 * t187 + t58 * t61 - t59 * t62) * m(5) + (-t141 * (-Ifges(4,1) * t139 - t219) / 0.2e1 + t139 * (-Ifges(4,2) * t141 - t220) / 0.2e1 - t151) * t246 + (-mrSges(4,1) * t194 - mrSges(4,2) * t193 + t198) * t118 + t247 * mrSges(6,3) + (m(6) * (t146 + t230) + t250 + t252) * (t188 + pkin(8)) + (t222 * t175 + (t233 * t53 + t251 * t234) * qJD(4) * m(6) + m(5) * (t233 * t17 - t234 * t18 + (-t233 * t58 + t234 * t59) * qJD(4)) + t248 * t173) * pkin(3); -t58 * t90 - t25 * t55 - t24 * t56 - pkin(4) * t10 + t146 * mrSges(6,3) + (t214 - t222) * t59 + t144 + (t145 + t252) * pkin(8) + (-pkin(4) * t18 - t14 * t24 - t15 * t25 - t53 * t59) * m(6); t72 - t53 * (mrSges(6,1) * t88 + mrSges(6,2) * t87) + (Ifges(6,1) * t87 - t232) * t243 + t35 * t242 + (Ifges(6,5) * t87 - Ifges(6,6) * t88) * t241 - t14 * t55 + t15 * t56 + (t14 * t87 + t15 * t88) * mrSges(6,3) + (-Ifges(6,2) * t88 + t36 + t86) * t245 + t249;];
tauc = t1(:);
