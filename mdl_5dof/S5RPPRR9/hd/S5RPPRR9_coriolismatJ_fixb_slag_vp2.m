% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:20
% EndTime: 2019-12-31 18:02:23
% DurationCPUTime: 1.50s
% Computational Cost: add. (2941->249), mult. (5759->377), div. (0->0), fcn. (4819->6), ass. (0->157)
t149 = cos(qJ(4));
t230 = t149 / 0.2e1;
t146 = sin(qJ(5));
t139 = t146 ^ 2;
t148 = cos(qJ(5));
t141 = t148 ^ 2;
t250 = mrSges(6,3) * (t141 / 0.2e1 + t139 / 0.2e1);
t187 = t139 + t141;
t249 = t187 * mrSges(6,3);
t144 = cos(pkin(8));
t248 = t144 * qJ(2);
t138 = Ifges(6,4) * t148;
t115 = Ifges(6,1) * t146 + t138;
t143 = sin(pkin(8));
t150 = -pkin(1) - pkin(2);
t104 = t143 * t150 - pkin(6) + t248;
t225 = t149 * pkin(7);
t147 = sin(qJ(4));
t228 = t147 * pkin(4);
t120 = t225 - t228;
t196 = t146 * t147;
t39 = t104 * t196 + t148 * t120;
t194 = t147 * t148;
t40 = -t104 * t194 + t146 * t120;
t247 = -t39 * t146 + t40 * t148;
t176 = Ifges(6,2) * t146 - t138;
t140 = t147 ^ 2;
t142 = t149 ^ 2;
t190 = t149 * t143;
t94 = -t146 * t144 + t148 * t190;
t207 = t94 * t148;
t92 = -t148 * t144 - t146 * t190;
t209 = t92 * t146;
t174 = t207 - t209;
t185 = -0.1e1 + t187;
t22 = t174 * t149 + (-t185 * t140 - t142) * t143;
t240 = t22 / 0.2e1;
t202 = t143 * t147;
t26 = (-t174 + t190) * t202;
t105 = -t149 * mrSges(6,2) + mrSges(6,3) * t196;
t107 = t149 * mrSges(6,1) + mrSges(6,3) * t194;
t198 = t146 * t107;
t231 = t148 / 0.2e1;
t214 = t148 * mrSges(6,2);
t216 = t146 * mrSges(6,1);
t112 = t214 + t216;
t97 = t149 * t112;
t239 = t97 / 0.2e1;
t154 = t239 - t198 / 0.2e1 + t105 * t231;
t195 = t146 * t149;
t134 = t143 * qJ(2);
t179 = t144 * t150 - pkin(3) - t134;
t226 = t149 * pkin(4);
t227 = t147 * pkin(7);
t64 = -t179 + t226 + t227;
t34 = -t104 * t195 + t148 * t64;
t191 = t148 * t149;
t35 = t104 * t191 + t146 * t64;
t175 = -t146 * t34 + t148 * t35;
t156 = 0.2e1 * t104 * t149 - t175;
t95 = t146 * t143 + t144 * t191;
t206 = t95 * t148;
t192 = t148 * t143;
t93 = -t144 * t195 + t192;
t208 = t93 * t146;
t157 = (-t208 / 0.2e1 + t206 / 0.2e1) * mrSges(6,3);
t173 = t206 - t208;
t160 = t173 * pkin(7);
t106 = t147 * mrSges(6,2) + mrSges(6,3) * t195;
t108 = -t147 * mrSges(6,1) + mrSges(6,3) * t191;
t162 = -t92 * t108 / 0.2e1 - t94 * t106 / 0.2e1;
t178 = t92 * t39 + t94 * t40;
t199 = t144 * t149;
t182 = mrSges(5,2) * t199;
t111 = -t148 * mrSges(6,1) + t146 * mrSges(6,2);
t236 = t111 / 0.2e1;
t243 = m(6) * pkin(4);
t244 = m(6) / 0.2e1;
t245 = -m(6) / 0.2e1;
t6 = -t182 + t157 + t178 * t245 + t160 * t244 + ((t236 - mrSges(5,1) - t243 / 0.2e1) * t144 + (t156 * t245 + t154 + t239) * t143) * t147 + t162;
t246 = (t26 * qJD(2) + qJD(3) * t240) * m(6) - t6 * qJD(1);
t242 = -mrSges(6,1) / 0.2e1;
t241 = mrSges(6,2) / 0.2e1;
t238 = m(6) * (t173 - t199) * t147;
t237 = -t107 / 0.2e1;
t235 = -t112 / 0.2e1;
t234 = t146 / 0.2e1;
t233 = -t147 / 0.2e1;
t232 = -t148 / 0.2e1;
t229 = m(6) * t104;
t223 = m(6) * qJD(4);
t224 = t223 * t240;
t221 = Ifges(6,4) * t146;
t137 = Ifges(6,5) * t148;
t220 = Ifges(6,5) * t149;
t218 = Ifges(6,6) * t146;
t217 = Ifges(6,6) * t149;
t215 = t147 * mrSges(5,2);
t96 = t147 * t111;
t152 = (t143 * t96 / 0.2e1 + (t207 / 0.2e1 - t209 / 0.2e1) * mrSges(6,3)) * t147 + t92 * t105 / 0.2e1 + t94 * t237;
t167 = t95 * t241 + t93 * t242;
t9 = t152 + t167;
t210 = t9 * qJD(1);
t205 = -mrSges(5,1) + t111;
t161 = m(6) * t247;
t163 = t214 / 0.2e1 + t216 / 0.2e1;
t193 = t148 * t106;
t197 = t146 * t108;
t11 = (t140 - t142) * t229 / 0.2e1 + (t175 * t244 + t154) * t149 + (t161 / 0.2e1 + t193 / 0.2e1 - t197 / 0.2e1 - t163 * t147) * t147;
t204 = t11 * qJD(1);
t203 = t140 * t144;
t201 = t144 * t143;
t15 = t96 * t230 + (t105 * t234 + t107 * t231) * t147 - t140 * t250;
t189 = t15 * qJD(1);
t188 = Ifges(6,5) * t196 + Ifges(6,6) * t194;
t186 = qJD(4) * t149;
t184 = -t238 / 0.2e1;
t183 = t238 / 0.2e1;
t181 = t202 / 0.2e1;
t177 = t149 * mrSges(5,1) - t215;
t116 = Ifges(6,1) * t148 - t221;
t158 = t137 / 0.2e1 - t218 / 0.2e1 - Ifges(5,4);
t77 = t147 * t176 + t217;
t78 = -Ifges(6,6) * t147 + t149 * t176;
t131 = Ifges(6,4) * t196;
t79 = -Ifges(6,1) * t194 + t131 + t220;
t80 = -Ifges(6,5) * t147 - t116 * t149;
t1 = t39 * t107 + t34 * t108 + t40 * t105 + t35 * t106 + m(6) * (t34 * t39 + t35 * t40) + (t179 * mrSges(5,2) - t158 * t149 + t79 * t232 + t77 * t234) * t149 + (-t104 * t97 + t80 * t232 + t78 * t234 + t179 * mrSges(5,1) + t158 * t147 + (-Ifges(6,3) + Ifges(5,1) - Ifges(5,2) + (-t112 + t229) * t104) * t149) * t147;
t172 = t1 * qJD(1) + t11 * qJD(3);
t98 = Ifges(6,2) * t194 + t131;
t99 = t147 * t115;
t5 = t188 * t230 - t35 * t107 + t34 * t105 + (t104 * t96 + (-t99 / 0.2e1 + t77 / 0.2e1 + t35 * mrSges(6,3)) * t148 + (t79 / 0.2e1 + t98 / 0.2e1 - t34 * mrSges(6,3)) * t146) * t147;
t171 = t5 * qJD(1) - t15 * qJD(3);
t52 = t185 * t149 * t147;
t170 = -t22 * qJD(2) / 0.2e1 - t52 * qJD(3);
t168 = t40 * t241 + t39 * t242;
t166 = -t77 / 0.4e1 + t99 / 0.4e1 - pkin(7) * t105 / 0.2e1;
t165 = pkin(7) * t237 + t79 / 0.4e1 + t98 / 0.4e1;
t164 = -pkin(4) * t96 / 0.2e1 + t149 * t137 / 0.4e1;
t10 = mrSges(5,3) * t203 - t93 * t107 - t95 * t105 - m(6) * (t104 * t203 + t34 * t93 + t35 * t95) - m(3) * qJ(2) - mrSges(3,3) + (-m(4) * t134 + m(5) * t179 - mrSges(4,1) - t177) * t143 + (t142 * mrSges(5,3) + t140 * t112 - m(5) * (t140 + t142) * t104 - mrSges(4,2) - m(4) * t248) * t144;
t159 = t10 * qJD(1) + qJD(3) * t184;
t155 = t235 + t163;
t17 = t155 * t202;
t113 = Ifges(6,2) * t148 + t221;
t24 = -pkin(4) * t112 + (t115 / 0.2e1 - t176 / 0.2e1) * t148 + (t116 / 0.2e1 - t113 / 0.2e1) * t146;
t151 = t104 * t112 / 0.2e1 + (-t116 / 0.4e1 + t113 / 0.4e1) * t148 + (t115 / 0.4e1 - t176 / 0.4e1) * t146 + pkin(7) * t250;
t3 = (t220 / 0.2e1 + t165) * t148 + (-0.3e1 / 0.4e1 * t217 + t166) * t146 + (Ifges(6,3) / 0.2e1 + t151) * t147 + t164 + t168;
t45 = t155 * t149;
t153 = t3 * qJD(1) - t17 * qJD(2) + t45 * qJD(3) + t24 * qJD(4);
t119 = t140 * t201;
t46 = (-t163 + t235) * t149;
t18 = t112 * t181 + t163 * t202;
t8 = t152 - t167;
t7 = -t97 * t202 + t105 * t192 * t233 + t181 * t198 - t182 / 0.2e1 + t157 - t162 + (t156 * t202 + t160 + t178) * t244 + (mrSges(5,2) * t230 + (-pkin(4) * t244 + t236) * t147) * t144;
t4 = -Ifges(6,5) * t191 / 0.2e1 + Ifges(6,6) * t195 / 0.2e1 + Ifges(6,3) * t233 + (-t217 / 0.4e1 + t166) * t146 + t165 * t148 + t151 * t147 + t164 - t168;
t2 = qJD(2) * t183 + t11 * qJD(4) - t15 * qJD(5);
t12 = [-qJD(2) * t10 + qJD(4) * t1 + qJD(5) * t5, 0.2e1 * ((t92 * t93 + t94 * t95 + t119) * t244 + m(5) * (t119 + (t142 - 0.1e1) * t201) / 0.2e1) * qJD(2) + t7 * qJD(4) + t8 * qJD(5) - t159, t2, t7 * qJD(2) + t4 * qJD(5) + (t115 * t232 + t113 * t234 - Ifges(5,5) + (t205 - t243) * t104) * t186 + t172 + (t104 * t215 + pkin(4) * t97 + (Ifges(6,5) * t146 + Ifges(6,6) * t148) * t233 + t80 * t234 + t78 * t231 + Ifges(5,6) * t147 + (t161 + t193 - t197) * pkin(7) + t247 * mrSges(6,3)) * qJD(4), t8 * qJD(2) + t4 * qJD(4) + (-t35 * mrSges(6,1) - t34 * mrSges(6,2) + t188) * qJD(5) + t171; -t6 * qJD(4) + t9 * qJD(5) + t159, t26 * t223, qJD(1) * t184 + t224, t18 * qJD(5) + (t149 * t111 + m(6) * (-t187 * t227 - t226) - t147 * t249 - t177) * qJD(4) * t143 + t246, t210 + t18 * qJD(4) + (-t94 * mrSges(6,1) - t92 * mrSges(6,2)) * qJD(5); t2, qJD(1) * t183 + t224, t52 * t223, t204 + t46 * qJD(5) + t205 * qJD(4) * t147 + (-mrSges(5,2) + t249) * t186 + ((t187 * t225 - t228) * qJD(4) - t170) * m(6), t46 * qJD(4) + t96 * qJD(5) - t189; qJD(2) * t6 + qJD(5) * t3 - t172, -t17 * qJD(5) - t246, m(6) * t170 + t45 * qJD(5) - t204, t24 * qJD(5), (pkin(7) * t111 + t137 - t218) * qJD(5) + t153; -qJD(2) * t9 - qJD(4) * t3 - t171, t17 * qJD(4) - t210, -t45 * qJD(4) + t189, -t153, 0;];
Cq = t12;
