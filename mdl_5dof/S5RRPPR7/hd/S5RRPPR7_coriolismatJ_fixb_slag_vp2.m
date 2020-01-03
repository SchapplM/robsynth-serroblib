% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:35:00
% EndTime: 2019-12-31 19:35:05
% DurationCPUTime: 2.10s
% Computational Cost: add. (4431->258), mult. (8734->350), div. (0->0), fcn. (9050->6), ass. (0->150)
t146 = sin(qJ(5));
t143 = t146 ^ 2;
t148 = cos(qJ(5));
t144 = t148 ^ 2;
t186 = t143 + t144;
t250 = mrSges(5,2) - mrSges(4,1);
t145 = sin(pkin(8));
t147 = sin(qJ(2));
t196 = cos(pkin(8));
t229 = cos(qJ(2));
t125 = t145 * t147 - t196 * t229;
t130 = (-qJ(3) - pkin(6)) * t147;
t184 = t229 * pkin(6);
t132 = t229 * qJ(3) + t184;
t246 = t145 * t130 + t196 * t132;
t248 = -t125 * pkin(4) + t246;
t127 = t145 * t229 + t196 * t147;
t142 = t147 * pkin(2);
t174 = qJ(4) * t125 + t142;
t237 = pkin(3) + pkin(7);
t52 = t237 * t127 + t174;
t32 = -t146 * t52 + t148 * t248;
t203 = t148 * t32;
t33 = t146 * t248 + t148 * t52;
t210 = t146 * t33;
t249 = -t210 - t203;
t225 = mrSges(4,3) + mrSges(5,1);
t176 = t144 / 0.2e1 + t143 / 0.2e1;
t247 = mrSges(6,3) * t176;
t221 = Ifges(6,4) * t148;
t133 = -Ifges(6,2) * t146 + t221;
t222 = Ifges(6,4) * t146;
t134 = Ifges(6,1) * t148 - t222;
t230 = t148 / 0.2e1;
t232 = t146 / 0.2e1;
t245 = t133 * t230 + t134 * t232;
t244 = -t186 * mrSges(6,3) / 0.2e1;
t217 = Ifges(6,6) * t148;
t220 = Ifges(6,5) * t146;
t243 = -Ifges(5,6) - Ifges(4,4) + t220 / 0.2e1 + t217 / 0.2e1;
t242 = 0.2e1 * t127;
t241 = m(5) / 0.2e1;
t240 = m(6) / 0.2e1;
t239 = m(4) * pkin(2);
t192 = t125 * t148;
t114 = Ifges(6,4) * t192;
t236 = -t114 / 0.4e1;
t235 = -t125 / 0.2e1;
t227 = pkin(2) * t145;
t136 = qJ(4) + t227;
t234 = t136 / 0.2e1;
t233 = -t146 / 0.2e1;
t231 = -t148 / 0.2e1;
t228 = m(5) * t127;
t219 = Ifges(6,2) * t148;
t218 = Ifges(6,6) * t127;
t140 = -t229 * pkin(2) - pkin(1);
t155 = -t127 * qJ(4) + t140;
t49 = t237 * t125 + t155;
t94 = t196 * t130 - t145 * t132;
t61 = -t127 * pkin(4) + t94;
t30 = -t146 * t49 - t148 * t61;
t31 = -t146 * t61 + t148 * t49;
t167 = t146 * t31 + t148 * t30;
t197 = t248 * t125;
t193 = t125 * t146;
t81 = t127 * mrSges(6,1) - mrSges(6,3) * t193;
t201 = t148 * t81;
t83 = -t127 * mrSges(6,2) + mrSges(6,3) * t192;
t207 = t146 * t83;
t205 = t148 * mrSges(6,1);
t211 = t146 * mrSges(6,2);
t129 = -t205 + t211;
t77 = t129 * t125;
t5 = (t225 * t125 - t77) * t125 + (t225 * t127 + t201 + t207) * t127 + m(6) * (t167 * t127 - t197) + (m(5) + m(4)) * (-t125 * t246 - t127 * t94);
t216 = qJD(1) * t5;
t190 = t127 * t148;
t208 = t146 * t81;
t79 = t125 * pkin(3) + t155;
t9 = t83 * t190 + t79 * t228 + (-mrSges(5,2) * t125 - t208 - m(6) * (t146 * t30 - t148 * t31) - mrSges(5,3) * t127) * t127;
t215 = qJD(1) * t9;
t117 = t125 * mrSges(5,3);
t118 = t125 * mrSges(4,2);
t171 = t219 + t222;
t45 = t171 * t125 + t218;
t202 = t148 * t45;
t213 = t127 * Ifges(6,5);
t47 = Ifges(6,1) * t193 + t114 + t213;
t209 = t146 * t47;
t46 = -Ifges(6,6) * t125 + t171 * t127;
t172 = Ifges(6,1) * t146 + t221;
t48 = -Ifges(6,5) * t125 + t172 * t127;
t78 = t129 * t127;
t80 = pkin(3) * t127 + t174;
t191 = t127 * t146;
t82 = -t125 * mrSges(6,1) - mrSges(6,3) * t191;
t84 = t125 * mrSges(6,2) + mrSges(6,3) * t190;
t1 = -t140 * t118 + t30 * t82 + t31 * t84 + t32 * t81 + t33 * t83 + t61 * t77 + t248 * t78 + m(6) * (t248 * t61 + t30 * t32 + t31 * t33) + (-t80 * mrSges(5,2) - t243 * t125 + t46 * t230 + t48 * t232) * t125 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t147 + (m(4) * t140 + mrSges(4,1) * t125) * pkin(2)) * t147 + (t209 / 0.2e1 + mrSges(4,2) * t142 + t202 / 0.2e1 - t80 * mrSges(5,3) + t140 * mrSges(4,1) + t243 * t127 + (Ifges(5,3) - Ifges(4,1) - Ifges(5,2) + Ifges(4,2) - Ifges(6,3)) * t125) * t127 + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) - Ifges(3,2)) * t147 + Ifges(3,4) * t229) * t229 + (m(5) * t80 - t127 * mrSges(5,2) + t117) * t79;
t214 = t1 * qJD(1);
t212 = t146 * mrSges(6,1);
t206 = t146 * t84;
t204 = t148 * mrSges(6,2);
t200 = t148 * t82;
t131 = t204 + t212;
t170 = Ifges(6,5) * t148 - Ifges(6,6) * t146;
t2 = -t131 * t197 + t31 * t81 - t30 * t83 + (t45 * t232 - t127 * t170 / 0.2e1 + (t134 * t233 + t219 * t232) * t125 + t167 * mrSges(6,3) + (t114 + t47) * t231) * t125;
t199 = t2 * qJD(1);
t179 = t196 * pkin(2);
t139 = -t179 - pkin(3);
t135 = -pkin(7) + t139;
t175 = t186 * t127;
t185 = t239 / 0.2e1;
t194 = t125 * t136;
t149 = -t194 * t241 + (t135 * t175 - t194) * t240 + t131 * t235 - t125 * t145 * t185 + (t139 * t241 - t196 * t185 + t244) * t127;
t150 = t80 * t241 + (-t146 * t32 + t148 * t33) * t240 + t82 * t233 + t84 * t230 + t147 * t185;
t6 = t250 * t127 - t117 + t118 + t149 - t150;
t198 = t6 * qJD(1);
t153 = (t204 / 0.2e1 + t212 / 0.2e1) * t127;
t159 = t208 / 0.2e1 + t83 * t231;
t12 = t125 * t247 + t153 + t159;
t195 = t12 * qJD(1);
t160 = t211 / 0.2e1 - t205 / 0.2e1;
t154 = t160 * t127;
t158 = -t207 / 0.2e1 - t201 / 0.2e1;
t14 = -t154 - t158;
t189 = t14 * qJD(1);
t25 = (-m(5) / 0.2e1 - t176 * m(6)) * t242;
t187 = t25 * qJD(1);
t183 = -Ifges(6,2) / 0.4e1 + Ifges(6,1) / 0.4e1;
t180 = -t248 * t129 / 0.2e1;
t169 = -t217 - t220;
t36 = t136 * t129 + t171 * t233 + t172 * t230 + t245;
t151 = t134 / 0.4e1 + mrSges(6,2) * t234 + t183 * t148;
t152 = -t133 / 0.4e1 + mrSges(6,1) * t234 - t183 * t146;
t162 = -t32 * mrSges(6,1) / 0.2e1 + t33 * mrSges(6,2) / 0.2e1;
t163 = t135 * t247;
t4 = t180 + (Ifges(6,3) / 0.2e1 - t163) * t125 + (-t135 * t81 / 0.2e1 - t47 / 0.4e1 + t236 - 0.3e1 / 0.4e1 * t213 + t152 * t125) * t146 + (t135 * t83 / 0.2e1 - t45 / 0.4e1 - 0.3e1 / 0.4e1 * t218 + (-0.3e1 / 0.4e1 * t222 + t151) * t125) * t148 + t162;
t165 = -t4 * qJD(1) + t36 * qJD(2);
t157 = -t200 / 0.2e1 - t206 / 0.2e1;
t10 = t160 * t125 + 0.2e1 * (t248 / 0.4e1 - t210 / 0.4e1 - t203 / 0.4e1) * m(6) + t157;
t85 = mrSges(5,3) + 0.4e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * t136 + t131;
t164 = qJD(1) * t10 + qJD(2) * t85;
t156 = t125 * t170;
t24 = t228 / 0.2e1 + t175 * t240 + (-m(6) * t186 / 0.4e1 - m(5) / 0.4e1) * t242;
t15 = -t154 + t158;
t13 = t244 * t125 + t153 - t159;
t8 = (-mrSges(5,1) + t160) * t125 - t157 + 0.2e1 * t246 * t241 + (-t249 + t248) * t240;
t7 = t149 + t150;
t3 = -t209 / 0.4e1 - t202 / 0.4e1 + t146 * t236 + t127 * t169 / 0.4e1 + t180 + Ifges(6,6) * t190 / 0.2e1 + Ifges(6,5) * t191 / 0.2e1 + Ifges(6,3) * t235 - t159 * t135 + (-t163 + t151 * t148 + (-0.3e1 / 0.4e1 * t221 + t152) * t146) * t125 - t162;
t11 = [qJD(2) * t1 + qJD(3) * t5 - qJD(4) * t9 - qJD(5) * t2, t7 * qJD(3) + t8 * qJD(4) + t3 * qJD(5) + t214 + (-t156 / 0.2e1 - mrSges(3,1) * t184 + t46 * t233 + t48 * t230 + Ifges(3,5) * t229 + t61 * t131 + (-t139 * mrSges(5,1) + mrSges(4,3) * t179 + Ifges(5,4) - Ifges(4,5)) * t125 + (-mrSges(4,3) * t227 + Ifges(5,5) - Ifges(4,6) + t245) * t127 + (pkin(6) * mrSges(3,2) - Ifges(3,6)) * t147 + (t145 * t239 - mrSges(4,2) + mrSges(5,3)) * t94 + (m(5) * t139 - t196 * t239 + t250) * t246 + (m(5) * t94 + m(6) * t61 - t127 * mrSges(5,1) + t78) * t136 + (-m(6) * t249 + t200 + t206) * t135 + t249 * mrSges(6,3)) * qJD(2), qJD(2) * t7 + qJD(4) * t24 + qJD(5) * t15 + t216, qJD(2) * t8 + qJD(3) * t24 + qJD(5) * t13 - t215, -t199 + t3 * qJD(2) + t15 * qJD(3) + t13 * qJD(4) + (-t31 * mrSges(6,1) - t30 * mrSges(6,2) + t156) * qJD(5); qJD(3) * t6 + qJD(4) * t10 + qJD(5) * t4 - t214, qJD(4) * t85 - qJD(5) * t36, t198, t164, (-t131 * t135 + t169) * qJD(5) - t165; -qJD(2) * t6 + qJD(4) * t25 - qJD(5) * t14 - t216, -t198, 0, t187, qJD(5) * t129 - t189; -qJD(2) * t10 - qJD(3) * t25 - qJD(5) * t12 + t215, -t164, -t187, 0, -qJD(5) * t131 - t195; -qJD(2) * t4 + qJD(3) * t14 + qJD(4) * t12 + t199, t165, t189, t195, 0;];
Cq = t11;
