% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:44
% EndTime: 2019-12-31 19:49:47
% DurationCPUTime: 0.71s
% Computational Cost: add. (848->151), mult. (1684->200), div. (0->0), fcn. (889->6), ass. (0->83)
t171 = cos(pkin(8));
t175 = cos(qJ(2));
t170 = sin(pkin(8));
t173 = sin(qJ(2));
t212 = t170 * t173;
t219 = pkin(1) * qJD(2);
t147 = (t171 * t175 - t212) * t219;
t141 = qJD(1) * t147;
t228 = qJD(3) * qJD(4) + t141;
t172 = sin(qJ(4));
t168 = t172 ^ 2;
t174 = cos(qJ(4));
t169 = t174 ^ 2;
t210 = t172 * t174;
t227 = MDP(8) * t210 - (t168 - t169) * MDP(9);
t167 = qJD(1) + qJD(2);
t220 = pkin(1) * qJD(1);
t190 = t175 * t220;
t153 = pkin(2) * t167 + t190;
t191 = t173 * t220;
t156 = t171 * t191;
t134 = t170 * t153 + t156;
t130 = pkin(7) * t167 + t134;
t218 = t130 * t172;
t123 = qJD(3) * t174 - t218;
t226 = qJD(5) - t123;
t121 = -qJD(4) * pkin(4) + t226;
t124 = qJD(3) * t172 + t130 * t174;
t122 = qJD(4) * qJ(5) + t124;
t224 = t173 * MDP(5) + t175 * MDP(6);
t223 = (t168 + t169) * t167;
t193 = MDP(13) + MDP(15);
t222 = -MDP(14) + MDP(17);
t221 = pkin(2) * t171;
t164 = pkin(1) * t175 + pkin(2);
t211 = t171 * t173;
t205 = pkin(1) * t211 + t170 * t164;
t143 = pkin(7) + t205;
t176 = qJD(4) ^ 2;
t217 = t143 * t176;
t144 = t170 * t190 + t156;
t216 = t144 * t167;
t159 = pkin(2) * t170 + pkin(7);
t215 = t159 * t176;
t155 = t170 * t191;
t146 = t171 * t190 - t155;
t201 = qJD(4) * t172;
t207 = t146 * t201 + t174 * t216;
t206 = t228 * t174;
t203 = MDP(18) * t159;
t200 = qJD(4) * t174;
t199 = qJD(5) * t172;
t198 = t143 * MDP(18);
t197 = t176 * MDP(11);
t114 = (qJD(5) - t218) * qJD(4) + t206;
t116 = t130 * t200 + t228 * t172;
t188 = t114 * t174 + t116 * t172 + t121 * t200;
t186 = t123 + t218;
t180 = -pkin(4) * t174 - qJ(5) * t172 - pkin(3);
t183 = -pkin(1) * t212 + t164 * t171;
t131 = t180 - t183;
t185 = t131 * t167 - t147;
t184 = (-pkin(3) - t183) * t167 - t147;
t133 = t153 * t171 - t155;
t182 = pkin(4) * t172 - qJ(5) * t174;
t145 = (t170 * t175 + t211) * t219;
t181 = t145 * t167 + t217;
t140 = qJD(1) * t145;
t118 = t140 + (qJD(4) * t182 - t199) * t167;
t148 = pkin(4) * t201 - qJ(5) * t200 - t199;
t125 = t145 + t148;
t179 = -t125 * t167 - t118 - t217;
t178 = -t148 * t167 - t118 - t215;
t129 = -pkin(3) * t167 - t133;
t177 = t176 * t174 * MDP(10) + (t129 * t200 + t140 * t172) * MDP(14) + 0.2e1 * t227 * qJD(4) * t167;
t166 = t167 ^ 2;
t160 = -pkin(3) - t221;
t151 = t180 - t221;
t149 = t182 * t167;
t126 = t129 * t201;
t120 = t167 * t180 - t133;
t117 = t120 * t201;
t1 = [(-t133 * t145 + t134 * t147 - t140 * t183 + t141 * t205) * MDP(7) + t126 * MDP(13) + t117 * MDP(15) + (t223 * t147 + t188) * MDP(16) + (t118 * t131 + t120 * t125) * MDP(18) + (-t197 + t181 * MDP(14) + t179 * MDP(17) + (t116 * t143 + t121 * t147) * MDP(18) + (t184 * MDP(13) + t185 * MDP(15) + (-MDP(16) - t198) * t122) * qJD(4)) * t172 + ((-t140 - t181) * MDP(13) + t179 * MDP(15) + (t114 * t143 + t122 * t147) * MDP(18) + (t184 * MDP(14) + (-t120 - t185) * MDP(17) + t121 * t198) * qJD(4)) * t174 + t177 + t224 * t219 * (-qJD(1) - t167); (t133 * t144 - t134 * t146 + (-t140 * t171 + t141 * t170) * pkin(2)) * MDP(7) + (t126 + t207) * MDP(13) + (t117 + t207) * MDP(15) + (-t223 * t146 + t188) * MDP(16) + (t118 * t151 + (-t144 + t148) * t120) * MDP(18) + (-t197 + (t215 - t216) * MDP(14) + (t178 + t216) * MDP(17) + (t116 * t159 - t121 * t146) * MDP(18) + ((MDP(13) * t160 + MDP(15) * t151) * t167 + (-MDP(16) - t203) * t122) * qJD(4)) * t172 + ((-t140 - t215) * MDP(13) + t178 * MDP(15) + (t114 * t159 - t122 * t146) * MDP(18) + ((t160 * t167 + t146) * MDP(14) + (-t151 * t167 - t120 - t146) * MDP(17) + t121 * t203) * qJD(4)) * t174 + t177 + t224 * t220 * (-qJD(2) + t167); (t114 * t172 - t116 * t174 + (t121 * t172 + t122 * t174) * qJD(4)) * MDP(18) + (-t193 * t172 + t222 * t174) * t176; (qJ(5) * t114 - t120 * t149 - t121 * t124 + t226 * t122) * MDP(18) - t227 * t166 + (t186 * MDP(14) + (0.2e1 * qJD(5) - t186) * MDP(17) + t193 * t124) * qJD(4) + ((-t129 * MDP(13) - t120 * MDP(15) + t149 * MDP(17)) * t172 + (-t129 * MDP(14) + t149 * MDP(15) + t120 * MDP(17)) * t174) * t167 + t222 * t206 + (-MDP(18) * pkin(4) - t193) * t116; -t166 * MDP(15) * t210 + (-t166 * t168 - t176) * MDP(17) + (t120 * t167 * t172 - qJD(4) * t122 + t116) * MDP(18);];
tauc = t1;
