% Calculate Coriolis joint torque vector for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:24
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:24:05
% EndTime: 2021-01-15 11:24:08
% DurationCPUTime: 0.91s
% Computational Cost: add. (1050->182), mult. (2157->233), div. (0->0), fcn. (1235->4), ass. (0->80)
t199 = sin(qJ(3));
t200 = cos(qJ(3));
t251 = -t200 * t199 * MDP(7) + MDP(8) * (t199 ^ 2 - t200 ^ 2) + qJ(2) * (t200 * MDP(12) - t199 * MDP(13));
t250 = MDP(16) + MDP(19);
t247 = MDP(12) * t199 + MDP(13) * t200;
t197 = sin(pkin(7));
t198 = cos(pkin(7));
t173 = t197 * t200 + t198 * t199;
t246 = t173 * qJD(1);
t245 = MDP(14) + MDP(18);
t244 = -MDP(15) + MDP(20);
t201 = -pkin(1) - pkin(6);
t236 = qJ(4) - t201;
t218 = t236 * t200;
t225 = qJD(4) * t199;
t160 = -qJD(3) * t218 - t225;
t224 = qJD(4) * t200;
t227 = qJD(3) * t199;
t207 = t227 * t236 - t224;
t138 = t160 * t197 - t198 * t207;
t139 = t198 * t160 + t197 * t207;
t177 = t236 * t199;
t150 = -t177 * t197 + t198 * t218;
t151 = -t198 * t177 - t197 * t218;
t226 = qJD(3) * t200;
t219 = qJD(1) * t226;
t178 = t198 * t219;
t231 = qJD(1) * t199;
t220 = t197 * t231;
t158 = qJD(3) * t220 - t178;
t159 = qJD(3) * t246;
t230 = qJD(1) * t200;
t169 = t198 * t230 - t220;
t242 = t138 * t169 - t139 * t246 - t150 * t159 + t151 * t158;
t164 = t169 ^ 2;
t179 = qJD(1) * t201 + qJD(2);
t152 = t179 * t226 + (-qJ(4) * t226 - t225) * qJD(1);
t204 = -t179 * t227 + (qJ(4) * t227 - t224) * qJD(1);
t128 = t152 * t197 - t198 * t204;
t241 = t128 * t150;
t172 = -t197 * t199 + t198 * t200;
t240 = t128 * t172;
t161 = (-qJ(4) * qJD(1) + t179) * t199;
t237 = t161 * t197;
t155 = t198 * t161;
t129 = t198 * t152 + t197 * t204;
t162 = -qJ(4) * t230 + t200 * t179;
t157 = qJD(3) * pkin(3) + t162;
t136 = t197 * t157 + t155;
t175 = pkin(3) * t219 + qJD(1) * qJD(2);
t188 = t199 * pkin(3) + qJ(2);
t229 = qJD(3) * t138;
t228 = qJD(3) * t139;
t141 = t162 * t198 - t237;
t223 = qJD(5) - t141;
t180 = pkin(3) * t226 + qJD(2);
t222 = qJD(3) * qJD(5);
t176 = pkin(3) * t231 + qJD(1) * qJ(2) + qJD(4);
t217 = -qJ(2) * MDP(6) - MDP(5);
t215 = qJD(3) * t141 - t129;
t135 = t157 * t198 - t237;
t209 = -pkin(4) * t158 + qJ(5) * t159 + t175;
t127 = t222 + t129;
t132 = -qJD(3) * pkin(4) + qJD(5) - t135;
t133 = qJD(3) * qJ(5) + t136;
t167 = t197 * t227 - t198 * t226;
t168 = -t197 * t226 - t198 * t227;
t206 = t127 * t173 - t132 * t168 - t133 * t167 - t240;
t205 = t129 * t173 + t135 * t168 - t136 * t167 - t240;
t203 = qJD(1) ^ 2;
t202 = qJD(3) ^ 2;
t189 = -pkin(3) * t198 - pkin(4);
t186 = pkin(3) * t197 + qJ(5);
t144 = pkin(4) * t173 - qJ(5) * t172 + t188;
t142 = pkin(3) * t230 + pkin(4) * t169 + qJ(5) * t246;
t140 = t162 * t197 + t155;
t137 = pkin(4) * t246 - qJ(5) * t169 + t176;
t131 = -pkin(4) * t167 - qJ(5) * t168 - qJD(5) * t172 + t180;
t126 = -qJD(5) * t169 + t209;
t1 = [(-t158 * t188 - t167 * t176 + t173 * t175 + t180 * t246 - t229) * MDP(14) + (-t159 * t188 + t168 * t176 + t169 * t180 + t172 * t175 - t228) * MDP(15) + (-t205 + t242) * MDP(16) + (t129 * t151 - t135 * t138 + t136 * t139 + t175 * t188 + t176 * t180 + t241) * MDP(17) + (t126 * t173 + t131 * t246 - t137 * t167 - t144 * t158 - t229) * MDP(18) + (-t206 + t242) * MDP(19) + (-t126 * t172 - t131 * t169 - t137 * t168 + t144 * t159 + t228) * MDP(20) + (t126 * t144 + t127 * t151 + t131 * t137 + t132 * t138 + t133 * t139 + t241) * MDP(21) + ((-MDP(13) * t201 - MDP(10)) * t200 + (-MDP(12) * t201 - MDP(9)) * t199) * t202 + (0.2e1 * (-t217 + t247) * qJD(2) + 0.2e1 * t251 * qJD(3)) * qJD(1); (-qJD(1) * t176 + t205) * MDP(17) + (-qJD(1) * t137 + t206) * MDP(21) + t217 * t203 + t247 * (-t202 - t203) + t245 * (-qJD(1) * t246 + qJD(3) * t168) + t244 * (qJD(1) * t169 - qJD(3) * t167) + t250 * (t173 * t158 + t159 * t172 + t167 * t246 - t168 * t169); t215 * MDP(15) + (t135 * t140 - t136 * t141) * MDP(17) + (t158 * t186 - t159 * t189) * MDP(19) + (-t215 + 0.2e1 * t222) * MDP(20) + (t127 * t186 + t128 * t189 - t132 * t140 + t133 * t223 - t137 * t142) * MDP(21) + (-t176 * MDP(14) + (t136 - t140) * MDP(16) - t137 * MDP(18) + (t133 - t140) * MDP(19) + t142 * MDP(20)) * t169 + (t176 * MDP(15) + (-t135 + t141) * MDP(16) - t142 * MDP(18) + (t132 - t223) * MDP(19) - t137 * MDP(20)) * t246 + ((t158 * t197 + t159 * t198) * MDP(16) + (-t128 * t198 + t129 * t197) * MDP(17) + (-MDP(14) * t246 - t169 * MDP(15) - t176 * MDP(17)) * t230) * pkin(3) - t251 * t203 + t245 * (qJD(3) * t140 - t128); (t135 * t169 + t136 * t246 + t175) * MDP(17) + (t133 * t246 + (-qJD(5) - t132) * t169 + t209) * MDP(21) + t245 * t178 + (0.2e1 * t244 * t246 + t245 * (t169 - t220)) * qJD(3) + t250 * (-t246 ^ 2 - t164); t169 * t246 * MDP(18) + (-t164 - t202) * MDP(20) + (t137 * t169 + t128) * MDP(21) + ((-t197 * t230 - t198 * t231 + t246) * MDP(19) - t133 * MDP(21)) * qJD(3);];
tauc = t1;
