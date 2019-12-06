% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:22
% EndTime: 2019-12-05 16:46:25
% DurationCPUTime: 0.90s
% Computational Cost: add. (897->165), mult. (1743->215), div. (0->0), fcn. (1090->6), ass. (0->92)
t236 = 2 * qJD(4);
t180 = sin(qJ(4));
t178 = t180 ^ 2;
t183 = cos(qJ(4));
t179 = t183 ^ 2;
t237 = t178 + t179;
t246 = -t237 * MDP(16) + MDP(7);
t181 = sin(qJ(3));
t182 = sin(qJ(2));
t217 = qJD(1) * t182;
t171 = t181 * t217;
t184 = cos(qJ(3));
t185 = cos(qJ(2));
t209 = t185 * qJD(1);
t158 = t184 * t209 - t171;
t215 = qJD(3) * t184;
t204 = pkin(2) * t215;
t238 = t158 - t204;
t208 = qJD(1) * qJD(2);
t245 = qJD(3) * t217 + t182 * t208;
t223 = t180 * t183;
t244 = MDP(8) * t223 - (t178 - t179) * MDP(9);
t162 = t181 * t182 - t184 * t185;
t177 = qJD(2) + qJD(3);
t242 = t162 * t177;
t170 = qJD(2) * pkin(2) + t209;
t199 = t185 * t208;
t219 = t245 * t181;
t135 = (qJD(3) * t170 + t199) * t184 - t219;
t132 = t183 * t135;
t154 = t170 * t181 + t184 * t217;
t147 = pkin(7) * t177 + t154;
t224 = t180 * t147;
t126 = t132 + (qJD(5) - t224) * qJD(4);
t131 = t180 * t135;
t212 = qJD(4) * t183;
t128 = t147 * t212 + t131;
t240 = t126 * t183 + t128 * t180;
t210 = t180 * qJD(5);
t213 = qJD(4) * t180;
t155 = pkin(4) * t213 - qJ(5) * t212 - t210;
t216 = qJD(3) * t181;
t203 = pkin(2) * t216;
t149 = t155 + t203;
t163 = t181 * t185 + t182 * t184;
t157 = t163 * qJD(1);
t239 = t149 - t157;
t186 = qJD(4) ^ 2;
t235 = pkin(7) * t186;
t234 = t177 * pkin(3);
t233 = t184 * pkin(2);
t232 = MDP(18) * pkin(7);
t231 = qJD(4) * pkin(4);
t141 = t177 * t163;
t230 = t141 * t177;
t153 = t184 * t170 - t171;
t229 = t153 * t177;
t228 = t154 * t177;
t164 = -pkin(4) * t183 - qJ(5) * t180 - pkin(3);
t227 = t164 * t177;
t173 = pkin(2) * t181 + pkin(7);
t226 = t173 * t186;
t225 = t177 * t183;
t222 = t183 * t147;
t221 = t153 * t213 + t154 * t225;
t220 = t157 * t225 + t158 * t213;
t214 = qJD(4) * qJ(5);
t211 = t173 * MDP(18);
t207 = MDP(13) + MDP(15);
t206 = MDP(14) - MDP(17);
t197 = qJD(5) + t224;
t138 = t197 - t231;
t201 = t138 * t212 + t240;
t136 = t170 * t216 + t181 * t199 + t184 * t245;
t196 = -t157 + t203;
t195 = pkin(4) * t180 - qJ(5) * t183;
t194 = t163 * t186 + t230;
t193 = t242 * t236;
t129 = (qJD(4) * t195 - t210) * t177 + t136;
t192 = -t155 * t177 - t129 - t235;
t159 = t164 - t233;
t191 = t159 * t177 - t204;
t190 = (-pkin(3) - t233) * t177 - t204;
t146 = -t153 - t234;
t188 = t186 * t183 * MDP(10) + (t136 * t180 + t146 * t212) * MDP(14) + t244 * t177 * t236;
t176 = t177 ^ 2;
t156 = t195 * t177;
t143 = t146 * t213;
t139 = t214 + t222;
t134 = -t153 + t227;
t130 = t134 * t213;
t1 = [-MDP(6) * t230 + (-MDP(3) * t182 - MDP(4) * t185) * qJD(2) ^ 2 + t207 * (t180 * t193 - t183 * t194) + t206 * (t180 * t194 + t183 * t193) + t246 * t177 * t242 + (t129 * t162 + t134 * t141 + ((t138 * t183 - t139 * t180) * qJD(4) + t240) * t163 - (t138 * t180 + t139 * t183) * t242) * MDP(18); -t136 * MDP(6) + (-t170 * t215 - t184 * t199 + t219) * MDP(7) + (t143 + t220) * MDP(13) + (t130 + t220) * MDP(15) + t201 * MDP(16) + (t129 * t159 + t134 * t239) * MDP(18) + (-t196 * MDP(6) + t246 * t238) * t177 + (-t129 * MDP(17) + (t128 * t173 - t138 * t238) * MDP(18) + (t173 * t206 - MDP(11)) * t186 + (t196 * MDP(14) - MDP(17) * t239) * t177 + (t190 * MDP(13) + t191 * MDP(15) + (-MDP(16) - t211) * t139) * qJD(4)) * t180 + ((-t177 * t203 - t136 - t226) * MDP(13) + (-t149 * t177 - t129 - t226) * MDP(15) + (t126 * t173 - t139 * t238) * MDP(18) + ((t158 + t190) * MDP(14) + (-t134 - t158 - t191) * MDP(17) + t138 * t211) * qJD(4)) * t183 + t188; (-t136 + t228) * MDP(6) + (-t135 + t229) * MDP(7) + (t143 + t221) * MDP(13) + (t130 + t221) * MDP(15) + (-t229 * t237 + t201) * MDP(16) + (t129 * t164 + (-t154 + t155) * t134) * MDP(18) + (-t186 * MDP(11) + (-t228 + t235) * MDP(14) + (t192 + t228) * MDP(17) + (pkin(7) * t128 - t138 * t153) * MDP(18) + ((-MDP(13) * pkin(3) + MDP(15) * t164) * t177 + (-MDP(16) - t232) * t139) * qJD(4)) * t180 + ((-t136 - t235) * MDP(13) + t192 * MDP(15) + (pkin(7) * t126 - t139 * t153) * MDP(18) + ((t153 - t234) * MDP(14) + (-t134 - t153 - t227) * MDP(17) + t138 * t232) * qJD(4)) * t183 + t188; -t132 * MDP(14) + (qJD(5) * t236 + t132) * MDP(17) + (-t128 * pkin(4) + t126 * qJ(5) - t134 * t156 - t138 * t222 + t139 * t197) * MDP(18) - t207 * t131 - t244 * t176 + ((-t146 * MDP(13) - t134 * MDP(15) + (t139 - t214) * MDP(16) + t156 * MDP(17)) * t180 + (-t146 * MDP(14) + t156 * MDP(15) + (qJD(5) - t138 - t231) * MDP(16) + t134 * MDP(17)) * t183) * t177; -t176 * MDP(15) * t223 + (-t176 * t178 - t186) * MDP(17) + (t134 * t177 * t180 - qJD(4) * t139 + t128) * MDP(18);];
tauc = t1;
