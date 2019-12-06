% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:07
% EndTime: 2019-12-05 15:43:11
% DurationCPUTime: 1.00s
% Computational Cost: add. (728->141), mult. (1903->205), div. (0->0), fcn. (1468->6), ass. (0->73)
t167 = sin(pkin(9));
t170 = sin(qJ(4));
t168 = cos(pkin(9));
t172 = cos(qJ(4));
t201 = t168 * t172;
t180 = t167 * t170 - t201;
t146 = t180 * qJD(2);
t171 = cos(qJ(5));
t139 = t171 * t146;
t193 = t172 * qJD(4);
t197 = qJD(2) * t168;
t158 = t193 * t197;
t196 = qJD(2) * t170;
t191 = t167 * t196;
t142 = -qJD(4) * t191 + t158;
t152 = t167 * t172 + t168 * t170;
t149 = t152 * qJD(4);
t143 = qJD(2) * t149;
t147 = t152 * qJD(2);
t169 = sin(qJ(5));
t195 = qJD(5) * t169;
t111 = -qJD(5) * t139 + t171 * t142 - t169 * t143 - t147 * t195;
t184 = -t146 * t169 + t171 * t147;
t112 = qJD(5) * t184 + t142 * t169 + t171 * t143;
t126 = t147 * t169 + t139;
t166 = qJD(4) + qJD(5);
t203 = t126 * t166;
t204 = t184 * t166;
t215 = t126 * t184 * MDP(16) + (-t112 + t204) * MDP(19) + (-t126 ^ 2 + t184 ^ 2) * MDP(17) + (t111 + t203) * MDP(18);
t187 = (t167 ^ 2 + t168 ^ 2) * qJD(2);
t214 = MDP(7) * t187;
t178 = t152 * qJD(3);
t177 = qJD(2) * t178;
t205 = pkin(6) + qJ(3);
t156 = t205 * t167;
t163 = t168 * qJD(1);
t144 = -qJD(2) * t156 + t163;
t154 = qJ(3) * t197 + t167 * qJD(1);
t145 = pkin(6) * t197 + t154;
t185 = -t144 * t170 - t145 * t172;
t114 = -pkin(7) * t142 + qJD(4) * t185 - t177;
t122 = -pkin(7) * t146 - t185;
t161 = -pkin(3) * t168 - pkin(2);
t155 = qJD(2) * t161 + qJD(3);
t132 = pkin(4) * t146 + t155;
t213 = t132 * t126 + t122 * t195 + (-t122 * t166 - t114) * t169;
t210 = t172 * t144 - t145 * t170;
t209 = qJD(3) * t146;
t208 = qJD(5) - t166;
t113 = -pkin(7) * t143 + t210 * qJD(4) - t209;
t207 = -t169 * t113 + t171 * t114 - t132 * t184;
t206 = pkin(4) * t147;
t200 = t171 * t122;
t198 = qJD(2) * t167;
t121 = -pkin(7) * t147 + t210;
t118 = qJD(4) * pkin(4) + t121;
t190 = -pkin(4) * t166 - t118;
t183 = -t152 * t169 - t171 * t180;
t131 = t152 * t171 - t169 * t180;
t182 = (-qJ(3) * t198 + t163) * t167 - t154 * t168;
t157 = t205 * t168;
t181 = t156 * t170 - t157 * t172;
t176 = -t156 * t193 + qJD(3) * t201 + (-qJD(3) * t167 - qJD(4) * t157) * t170;
t174 = qJD(4) * t181 - t178;
t148 = t180 * qJD(4);
t137 = pkin(4) * t180 + t161;
t124 = -pkin(7) * t180 - t181;
t123 = -pkin(7) * t152 - t156 * t172 - t170 * t157;
t120 = pkin(7) * t148 + t174;
t119 = -pkin(7) * t149 + t176;
t116 = qJD(5) * t131 - t148 * t169 + t171 * t149;
t115 = qJD(5) * t183 - t148 * t171 - t149 * t169;
t1 = [(-t116 * MDP(21) - t115 * MDP(22)) * t166 + (-t149 * MDP(14) + t148 * MDP(15)) * qJD(4); (t142 * t152 - t147 * t148) * MDP(9) + (-t142 * t180 - t152 * t143 + t148 * t146 - t147 * t149) * MDP(10) + (t161 * t143 + t155 * t149) * MDP(14) + (t161 * t142 - t155 * t148) * MDP(15) + (t111 * t131 + t115 * t184) * MDP(16) + (t111 * t183 - t112 * t131 - t115 * t126 - t116 * t184) * MDP(17) + (t137 * t112 + t132 * t116 + (t126 * t149 - t143 * t183) * pkin(4)) * MDP(21) + (t137 * t111 + t132 * t115 + (t131 * t143 + t149 * t184) * pkin(4)) * MDP(22) + (t115 * MDP(18) - t116 * MDP(19) + (-t119 * t169 + t120 * t171 + (-t123 * t169 - t124 * t171) * qJD(5)) * MDP(21) + (-t119 * t171 - t120 * t169 - (t123 * t171 - t124 * t169) * qJD(5)) * MDP(22)) * t166 + (0.2e1 * t214 + (qJ(3) * t187 - t182) * MDP(8)) * qJD(3) + (-t148 * MDP(11) - t149 * MDP(12) + MDP(14) * t174 - MDP(15) * t176) * qJD(4); t158 * MDP(15) + (t112 + t204) * MDP(21) + (t111 - t203) * MDP(22) + ((t168 * t196 + t172 * t198 + t147) * MDP(14) + (-t146 - t191) * MDP(15)) * qJD(4) + (t182 * MDP(8) - t214) * qJD(2); t147 * t146 * MDP(9) + (-t146 ^ 2 + t147 ^ 2) * MDP(10) + (t158 + (t146 - t191) * qJD(4)) * MDP(11) + (-t155 * t147 - t177) * MDP(14) + (t155 * t146 + t209) * MDP(15) + (-t126 * t206 - (-t121 * t169 - t200) * t166 + (t169 * t190 - t200) * qJD(5) + t207) * MDP(21) + (-t184 * t206 + (qJD(5) * t190 + t121 * t166 - t113) * t171 + t213) * MDP(22) + t215; (t208 * (-t118 * t169 - t200) + t207) * MDP(21) + ((-t208 * t118 - t113) * t171 + t213) * MDP(22) + t215;];
tauc = t1;
