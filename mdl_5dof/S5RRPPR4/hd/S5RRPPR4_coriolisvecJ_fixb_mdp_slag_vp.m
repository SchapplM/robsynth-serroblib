% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:49
% EndTime: 2019-12-31 19:27:51
% DurationCPUTime: 0.41s
% Computational Cost: add. (469->110), mult. (779->148), div. (0->0), fcn. (357->6), ass. (0->64)
t173 = 2 * qJD(5);
t126 = qJD(1) + qJD(2);
t130 = sin(pkin(8));
t131 = cos(pkin(8));
t135 = cos(qJ(2));
t169 = pkin(1) * qJD(1);
t155 = t135 * t169;
t133 = sin(qJ(2));
t156 = t133 * t169;
t99 = t130 * t155 - t131 * t156;
t146 = (t130 * qJD(3) - t99) * t126;
t168 = pkin(1) * qJD(2);
t154 = qJD(1) * t168;
t118 = t135 * t154;
t123 = t126 * qJD(3);
t105 = t118 + t123;
t149 = t133 * t154;
t90 = t130 * t105 - t131 * t149;
t172 = t146 + t90;
t171 = MDP(5) + MDP(7);
t132 = sin(qJ(5));
t134 = cos(qJ(5));
t170 = (t132 ^ 2 - t134 ^ 2) * MDP(14);
t136 = -pkin(2) - pkin(3);
t157 = t135 * t168;
t117 = qJD(3) + t157;
t158 = t133 * t168;
t96 = t130 * t117 - t131 * t158;
t167 = t96 * t126;
t137 = qJD(5) ^ 2;
t164 = t131 * qJ(3) + t130 * t136;
t166 = (-pkin(7) + t164) * t137;
t91 = t131 * t105 + t130 * t149;
t153 = -t135 * pkin(1) - pkin(2);
t119 = -pkin(3) + t153;
t120 = t133 * pkin(1) + qJ(3);
t165 = t130 * t119 + t131 * t120;
t162 = MDP(12) * t131;
t161 = MDP(13) * t134;
t160 = t134 * MDP(19);
t159 = t126 * t173;
t110 = t126 * qJ(3) + t156;
t148 = qJD(3) - t155;
t98 = t136 * t126 + t148;
t87 = -t130 * t110 + t131 * t98;
t85 = t126 * pkin(4) - t87;
t152 = t85 * t126 - t91;
t100 = (t130 * t133 + t131 * t135) * t169;
t150 = t131 * qJD(3) - t100;
t147 = -t137 * (-pkin(7) + t165) + t167;
t145 = -t130 * qJ(3) + t131 * t136;
t144 = t131 * t119 - t130 * t120;
t143 = -MDP(18) * t134 + MDP(19) * t132;
t142 = t132 * MDP(18) + t160;
t141 = t126 * t155 - t118;
t97 = t131 * t117 + t130 * t158;
t140 = qJD(5) * (-t126 * (pkin(4) - t144) - t85 - t97);
t139 = -t137 * t134 * MDP(15) - t159 * t170 + ((MDP(16) * t137) + t159 * t161) * t132;
t138 = qJD(5) * (-(pkin(4) - t145) * t126 - t150 - t85);
t125 = t126 ^ 2;
t107 = -t126 * pkin(2) + t148;
t89 = t90 * t134;
t88 = t131 * t110 + t130 * t98;
t1 = [(-t126 * t157 - t118) * MDP(6) + (t117 * t126 + t105) * MDP(8) + (t105 * t120 + t110 * t117 + (t153 * qJD(1) + t107) * t158) * MDP(9) + (t90 + t167) * MDP(10) + (t97 * t126 + t91) * MDP(11) + (-t90 * t144 + t91 * t165 - t87 * t96 + t88 * t97) * MDP(12) + (t132 * t140 + t147 * t134 + t89) * MDP(18) + ((-t147 - t90) * t132 + t134 * t140) * MDP(19) + t139 + t171 * (-qJD(1) - t126) * t158; t141 * MDP(6) + (0.2e1 * t123 - t141) * MDP(8) + (t105 * qJ(3) + t110 * qJD(3) + (-t110 * t135 + (-pkin(2) * qJD(2) - t107) * t133) * t169) * MDP(9) + t172 * MDP(10) + (t150 * t126 + t91) * MDP(11) + (t91 * t164 - t90 * t145 - t88 * t100 + t87 * t99 + (-t87 * t130 + t88 * t131) * qJD(3)) * MDP(12) + (t89 + (t146 - t166) * t134 + t132 * t138) * MDP(18) + ((t166 - t172) * t132 + t134 * t138) * MDP(19) + t139 + t171 * (-qJD(2) + t126) * t156; MDP(9) * t149 - t90 * t162 + (t91 * MDP(12) + t143 * t137) * t130 + (-t110 * MDP(9) - t88 * t162 + (-t131 * MDP(11) - MDP(8)) * t126 + (t87 * MDP(12) + (-MDP(10) + t143) * t126) * t130 + t142 * t131 * t173) * t126; -t142 * t137; t152 * t160 + t125 * t170 + (t152 * MDP(18) - t125 * t161) * t132;];
tauc = t1;
