% Calculate Coriolis joint torque vector for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:56
% EndTime: 2022-01-23 09:20:57
% DurationCPUTime: 0.57s
% Computational Cost: add. (617->107), mult. (1292->172), div. (0->0), fcn. (699->8), ass. (0->69)
t141 = qJD(1) + qJD(3);
t136 = cos(pkin(8)) * pkin(1) + pkin(2);
t132 = t136 * qJD(1);
t149 = sin(qJ(3));
t199 = pkin(1) * sin(pkin(8));
t171 = qJD(3) * t199;
t165 = qJD(1) * t171;
t151 = cos(qJ(3));
t177 = qJD(3) * t151;
t160 = t132 * t177 - t149 * t165;
t113 = t141 * qJD(4) + t160;
t144 = sin(pkin(9));
t139 = t144 ^ 2;
t146 = cos(pkin(9));
t180 = t146 ^ 2 + t139;
t206 = t180 * t113;
t205 = t146 * MDP(8) + MDP(6);
t148 = sin(qJ(5));
t150 = cos(qJ(5));
t189 = t139 * t150;
t204 = -t148 * MDP(11) * t189 + (t148 ^ 2 - t150 ^ 2) * MDP(12) * t139;
t200 = -t136 * t151 + t149 * t199;
t172 = qJD(1) * t199;
t120 = t151 * t132 - t149 * t172;
t159 = qJD(4) - t120;
t183 = t149 * t132;
t121 = t151 * t172 + t183;
t116 = t141 * qJ(4) + t121;
t105 = -t146 * qJD(2) + t116 * t144;
t198 = t105 * t144;
t155 = t136 * t177 - t149 * t171;
t122 = qJD(4) + t155;
t196 = t122 * t141;
t187 = t141 * t146;
t133 = -qJD(5) + t187;
t194 = t133 * t148;
t193 = t133 * t150;
t191 = t139 * t141;
t190 = t139 * t148;
t188 = t141 * t144;
t185 = t146 * t148;
t184 = t146 * t150;
t129 = -t146 * pkin(4) - t144 * pkin(7) - pkin(3);
t103 = t129 * t141 + t159;
t106 = qJD(2) * t144 + t116 * t146;
t119 = qJD(3) * t183 + t151 * t165;
t182 = (t113 * t184 + t148 * t119 + (t103 * t150 - t106 * t148) * qJD(5)) * t146 + t113 * t189;
t175 = qJD(5) * t148;
t174 = qJD(5) * t150;
t173 = qJD(5) + t133;
t170 = t105 * t188;
t169 = t141 * t189;
t168 = t144 * t175;
t167 = t144 * t174;
t127 = t168 * t187;
t164 = (t133 * t168 + t127) * MDP(13) + (t133 + t187) * MDP(14) * t167 + 0.2e1 * t204 * qJD(5) * t141;
t158 = -t103 * t148 - t106 * t150;
t102 = t158 * qJD(5) - t113 * t185 + t150 * t119;
t162 = -t102 * t146 + t105 * t167 + t113 * t190;
t157 = t106 * t146 + t198;
t156 = t133 * t146 + t191;
t154 = t149 * t136 + t151 * t199;
t153 = qJ(4) * t174 + t159 * t148;
t138 = t141 ^ 2;
t124 = qJ(4) + t154;
t123 = t154 * qJD(3);
t117 = t129 + t200;
t114 = -t141 * pkin(3) + t159;
t1 = [(-t155 * t141 - t160) * MDP(7) + (t180 * t196 + t206) * MDP(9) + (t119 * (-pkin(3) + t200) + t114 * t123 + t157 * t122 + t124 * t206) * MDP(10) + (-(-t122 * t185 + t123 * t150) * t133 + t190 * t196 + (-(-t117 * t148 - t124 * t184) * t133 + t124 * t169) * qJD(5) + t162) * MDP(16) + ((t122 * t184 + t123 * t148) * t133 + t122 * t169 + (t117 * t193 + (-t156 * t124 - t198) * t148) * qJD(5) + t182) * MDP(17) + t164 + t205 * (-t123 * t141 - t119); t127 * MDP(17) + (-MDP(17) * t194 + (t133 - t187) * MDP(16) * t150) * t144 * qJD(5); (t120 * t141 - t160) * MDP(7) + (t141 * t159 * t180 + t206) * MDP(9) + (-pkin(3) * t119 + qJ(4) * t206 - t114 * t121 + t157 * t159) * MDP(10) + ((t150 * t121 + t129 * t175 + t153 * t146) * t133 + t153 * t191 + t162) * MDP(16) + (-t121 * t194 + t159 * t156 * t150 + (t129 * t193 + (-t156 * qJ(4) - t198) * t148) * qJD(5) + t182) * MDP(17) + t164 + t205 * (t121 * t141 - t119); (-t157 * t141 + t119) * MDP(10) - t180 * MDP(9) * t138 + (t148 * MDP(16) + t150 * MDP(17)) * (-t133 ^ 2 - t139 * t138); (t158 * t133 - t150 * t170 + t102) * MDP(16) + ((-t173 * t103 - t146 * t113) * t150 + (t173 * t106 - t119 + t170) * t148) * MDP(17) - (t148 * MDP(13) + t150 * MDP(14)) * t173 * t188 - t204 * t138;];
tauc = t1;
