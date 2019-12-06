% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:52
% EndTime: 2019-12-05 17:04:53
% DurationCPUTime: 0.41s
% Computational Cost: add. (337->89), mult. (735->133), div. (0->0), fcn. (364->6), ass. (0->58)
t104 = sin(qJ(5));
t107 = cos(qJ(5));
t135 = t107 * MDP(17);
t153 = t104 * MDP(16) + t135;
t154 = t153 * qJD(5);
t132 = -qJD(3) - qJD(4);
t152 = (t104 ^ 2 - t107 ^ 2) * MDP(12);
t108 = cos(qJ(4));
t105 = sin(qJ(4));
t106 = sin(qJ(3));
t149 = pkin(2) * qJD(2);
t129 = t106 * t149;
t124 = t105 * t129;
t101 = qJD(2) + qJD(3);
t109 = cos(qJ(3));
t128 = t109 * t149;
t93 = t101 * pkin(3) + t128;
t85 = -t108 * t93 + t124;
t145 = qJD(5) * t85;
t143 = t106 * t108;
t118 = t105 * t109 + t143;
t137 = qJD(4) * t108;
t126 = t106 * t137;
t111 = (t118 * qJD(3) + t126) * pkin(2);
t138 = qJD(4) * t105;
t127 = t93 * t138;
t78 = qJD(2) * t111 + t127;
t151 = t78 * t104 + t107 * t145;
t150 = t132 * t124;
t100 = qJD(2) - t132;
t98 = t109 * pkin(2) + pkin(3);
t148 = (t98 * t138 + t111) * t100;
t147 = t85 * t100;
t146 = (t105 * t93 + t108 * t129) * t100;
t144 = t105 * t106;
t110 = qJD(5) ^ 2;
t142 = t107 * t110;
t141 = t108 * MDP(9);
t139 = MDP(11) * t107;
t134 = t108 * MDP(10);
t133 = t110 * MDP(14);
t131 = MDP(13) * t142 + 0.2e1 * (t139 * t104 - t152) * qJD(5) * t100;
t125 = -(qJD(3) * t128 + qJD(4) * t93) * t108 - t150 - t147;
t123 = MDP(9) * t126;
t121 = pkin(6) * t110 - t146;
t120 = t110 * (pkin(2) * t143 + t105 * t98 + pkin(6)) + t148;
t117 = t108 * t109 - t144;
t79 = t98 * t137 + (t117 * qJD(3) - t106 * t138) * pkin(2);
t119 = qJD(5) * (t100 * (pkin(2) * t144 - t108 * t98) - t79);
t116 = -t93 * t137 - t150;
t114 = -t105 * MDP(9) - MDP(7) - t134;
t113 = t107 * MDP(16) - t104 * MDP(17) + MDP(9);
t82 = t104 * t145;
t112 = t82 * MDP(16) + t151 * MDP(17) + t131;
t99 = t100 ^ 2;
t97 = t105 * pkin(3) + pkin(6);
t89 = t117 * t149;
t1 = [-t153 * t110; (-t127 - t148) * MDP(9) + (-t79 * t100 + t116) * MDP(10) + ((-t120 - t78) * MDP(16) + MDP(17) * t119) * t107 + (MDP(16) * t119 + t120 * MDP(17) - t133) * t104 + (-qJD(2) * t123 + ((-t106 * MDP(6) - t109 * MDP(7)) * t101 + ((-MDP(6) - t141) * t106 + t114 * t109) * qJD(2)) * qJD(3)) * pkin(2) + t112; -MDP(9) * t127 + t116 * MDP(10) - t104 * t133 + (-t78 * t107 - t97 * t142 + t82) * MDP(16) + (t110 * t104 * t97 + t151) * MDP(17) + (t89 * MDP(10) + t113 * t118 * t149 + (-t108 * t154 + (-t105 * t113 - t134) * qJD(4)) * pkin(3)) * t100 + ((t101 * MDP(7) + qJD(3) * t114) * t109 + ((-qJD(3) + t101) * MDP(6) + t132 * t141) * t106) * t149 + t131 + (-pkin(3) * t137 + t89) * t154; (-t127 + t146) * MDP(9) + (t116 - t147) * MDP(10) + ((-t121 - t78) * MDP(16) - MDP(17) * t145) * t107 + (-MDP(16) * t145 + t121 * MDP(17) - t133) * t104 + (-t123 + (-t118 * MDP(9) - t109 * t134) * qJD(3)) * t149 + t112; t99 * t152 + t125 * t135 + (t125 * MDP(16) - t99 * t139) * t104;];
tauc = t1;
