% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:34
% EndTime: 2019-12-31 17:01:35
% DurationCPUTime: 0.25s
% Computational Cost: add. (190->55), mult. (463->92), div. (0->0), fcn. (236->6), ass. (0->44)
t92 = sin(qJ(4));
t94 = cos(qJ(4));
t124 = (t92 ^ 2 - t94 ^ 2) * MDP(9);
t123 = t92 * MDP(13);
t93 = sin(qJ(2));
t95 = cos(qJ(2));
t122 = t93 * MDP(5) + t95 * MDP(6);
t90 = sin(pkin(7));
t121 = t90 * t93;
t91 = cos(pkin(7));
t120 = t91 * t93;
t85 = t95 * pkin(1) + pkin(2);
t119 = pkin(1) * t120 + t90 * t85;
t117 = MDP(8) * t94;
t116 = pkin(1) * qJD(1);
t115 = pkin(1) * qJD(2);
t106 = t93 * t116;
t105 = t95 * t116;
t87 = qJD(1) + qJD(2);
t78 = t87 * pkin(2) + t105;
t64 = -t90 * t106 + t91 * t78;
t62 = -t87 * pkin(3) - t64;
t112 = qJD(4) * t62;
t111 = t94 * MDP(14);
t96 = qJD(4) ^ 2;
t110 = t96 * MDP(11);
t98 = pkin(1) * (t91 * t95 - t121);
t75 = qJD(2) * t98;
t69 = qJD(1) * t75;
t104 = -t62 * t87 - t69;
t73 = (t90 * t95 + t120) * t115;
t103 = (pkin(6) + t119) * t96 + t73 * t87;
t80 = t91 * t106;
t72 = t90 * t105 + t80;
t102 = -t72 * t87 + (t90 * pkin(2) + pkin(6)) * t96;
t99 = -pkin(1) * t121 + t91 * t85;
t101 = qJD(4) * ((-pkin(3) - t99) * t87 - t75);
t74 = qJD(1) * t98;
t100 = qJD(4) * ((-t91 * pkin(2) - pkin(3)) * t87 + t74);
t68 = qJD(1) * t73;
t97 = t96 * t94 * MDP(10) + t112 * t123 + (t94 * t112 + t68 * t92) * MDP(14) + 0.2e1 * (t117 * t92 - t124) * qJD(4) * t87;
t86 = t87 ^ 2;
t65 = t90 * t78 + t80;
t1 = [(t69 * t119 - t64 * t73 + t65 * t75 - t68 * t99) * MDP(7) + ((-t103 - t68) * MDP(13) + MDP(14) * t101) * t94 + (MDP(13) * t101 + t103 * MDP(14) - t110) * t92 + t97 + t122 * t115 * (-qJD(1) - t87); (t64 * t72 - t65 * t74 + (-t68 * t91 + t69 * t90) * pkin(2)) * MDP(7) + ((-t102 - t68) * MDP(13) + MDP(14) * t100) * t94 + (MDP(13) * t100 + t102 * MDP(14) - t110) * t92 + t97 + t122 * t116 * (-qJD(2) + t87); (-t111 - t123) * t96; t104 * t111 + t86 * t124 + (t104 * MDP(13) - t86 * t117) * t92;];
tauc = t1;
