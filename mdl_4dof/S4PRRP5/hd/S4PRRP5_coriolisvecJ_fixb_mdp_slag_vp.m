% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:23
% EndTime: 2019-12-31 16:29:24
% DurationCPUTime: 0.29s
% Computational Cost: add. (175->67), mult. (457->111), div. (0->0), fcn. (227->4), ass. (0->41)
t71 = sin(qJ(3));
t69 = t71 ^ 2;
t73 = cos(qJ(3));
t70 = t73 ^ 2;
t100 = (t69 - t70) * MDP(6);
t72 = sin(qJ(2));
t88 = t72 * qJD(1);
t65 = qJD(2) * pkin(5) + t88;
t80 = qJ(4) * qJD(2) + t65;
t57 = t80 * t71;
t92 = qJD(3) * pkin(3);
t56 = -t57 + t92;
t58 = t80 * t73;
t78 = t56 * t71 - t58 * t73;
t99 = t78 * MDP(13);
t82 = (t69 + t70) * MDP(12);
t98 = -qJ(4) - pkin(5);
t97 = t56 + t57;
t89 = qJD(3) * t71;
t62 = (pkin(3) * t89 + t88) * qJD(2);
t75 = qJD(3) ^ 2;
t76 = qJD(2) ^ 2;
t94 = t75 + t76;
t93 = qJD(2) * pkin(2);
t84 = -t73 * pkin(3) - pkin(2);
t74 = cos(qJ(2));
t87 = t74 * qJD(1);
t91 = MDP(13) * (t84 * qJD(2) + qJD(4) - t87);
t90 = qJD(2) * t73;
t86 = qJ(4) * qJD(3);
t85 = pkin(3) * t91;
t83 = qJD(3) * t98;
t79 = qJD(4) + t87;
t77 = -0.2e1 * t93;
t64 = t98 * t73;
t63 = t98 * t71;
t60 = -t71 * qJD(4) + t73 * t83;
t59 = t73 * qJD(4) + t71 * t83;
t55 = -t73 * t65 * qJD(3) + (-t79 * t71 - t73 * t86) * qJD(2);
t54 = -t65 * t89 + (-t71 * t86 + t79 * t73) * qJD(2);
t1 = [(-t62 * MDP(13) + (-MDP(4) + t82) * t76 + (-t99 + 0.2e1 * (-t71 * MDP(10) - t73 * MDP(11)) * qJD(3)) * qJD(2)) * t74 + (qJD(2) * t91 - t76 * MDP(3) + (-t94 * MDP(10) + (-qJD(3) * t56 + t54) * MDP(13)) * t73 + (t94 * MDP(11) + (-qJD(3) * t58 - t55) * MDP(13)) * t71) * t72; (t54 * t73 + t59 * t90 + (-qJD(2) * t60 - t55) * t71) * MDP(12) + (-t54 * t64 + t55 * t63 + t56 * t60 + t58 * t59 + t62 * t84) * MDP(13) + (-t72 * t91 + (-qJD(2) * t82 + t99) * t74) * qJD(1) + (t73 * MDP(7) - t71 * MDP(8) + (-t73 * MDP(10) + t71 * MDP(11)) * pkin(5)) * t75 + (-0.2e1 * qJD(2) * t100 + (t77 * MDP(11) + (-qJD(2) * t63 - t56) * MDP(12)) * t73 + (0.2e1 * MDP(5) * t90 + t77 * MDP(10) + (qJD(2) * t64 - t58) * MDP(12) + t85) * t71) * qJD(3); (t55 * pkin(3) + t97 * t58) * MDP(13) + (-t71 * t73 * MDP(5) + t100) * t76 + ((MDP(10) * t93 - t85) * t71 + (t93 * MDP(11) + (-t92 + t97) * MDP(12)) * t73) * qJD(2); (t78 * qJD(2) + t62) * MDP(13) - t76 * t82;];
tauc = t1;
