% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRRP6
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
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:43
% EndTime: 2019-12-31 16:30:44
% DurationCPUTime: 0.32s
% Computational Cost: add. (214->89), mult. (541->121), div. (0->0), fcn. (262->4), ass. (0->48)
t110 = -2 * qJD(3);
t75 = sin(qJ(3));
t78 = cos(qJ(2));
t91 = qJD(1) * qJD(2);
t88 = t78 * t91;
t69 = t75 * t88;
t76 = sin(qJ(2));
t94 = t76 * qJD(1);
t71 = qJD(2) * pkin(5) + t94;
t77 = cos(qJ(3));
t96 = qJD(3) * t77;
t60 = t71 * t96 + t69;
t103 = t77 * t71;
t92 = qJD(3) * qJ(4);
t64 = t92 + t103;
t109 = -t64 * qJD(3) + t60;
t79 = qJD(3) ^ 2;
t80 = qJD(2) ^ 2;
t108 = (t79 + t80) * t76;
t73 = t75 ^ 2;
t74 = t77 ^ 2;
t107 = (t73 - t74) * MDP(6);
t106 = (t73 + t74) * MDP(13);
t105 = t75 * t71;
t104 = t75 * t77;
t99 = qJD(2) * pkin(2);
t67 = -t77 * pkin(3) - t75 * qJ(4) - pkin(2);
t98 = qJD(2) * t67;
t97 = qJD(2) * t75;
t93 = t78 * qJD(1);
t90 = MDP(10) + MDP(12);
t89 = MDP(11) - MDP(14);
t72 = -t93 - t99;
t87 = t72 - t99;
t61 = -t93 + t98;
t86 = t61 + t98;
t85 = pkin(5) * MDP(15) + MDP(13);
t84 = -qJD(3) * pkin(3) + qJD(4);
t83 = qJD(2) * t78 * t110;
t82 = pkin(3) * t75 - qJ(4) * t77;
t62 = qJD(3) * t82 - t75 * qJD(4);
t56 = (t62 + t94) * qJD(2);
t81 = -qJD(2) * t62 + t76 * t91 - t56;
t70 = t77 * t88;
t66 = t82 * qJD(2);
t63 = t84 + t105;
t57 = t70 + (qJD(4) - t105) * qJD(3);
t1 = [t90 * (-t77 * t108 + t75 * t83) + t89 * (t75 * t108 + t77 * t83) + ((-t56 + (t63 * t75 + t64 * t77) * qJD(2)) * MDP(15) + (-MDP(4) + t106) * t80) * t78 + (-t80 * MDP(3) + (qJD(2) * t61 + t109 * t75 + t57 * t77 + t63 * t96) * MDP(15)) * t76; (t56 * t67 + (t62 - t94) * t61) * MDP(15) + (-t93 * t106 + t107 * t110) * qJD(2) + (t81 * MDP(12) + t57 * MDP(13) + (pkin(5) * t57 - t64 * t93) * MDP(15) + (-pkin(5) * t90 + MDP(7)) * t79 + (0.2e1 * MDP(5) * t97 + (t87 + t93) * MDP(11) + (-t86 - t93) * MDP(14) + t85 * t63) * qJD(3)) * t77 + (t60 * MDP(13) + t81 * MDP(14) + (pkin(5) * t60 - t63 * t93) * MDP(15) + (pkin(5) * t89 - MDP(8)) * t79 + (MDP(10) * t87 + MDP(12) * t86 - t64 * t85 + t90 * t93) * qJD(3)) * t75; -t70 * MDP(11) + (0.2e1 * qJD(3) * qJD(4) + t70) * MDP(14) + (-t63 * t103 - t60 * pkin(3) + t57 * qJ(4) - t61 * t66 + (qJD(4) + t105) * t64) * MDP(15) - t90 * t69 + (-MDP(5) * t104 + t107) * t80 + ((-t72 * MDP(10) - t61 * MDP(12) + (t64 - t92) * MDP(13) + t66 * MDP(14)) * t75 + (-t72 * MDP(11) + t66 * MDP(12) + (-t63 + t84) * MDP(13) + t61 * MDP(14)) * t77) * qJD(2); -t80 * MDP(12) * t104 + (-t73 * t80 - t79) * MDP(14) + (t61 * t97 + t109) * MDP(15);];
tauc = t1;
