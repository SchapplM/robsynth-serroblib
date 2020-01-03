% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:41
% EndTime: 2019-12-31 16:41:43
% DurationCPUTime: 0.34s
% Computational Cost: add. (168->62), mult. (411->98), div. (0->0), fcn. (238->4), ass. (0->34)
t71 = sin(pkin(6));
t72 = cos(pkin(6));
t88 = t71 ^ 2 + t72 ^ 2;
t93 = qJD(3) * t88;
t96 = qJ(2) * MDP(6) + t71 * MDP(7) + t72 * MDP(8) + MDP(5);
t74 = sin(qJ(4));
t75 = cos(qJ(4));
t89 = t75 * t72;
t95 = -t74 * t71 + t89;
t73 = -pkin(1) - qJ(3);
t94 = qJD(1) * t73;
t91 = -pkin(5) + t73;
t57 = t75 * t71 + t74 * t72;
t84 = qJD(1) * t57;
t83 = qJD(1) * t71;
t70 = qJD(1) * qJ(2);
t67 = qJD(3) + t70;
t82 = t74 * t83;
t81 = qJD(1) * t89;
t55 = t57 * qJD(4);
t78 = t57 * qJD(3);
t77 = t95 * qJD(3);
t76 = qJD(1) ^ 2;
t66 = t71 * pkin(3) + qJ(2);
t63 = qJD(2) + t94;
t62 = qJD(4) * t82;
t61 = pkin(3) * t83 + t67;
t60 = t91 * t72;
t59 = t91 * t71;
t56 = t95 * qJD(4);
t54 = t81 - t82;
t49 = qJD(4) * t81 - t62;
t48 = qJD(1) * t55;
t1 = [((t67 + t70) * qJD(2) + (-t63 - t94) * t93) * MDP(10) + (-t48 * t95 - t54 * t55) * MDP(11) + (t48 * t57 - t49 * t95 - t54 * t56 + t55 * t84) * MDP(12) - t55 * qJD(4) * MDP(13) - t56 * qJD(4) * MDP(14) + (t66 * t49 + t61 * t56 + ((-t59 * t75 - t60 * t74) * qJD(4) - t77) * qJD(4) + 0.2e1 * t84 * qJD(2)) * MDP(16) + (-t66 * t48 - t61 * t55 + ((t59 * t74 - t60 * t75) * qJD(4) + t78) * qJD(4) + (qJD(1) * t95 + t54) * qJD(2)) * MDP(17) + 0.2e1 * (MDP(9) * t93 + t96 * qJD(2)) * qJD(1); (-t55 * MDP(16) - t56 * MDP(17)) * qJD(4) + ((-t67 - t93) * MDP(10) - t84 * MDP(16) - t54 * MDP(17)) * qJD(1) - t96 * t76; (t54 * qJD(4) - t62) * MDP(16) - t84 * qJD(4) * MDP(17) - t88 * MDP(9) * t76 + ((t88 * t63 + qJD(2)) * MDP(10) + (MDP(16) * t89 - t57 * MDP(17)) * qJD(4)) * qJD(1); t54 * t84 * MDP(11) + (t54 ^ 2 - t84 ^ 2) * MDP(12) + (t62 + (t54 - t81) * qJD(4)) * MDP(14) + (-qJD(1) * t77 - t61 * t54) * MDP(16) + (qJD(1) * t78 + t61 * t84) * MDP(17);];
tauc = t1;
