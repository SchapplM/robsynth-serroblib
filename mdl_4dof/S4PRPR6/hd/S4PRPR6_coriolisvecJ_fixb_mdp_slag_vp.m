% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:42
% EndTime: 2019-12-31 16:24:43
% DurationCPUTime: 0.37s
% Computational Cost: add. (189->72), mult. (529->116), div. (0->0), fcn. (342->6), ass. (0->42)
t78 = qJD(2) ^ 2;
t72 = sin(pkin(7));
t73 = cos(pkin(7));
t94 = t72 ^ 2 + t73 ^ 2;
t100 = t94 * t78 * MDP(7);
t74 = sin(qJ(4));
t76 = cos(qJ(4));
t61 = t76 * t72 + t74 * t73;
t58 = t61 * qJD(4);
t75 = sin(qJ(2));
t90 = t75 * qJD(1);
t67 = qJD(2) * qJ(3) + t90;
t99 = t94 * t67 - t90;
t98 = t74 * t72;
t96 = t76 * t73;
t95 = pkin(5) + qJ(3);
t93 = qJD(2) * pkin(2);
t77 = cos(qJ(2));
t91 = qJD(1) * t77;
t56 = qJD(2) * t61;
t89 = t75 * qJD(4) ^ 2;
t87 = qJD(2) * t98;
t86 = qJD(2) * t96;
t69 = -pkin(3) * t73 - pkin(2);
t85 = t94 * t77;
t62 = (qJD(3) + t91) * qJD(2);
t84 = t94 * t62;
t83 = t94 * qJD(3);
t81 = qJD(3) - t91;
t60 = -t96 + t98;
t57 = t60 * qJD(4);
t80 = t77 * t58;
t79 = t77 * t57;
t66 = qJD(4) * t86;
t65 = t81 - t93;
t64 = t95 * t73;
t63 = t95 * t72;
t59 = t69 * qJD(2) + t81;
t54 = -t86 + t87;
t53 = qJD(2) * t58;
t52 = -qJD(4) * t87 + t66;
t1 = [t77 * t100 + (t75 * t84 + (t65 * t75 + t99 * t77) * qJD(2)) * MDP(8) + (-t77 * t53 + t60 * t89 + (t75 * t54 - t80) * qJD(2)) * MDP(14) + (-t77 * t52 + t61 * t89 + (t75 * t56 + t79) * qJD(2)) * MDP(15) + (-t77 * MDP(4) + (-t73 * MDP(5) + t72 * MDP(6) - MDP(3)) * t75) * t78; (t84 + (-qJD(1) * t85 + t83) * qJD(2)) * MDP(7) + (t67 * t83 + qJ(3) * t84 + ((-t65 - t93) * t75 - t67 * t85) * qJD(1)) * MDP(8) + (t52 * t61 - t56 * t57) * MDP(9) + (-t52 * t60 - t53 * t61 + t54 * t57 - t56 * t58) * MDP(10) - t57 * qJD(4) * MDP(11) - t58 * qJD(4) * MDP(12) + (t69 * t53 + t59 * t58 + ((t63 * t74 - t64 * t76) * qJD(4) - t61 * qJD(3)) * qJD(4) + (t80 + (qJD(2) * t60 - t54) * t75) * qJD(1)) * MDP(14) + (t69 * t52 - t59 * t57 + ((t63 * t76 + t64 * t74) * qJD(4) + t60 * qJD(3)) * qJD(4) - t79 * qJD(1)) * MDP(15); t56 * qJD(4) * MDP(14) + (-qJD(4) * t54 + t66) * MDP(15) - t100 + (-t99 * MDP(8) + (t61 * MDP(14) - MDP(15) * t98) * qJD(4)) * qJD(2); t56 * t54 * MDP(9) + (-t54 ^ 2 + t56 ^ 2) * MDP(10) + (t66 + (t54 - t87) * qJD(4)) * MDP(11) + (-t59 * t56 - t61 * t62) * MDP(14) + (t59 * t54 + t60 * t62) * MDP(15);];
tauc = t1;
