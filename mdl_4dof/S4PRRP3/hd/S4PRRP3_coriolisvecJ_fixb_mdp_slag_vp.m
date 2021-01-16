% Calculate Coriolis joint torque vector for
% S4PRRP3
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
%   see S4PRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:27
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:27:06
% EndTime: 2021-01-14 22:27:08
% DurationCPUTime: 0.32s
% Computational Cost: add. (172->75), mult. (448->112), div. (0->0), fcn. (192->2), ass. (0->43)
t88 = (pkin(2) * MDP(11));
t60 = sin(qJ(3));
t58 = t60 ^ 2;
t61 = cos(qJ(3));
t59 = t61 ^ 2;
t87 = (t58 - t59) * MDP(6);
t75 = qJD(2) * t60;
t55 = qJD(3) * pkin(5) * t75;
t86 = 2 * MDP(13);
t85 = pkin(3) * t58;
t84 = t61 * pkin(3);
t56 = -pkin(2) - t84;
t52 = qJD(2) * t56 + qJD(4);
t83 = t52 * t60;
t63 = qJD(2) ^ 2;
t82 = t63 * t61;
t81 = -qJ(4) - pkin(5);
t53 = t81 * t60;
t57 = t61 * qJD(1);
t48 = qJD(2) * t53 + t57;
t76 = qJD(3) * pkin(3);
t47 = t48 + t76;
t80 = t47 - t48;
t78 = MDP(15) * pkin(3);
t77 = pkin(2) * MDP(10);
t74 = qJD(2) * t61;
t73 = t60 * qJD(1);
t72 = t61 * qJD(4);
t71 = -qJD(4) - t52;
t70 = qJD(3) * qJD(1);
t69 = t60 * t82;
t67 = qJD(3) * t60 * qJ(4);
t54 = t81 * t61;
t66 = qJD(3) * t81;
t65 = t61 * t66;
t50 = -qJD(2) * t54 + t73;
t64 = t47 * t60 - t50 * t61;
t51 = -t60 * qJD(4) + t65;
t62 = qJD(3) ^ 2;
t49 = t60 * t66 + t72;
t46 = qJD(2) * t51 - t60 * t70;
t45 = t61 * t70 - t55 + (-t67 + t72) * qJD(2);
t1 = [(-t64 * qJD(3) + t45 * t60 + t46 * t61) * MDP(15) + ((-MDP(11) - MDP(13)) * t61 + (-MDP(10) - MDP(12)) * t60) * t62; (t45 * t61 - t46 * t60 + t49 * t74 - t51 * t75) * MDP(14) + (-t45 * t54 + t46 * t53 + t47 * t51 + t50 * t49) * MDP(15) + (t61 * MDP(7) - t60 * MDP(8) + (-MDP(10) * t61 + MDP(11) * t60) * pkin(5)) * t62 + ((t51 + t83) * MDP(12) + (t52 * t61 - t49) * MDP(13) + (-t47 * t61 - t50 * t60) * MDP(14) + t78 * t83 + (-0.2e1 * t87 + t85 * t86 + (t56 * MDP(13) - t53 * MDP(14) - (2 * t88)) * t61 + (0.2e1 * t61 * MDP(5) - 0.2e1 * t77 + (t56 - 0.2e1 * t84) * MDP(12) + t54 * MDP(14) + t56 * t78) * t60) * qJD(2)) * qJD(3); -MDP(5) * t69 + t82 * t88 + (pkin(3) * t69 + (t50 - t73) * qJD(3) + (t71 * t60 + t65) * qJD(2)) * MDP(12) + (t55 + (t48 - t57) * qJD(3) + (t71 * t61 + t67) * qJD(2)) * MDP(13) + (-t76 + t80) * MDP(14) * t74 + (t80 * t50 + (-t52 * t75 + t46) * pkin(3)) * MDP(15) + (-MDP(13) * t85 + t60 * t77 + t87) * t63; (-t58 - t59) * MDP(14) * t63 + (t64 * MDP(15) + (t61 * t86 + (0.2e1 * MDP(12) + t78) * t60) * qJD(3)) * qJD(2);];
tauc = t1;
