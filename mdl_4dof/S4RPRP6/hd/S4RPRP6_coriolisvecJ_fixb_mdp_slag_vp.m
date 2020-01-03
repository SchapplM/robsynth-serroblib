% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:15
% EndTime: 2019-12-31 16:46:15
% DurationCPUTime: 0.26s
% Computational Cost: add. (167->61), mult. (343->95), div. (0->0), fcn. (124->2), ass. (0->37)
t68 = -pkin(1) - pkin(5);
t60 = qJD(1) * t68 + qJD(2);
t91 = -qJ(4) * qJD(1) + t60;
t67 = cos(qJ(3));
t80 = qJD(3) * t67;
t90 = pkin(3) * t80 + qJD(2);
t66 = sin(qJ(3));
t64 = t66 ^ 2;
t65 = t67 ^ 2;
t89 = (t64 - t65) * MDP(8);
t55 = t91 * t67;
t84 = qJD(3) * pkin(3);
t51 = t55 + t84;
t88 = t51 - t55;
t87 = t90 * qJD(1);
t69 = qJD(3) ^ 2;
t70 = qJD(1) ^ 2;
t85 = -t69 - t70;
t83 = t67 * MDP(7);
t82 = qJ(4) - t68;
t81 = qJD(3) * t66;
t72 = pkin(3) * t66 + qJ(2);
t57 = qJD(1) * t72 + qJD(4);
t79 = t57 * qJD(1);
t78 = t66 * qJD(4);
t77 = t67 * MDP(12);
t76 = t67 * qJD(4);
t74 = qJ(4) * qJD(3);
t59 = t82 * t67;
t71 = -MDP(6) * qJ(2) - MDP(5);
t58 = t82 * t66;
t54 = t91 * t66;
t53 = -qJD(3) * t59 - t78;
t52 = t81 * t82 - t76;
t50 = t60 * t80 + (-t67 * t74 - t78) * qJD(1);
t49 = -t60 * t81 + (t66 * t74 - t76) * qJD(1);
t1 = [(-t49 * t67 - t50 * t66 + t51 * t81 - t54 * t80) * MDP(14) + (-t49 * t59 - t50 * t58 + t51 * t52 + t54 * t53 + t57 * t90 + t87 * t72) * MDP(15) + ((-MDP(13) * t68 - MDP(10)) * t67 + (-MDP(12) * t68 - MDP(9)) * t66) * t69 + ((-t52 * t67 - t53 * t66) * MDP(14) + 0.2e1 * (MDP(12) * t66 + MDP(13) * t67 - t71) * qJD(2) + (-0.2e1 * t66 * t83 + 0.2e1 * t89 + (t58 * t67 - t59 * t66) * MDP(14) + 0.2e1 * (-MDP(13) * t66 + t77) * qJ(2)) * qJD(3)) * qJD(1); -MDP(15) * t79 + t71 * t70 + (t85 * MDP(13) + (qJD(3) * t54 + t49) * MDP(15)) * t67 + (t85 * MDP(12) + (-qJD(3) * t51 + t50) * MDP(15)) * t66; (t88 * t54 + (-t67 * t79 + t49) * pkin(3)) * MDP(15) + (-qJ(2) * t77 - t89) * t70 + ((MDP(13) * qJ(2) + t83) * t70 + (t84 - t88) * MDP(14) * qJD(1)) * t66; ((t51 * t67 + t54 * t66) * qJD(1) + t87) * MDP(15) + (-t64 - t65) * MDP(14) * t70;];
tauc = t1;
