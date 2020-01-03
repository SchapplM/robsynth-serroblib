% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:57
% EndTime: 2019-12-31 16:37:58
% DurationCPUTime: 0.28s
% Computational Cost: add. (146->45), mult. (383->79), div. (0->0), fcn. (248->6), ass. (0->28)
t63 = sin(pkin(6)) * pkin(1) + qJ(3);
t68 = sin(pkin(7));
t70 = cos(pkin(7));
t79 = qJD(1) * (t68 ^ 2 + t70 ^ 2);
t88 = t63 * t79 * MDP(8);
t87 = MDP(7) * t79;
t72 = sin(qJ(4));
t73 = cos(qJ(4));
t59 = t68 * t73 + t70 * t72;
t53 = t59 * qJD(1);
t85 = pkin(5) + t63;
t83 = qJD(1) * t72;
t82 = qJD(1) * t73;
t80 = t68 * t83;
t77 = t68 * t72 - t70 * t73;
t60 = -cos(pkin(6)) * pkin(1) - pkin(2) - t70 * pkin(3);
t76 = t59 * qJD(3);
t55 = t59 * qJD(4);
t75 = t77 * qJD(3);
t62 = qJD(4) * t70 * t82;
t57 = t85 * t70;
t56 = t85 * t68;
t54 = t77 * qJD(4);
t52 = t77 * qJD(1);
t51 = qJD(1) * t60 + qJD(3);
t50 = qJD(1) * t55;
t49 = -qJD(4) * t80 + t62;
t1 = [(t49 * t59 - t53 * t54) * MDP(9) + (-t49 * t77 - t59 * t50 + t54 * t52 - t53 * t55) * MDP(10) - t54 * qJD(4) * MDP(11) - t55 * qJD(4) * MDP(12) + (t60 * t50 + t51 * t55 + ((t56 * t72 - t57 * t73) * qJD(4) - t76) * qJD(4)) * MDP(14) + (t60 * t49 - t51 * t54 + ((t56 * t73 + t57 * t72) * qJD(4) + t75) * qJD(4)) * MDP(15) + (0.2e1 * t87 + 0.2e1 * t88) * qJD(3); (-MDP(14) * t55 + MDP(15) * t54) * qJD(4); t62 * MDP(15) + ((t68 * t82 + t70 * t83 + t53) * MDP(14) + (-t52 - t80) * MDP(15)) * qJD(4) + (-t87 - t88) * qJD(1); t53 * t52 * MDP(9) + (-t52 ^ 2 + t53 ^ 2) * MDP(10) + (t62 + (t52 - t80) * qJD(4)) * MDP(11) + (-qJD(1) * t76 - t51 * t53) * MDP(14) + (qJD(1) * t75 + t51 * t52) * MDP(15);];
tauc = t1;
