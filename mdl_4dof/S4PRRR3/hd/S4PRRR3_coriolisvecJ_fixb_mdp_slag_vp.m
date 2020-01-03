% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:39
% EndTime: 2019-12-31 16:31:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (102->38), mult. (231->66), div. (0->0), fcn. (90->4), ass. (0->32)
t53 = sin(qJ(4));
t55 = cos(qJ(4));
t79 = (t53 ^ 2 - t55 ^ 2) * MDP(9);
t78 = t55 * MDP(13) + MDP(6);
t68 = t55 * MDP(14);
t60 = t53 * MDP(13) + t68;
t77 = t60 * qJD(4);
t57 = qJD(4) ^ 2;
t76 = t53 * t57;
t75 = t55 * t57;
t50 = qJD(2) + qJD(3);
t56 = cos(qJ(3));
t72 = pkin(2) * qJD(2);
t44 = -t50 * pkin(3) - t56 * t72;
t54 = sin(qJ(3));
t71 = pkin(2) * qJD(3);
t61 = qJD(2) * t71;
t67 = t55 * qJD(4);
t74 = t54 * t53 * t61 + t44 * t67;
t70 = t53 * qJD(4);
t66 = -qJD(2) - t50;
t65 = -qJD(3) + t50;
t64 = t50 * t70;
t63 = t50 * t67;
t62 = t53 * t50 * MDP(14);
t59 = 0.2e1 * t53 * MDP(8) * t63 - 0.2e1 * t50 * qJD(4) * t79 + MDP(10) * t75 - MDP(11) * t76;
t58 = -t44 * t50 - t56 * t61;
t49 = t50 ^ 2;
t48 = -t56 * pkin(2) - pkin(3);
t47 = t54 * pkin(2) + pkin(6);
t40 = t44 * t70;
t1 = [-t60 * t57; (-t47 * t75 + t48 * t64 + t40) * MDP(13) + (t47 * t76 + t48 * t63 + t74) * MDP(14) + ((t66 * MDP(7) - t77) * t56 + (t78 * t66 + t62) * t54) * t71 + t59; (-pkin(3) * t64 - pkin(6) * t75 + t40) * MDP(13) + (-pkin(3) * t63 + pkin(6) * t76 + t74) * MDP(14) + ((t65 * MDP(7) + t77) * t56 + (t78 * t65 - t62) * t54) * t72 + t59; t58 * t68 + t49 * t79 + (-t49 * t55 * MDP(8) + t58 * MDP(13)) * t53;];
tauc = t1;
