% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:07
% EndTime: 2019-12-05 14:58:08
% DurationCPUTime: 0.23s
% Computational Cost: add. (147->38), mult. (413->85), div. (0->0), fcn. (320->8), ass. (0->31)
t90 = qJD(4) * pkin(4);
t69 = sin(qJ(5));
t71 = cos(qJ(5));
t89 = (t69 ^ 2 - t71 ^ 2) * MDP(8);
t65 = sin(pkin(9));
t67 = cos(pkin(9));
t70 = sin(qJ(4));
t72 = cos(qJ(4));
t61 = t72 * t65 + t70 * t67;
t57 = t61 * qJD(4);
t66 = sin(pkin(8));
t52 = t66 * t57;
t73 = qJD(5) ^ 2;
t88 = t61 * t73;
t86 = t71 * MDP(7);
t85 = qJD(1) * t66;
t60 = t70 * t65 - t72 * t67;
t56 = t60 * qJD(4);
t84 = t66 * t56 * qJD(4);
t83 = t69 * qJD(5);
t82 = t71 * MDP(13);
t81 = t71 * qJD(5);
t77 = t70 * (t67 * qJD(2) - t65 * t85) + t72 * (t65 * qJD(2) + t67 * t85);
t80 = pkin(6) * t73 + t77 * qJD(4);
t79 = qJD(5) * t90;
t78 = t90 * qJD(4);
t76 = -t69 * MDP(12) - t82;
t75 = 0.2e1 * t52 + qJD(5) * cos(pkin(8));
t74 = qJD(4) ^ 2;
t55 = t60 * t66;
t1 = [MDP(5) * t84 + t52 * qJD(4) * MDP(6) + (t71 * t84 + (t55 * t81 + t69 * t75) * qJD(5)) * MDP(12) + (-t69 * t84 + (-t55 * t83 + t71 * t75) * qJD(5)) * MDP(13); (t56 * t83 - t71 * t88) * MDP(12) + (t56 * t81 + t69 * t88) * MDP(13) + (-t57 * MDP(5) + t56 * MDP(6) + (-t57 * t71 + t60 * t83) * MDP(12) + (t57 * t69 + t60 * t81) * MDP(13)) * qJD(4); t76 * t73; (-t80 * MDP(12) - MDP(13) * t79 + t73 * MDP(9)) * t71 + (-t73 * MDP(10) - MDP(12) * t79 + t80 * MDP(13)) * t69 + ((t71 * MDP(12) - t69 * MDP(13)) * t77 + (t76 * pkin(4) + 0.2e1 * t69 * t86 - 0.2e1 * t89) * qJD(5)) * qJD(4); t74 * t89 + t78 * t82 + (t78 * MDP(12) - t74 * t86) * t69;];
tauc = t1;
