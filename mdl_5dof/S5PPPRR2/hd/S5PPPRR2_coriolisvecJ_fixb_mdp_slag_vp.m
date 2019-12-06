% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPPRR2
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
%   see S5PPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:38
% EndTime: 2019-12-05 14:59:39
% DurationCPUTime: 0.24s
% Computational Cost: add. (145->42), mult. (381->92), div. (0->0), fcn. (276->8), ass. (0->33)
t99 = qJD(4) * pkin(4);
t73 = qJD(4) ^ 2;
t98 = 2 * qJD(4);
t68 = sin(qJ(5));
t70 = cos(qJ(5));
t97 = t70 * MDP(12) - t68 * MDP(13);
t96 = (t68 ^ 2 - t70 ^ 2) * MDP(8);
t64 = sin(pkin(9));
t69 = sin(qJ(4));
t95 = t64 * t69;
t71 = cos(qJ(4));
t94 = t64 * t71;
t65 = sin(pkin(8));
t66 = cos(pkin(9));
t93 = t65 * t66;
t90 = t70 * MDP(7);
t89 = t73 * MDP(6);
t67 = cos(pkin(8));
t58 = -t67 * t69 + t71 * t93;
t88 = qJD(5) * t58;
t87 = t58 * t73;
t86 = (t67 * t71 + t69 * t93) * qJD(4);
t83 = t70 * MDP(13);
t82 = t73 * t94;
t72 = qJD(5) ^ 2;
t77 = t71 * (qJD(1) * t93 + t64 * qJD(2)) + t69 * (-t67 * qJD(1) + qJD(3));
t81 = pkin(6) * t72 + t77 * qJD(4);
t80 = qJD(5) * t99;
t79 = t95 * t98;
t78 = t99 * qJD(4);
t75 = -t68 * MDP(12) - t83;
t74 = -qJD(5) * t64 * t65 + 0.2e1 * t86;
t1 = [-MDP(5) * t87 + t86 * qJD(4) * MDP(6) + (-t70 * t87 + (t74 * t68 - t70 * t88) * qJD(5)) * MDP(12) + (t68 * t87 + (t68 * t88 + t74 * t70) * qJD(5)) * MDP(13); -MDP(5) * t82 + t89 * t95 + (-t70 * t82 + (t68 * t79 + (t66 * t68 - t70 * t94) * qJD(5)) * qJD(5)) * MDP(12) + (t68 * t82 + (t70 * t79 + (t66 * t70 + t68 * t94) * qJD(5)) * qJD(5)) * MDP(13); (qJD(5) * t75 * t98 - t89) * t71 + (-t73 * MDP(5) - t97 * (t72 + t73)) * t69; (-t81 * MDP(12) - MDP(13) * t80 + t72 * MDP(9)) * t70 + (-t72 * MDP(10) - MDP(12) * t80 + t81 * MDP(13)) * t68 + (t97 * t77 + (t75 * pkin(4) + 0.2e1 * t68 * t90 - 0.2e1 * t96) * qJD(5)) * qJD(4); t73 * t96 + t78 * t83 + (t78 * MDP(12) - t73 * t90) * t68;];
tauc = t1;
