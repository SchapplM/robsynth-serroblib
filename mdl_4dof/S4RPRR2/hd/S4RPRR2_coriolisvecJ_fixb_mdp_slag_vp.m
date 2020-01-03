% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:13
% EndTime: 2019-12-31 16:48:14
% DurationCPUTime: 0.17s
% Computational Cost: add. (174->44), mult. (407->68), div. (0->0), fcn. (198->6), ass. (0->35)
t70 = cos(pkin(7)) * pkin(1) + pkin(2);
t69 = t70 * qJD(1);
t79 = sin(qJ(3));
t81 = cos(qJ(3));
t104 = pkin(1) * sin(pkin(7));
t93 = qJD(1) * t104;
t60 = t81 * t69 - t79 * t93;
t78 = sin(qJ(4));
t80 = cos(qJ(4));
t107 = (t78 ^ 2 - t80 ^ 2) * MDP(9);
t61 = t79 * t69 + t81 * t93;
t106 = t61 * qJD(3);
t105 = t78 * MDP(13);
t73 = qJD(1) + qJD(3);
t103 = t73 * pkin(3);
t102 = t61 * t73;
t86 = t81 * t104 + t79 * t70;
t101 = t86 * qJD(3) * t73;
t98 = MDP(8) * t80;
t56 = -t60 - t103;
t97 = qJD(4) * t56;
t96 = t80 * MDP(14);
t82 = qJD(4) ^ 2;
t95 = t82 * MDP(11);
t58 = t60 * qJD(3);
t92 = -t56 * t73 - t58;
t90 = pkin(6) * t82 - t102;
t89 = (pkin(6) + t86) * t82 + t101;
t88 = qJD(4) * (t60 - t103);
t85 = -t79 * t104 + t81 * t70;
t63 = t85 * qJD(3);
t87 = qJD(4) * ((-pkin(3) - t85) * t73 - t63);
t83 = t82 * t80 * MDP(10) + t97 * t105 + (t106 * t78 + t80 * t97) * MDP(14) + 0.2e1 * (t98 * t78 - t107) * qJD(4) * t73;
t72 = t73 ^ 2;
t1 = [(-t106 - t101) * MDP(6) + (-t63 * t73 - t58) * MDP(7) + ((-t106 - t89) * MDP(13) + MDP(14) * t87) * t80 + (MDP(13) * t87 + MDP(14) * t89 - t95) * t78 + t83; (-t96 - t105) * t82; (-t106 + t102) * MDP(6) + (t60 * t73 - t58) * MDP(7) + ((-t106 - t90) * MDP(13) + MDP(14) * t88) * t80 + (MDP(13) * t88 + MDP(14) * t90 - t95) * t78 + t83; t92 * t96 + t72 * t107 + (t92 * MDP(13) - t72 * t98) * t78;];
tauc = t1;
