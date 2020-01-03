% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:03
% EndTime: 2019-12-31 17:51:04
% DurationCPUTime: 0.33s
% Computational Cost: add. (287->68), mult. (542->104), div. (0->0), fcn. (231->4), ass. (0->41)
t66 = -cos(pkin(7)) * pkin(1) - pkin(2) - pkin(6);
t96 = qJ(5) - t66;
t104 = -t96 * qJD(1) + qJD(3);
t75 = sin(qJ(4));
t76 = cos(qJ(4));
t57 = t76 * qJD(2) + t104 * t75;
t106 = qJD(4) * t57;
t56 = -t75 * qJD(2) + t104 * t76;
t88 = t76 * qJD(4);
t105 = pkin(4) * t88 + qJD(3);
t71 = t75 ^ 2;
t72 = t76 ^ 2;
t103 = (t71 - t72) * MDP(9);
t86 = qJD(1) * qJD(5);
t54 = -t76 * t86 - t106;
t102 = (t54 + t106) * MDP(16);
t97 = qJD(4) * pkin(4);
t55 = t56 + t97;
t101 = t55 - t56;
t100 = t105 * qJD(1);
t77 = qJD(4) ^ 2;
t78 = qJD(1) ^ 2;
t98 = -t77 - t78;
t95 = qJD(1) * t75;
t94 = t75 * MDP(13);
t92 = t75 * qJD(4);
t91 = t76 * MDP(14);
t90 = t76 * qJD(1);
t84 = t76 * t75 * MDP(8);
t68 = sin(pkin(7)) * pkin(1) + qJ(3);
t63 = t96 * t76;
t53 = t56 * qJD(4) - t75 * t86;
t81 = (-qJD(4) * t55 + t53) * MDP(16);
t80 = t76 * MDP(13) - t75 * MDP(14);
t79 = t75 * pkin(4) + t68;
t65 = qJD(1) * t68;
t62 = t96 * t75;
t61 = t79 * qJD(1) + qJD(5);
t59 = -qJD(4) * t63 - t75 * qJD(5);
t58 = -t76 * qJD(5) + t96 * t92;
t1 = [(-t53 * t75 - t54 * t76 + t55 * t92 - t57 * t88) * MDP(15) + (t100 * t79 + t61 * t105 - t53 * t62 - t54 * t63 + t55 * t58 + t57 * t59) * MDP(16) + (qJD(3) * MDP(7) + t80 * qJD(4)) * t65 + (-t75 * MDP(10) - t76 * MDP(11) + (-t91 - t94) * t66) * t77 + ((-t58 * t76 - t59 * t75) * MDP(15) + (t68 * MDP(7) + (2 * MDP(6)) + 0.2e1 * t91 + 0.2e1 * t94) * qJD(3) + (-0.2e1 * t84 + 0.2e1 * t103 + (t62 * t76 - t63 * t75) * MDP(15) + t80 * t68) * qJD(4)) * qJD(1); (-t77 * MDP(13) + t81) * t76 + (t77 * MDP(14) - t102) * t75; -t78 * MDP(6) + (-t61 * MDP(16) - t65 * MDP(7)) * qJD(1) + (t98 * MDP(14) + t102) * t76 + (t98 * MDP(13) + t81) * t75; (t97 - t101) * MDP(15) * t95 + (t101 * t57 + (-t61 * t90 + t54) * pkin(4)) * MDP(16) + (t84 - t103) * t78 + (-MDP(13) * t90 + t95 * MDP(14)) * t65; ((t55 * t76 + t57 * t75) * qJD(1) + t100) * MDP(16) + (-t71 - t72) * MDP(15) * t78;];
tauc = t1;
