% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:48
% EndTime: 2019-12-31 16:42:49
% DurationCPUTime: 0.19s
% Computational Cost: add. (176->51), mult. (423->88), div. (0->0), fcn. (198->4), ass. (0->34)
t55 = sin(pkin(6)) * pkin(1) + pkin(5);
t77 = qJ(4) + t55;
t62 = sin(qJ(3));
t63 = cos(qJ(3));
t69 = t77 * qJD(1);
t46 = t63 * qJD(2) - t69 * t62;
t47 = t62 * qJD(2) + t69 * t63;
t84 = qJD(3) * t47;
t58 = t62 ^ 2;
t59 = t63 ^ 2;
t83 = -t62 * t63 * MDP(5) + (t58 - t59) * MDP(6);
t64 = qJD(3) ^ 2;
t82 = t62 * t64;
t81 = t63 * t64;
t78 = qJD(3) * pkin(3);
t45 = t46 + t78;
t80 = t45 - t46;
t56 = -cos(pkin(6)) * pkin(1) - pkin(2);
t54 = qJD(1) * t56;
t76 = qJD(1) * t63;
t75 = t62 * qJD(1);
t73 = qJD(1) * qJD(4);
t70 = qJD(3) * t77;
t67 = 0.2e1 * qJD(3) * t54;
t66 = (-t63 * pkin(3) + t56) * qJD(1);
t65 = qJD(1) ^ 2;
t52 = t77 * t63;
t51 = t77 * t62;
t50 = qJD(4) + t66;
t49 = -t62 * qJD(4) - t63 * t70;
t48 = t63 * qJD(4) - t62 * t70;
t44 = -t62 * t73 - t84;
t43 = t46 * qJD(3) + t63 * t73;
t1 = [MDP(7) * t81 - MDP(8) * t82 + (-t55 * t81 + t62 * t67) * MDP(10) + (t55 * t82 + t63 * t67) * MDP(11) + (t43 * t63 - t44 * t62 + (-t45 * t63 - t47 * t62) * qJD(3) + (t48 * t63 - t49 * t62 + (t51 * t63 - t52 * t62) * qJD(3)) * qJD(1)) * MDP(12) + (t43 * t52 - t44 * t51 + t45 * t49 + t47 * t48 + (t50 + t66) * t62 * t78) * MDP(13) - 0.2e1 * t83 * qJD(1) * qJD(3); (-t64 * MDP(11) + (t44 + t84) * MDP(13)) * t63 + (-t64 * MDP(10) + (-qJD(3) * t45 + t43) * MDP(13)) * t62; (-t78 + t80) * MDP(12) * t76 + (t80 * t47 + (-t50 * t75 + t44) * pkin(3)) * MDP(13) + t83 * t65 + (-MDP(10) * t75 - MDP(11) * t76) * t54; (-t58 - t59) * MDP(12) * t65 + (-t47 * t63 + (t45 + t78) * t62) * MDP(13) * qJD(1);];
tauc = t1;
