% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:13
% EndTime: 2019-03-08 18:36:14
% DurationCPUTime: 0.24s
% Computational Cost: add. (181->61), mult. (423->97), div. (0->0), fcn. (206->4), ass. (0->37)
t80 = qJD(2) + qJD(3);
t56 = qJD(1) + qJD(2);
t60 = cos(qJ(2));
t77 = pkin(1) * qJD(1);
t70 = t60 * t77;
t51 = t56 * pkin(2) + t70;
t59 = cos(qJ(3));
t57 = sin(qJ(3));
t58 = sin(qJ(2));
t71 = t58 * t77;
t67 = t57 * t71;
t78 = t80 * t67;
t41 = (qJD(2) * t70 + qJD(3) * t51) * t59 - t78;
t79 = t41 * t57;
t76 = t59 * MDP(8);
t75 = t59 * MDP(9);
t74 = qJD(3) * t57;
t73 = qJD(3) * t59;
t69 = t51 * t74;
t68 = MDP(8) * t73;
t66 = -t57 * t60 - t58 * t59;
t65 = -t57 * t58 + t59 * t60;
t64 = -t57 * MDP(8) - t75;
t46 = t59 * t51 - t67;
t63 = -t51 * t73 + t78;
t62 = -MDP(6) + t64;
t61 = (t66 * qJD(2) - t58 * t73) * pkin(1);
t55 = qJD(1) + t80;
t54 = t60 * pkin(1) + pkin(2);
t49 = t65 * t77;
t48 = t66 * t77;
t47 = t57 * t51 + t59 * t71;
t45 = t55 * pkin(3) + t46;
t44 = -t54 * t74 + t61;
t43 = t54 * t73 + (t65 * qJD(2) - t58 * t74) * pkin(1);
t42 = qJD(1) * t61 - t69;
t1 = [(t44 * t55 - t69) * MDP(8) + (-t43 * t55 + t63) * MDP(9) + (t54 * t79 + t47 * t43 + t42 * (t59 * t54 + pkin(3)) + t45 * t44) * MDP(10) + ((-qJD(1) * t68 + (t41 * t59 - t42 * t57) * MDP(10)) * t58 + ((-t58 * MDP(5) - t60 * MDP(6)) * t56 + ((-MDP(5) - t76) * t58 + t62 * t60) * qJD(1)) * qJD(2)) * pkin(1); -t48 * t55 * MDP(8) + (t49 * t55 + t78) * MDP(9) + (pkin(2) * t79 + t42 * (t59 * pkin(2) + pkin(3)) - t47 * t49 - t45 * t48) * MDP(10) + (t64 * t51 + ((-t45 * t57 + t47 * t59) * MDP(10) + t64 * t55) * pkin(2)) * qJD(3) + ((t56 * MDP(6) + t62 * qJD(2)) * t60 + ((-qJD(2) + t56) * MDP(5) - t80 * t76) * t58) * t77; (t47 * t55 - t69) * MDP(8) + (t46 * t55 + t63) * MDP(9) + (t42 * pkin(3) + (t45 - t46) * t47) * MDP(10) + (-t58 * t68 + (t66 * MDP(8) - t60 * t75) * qJD(2)) * t77; 0;];
tauc  = t1;
