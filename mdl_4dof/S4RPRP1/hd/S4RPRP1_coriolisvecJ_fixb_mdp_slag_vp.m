% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRP1
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
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:43
% EndTime: 2019-03-08 18:29:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (127->37), mult. (287->51), div. (0->0), fcn. (132->4), ass. (0->27)
t69 = MDP(6) + MDP(8);
t68 = pkin(1) * sin(pkin(6));
t49 = cos(pkin(6)) * pkin(1) + pkin(2);
t48 = t49 * qJD(1);
t54 = sin(qJ(3));
t67 = t54 * t48;
t62 = qJD(3) * t68;
t61 = qJD(1) * t62;
t55 = cos(qJ(3));
t65 = qJD(3) * t55;
t66 = t48 * t65 - t54 * t61;
t38 = qJD(3) * t67 + t55 * t61;
t63 = qJD(1) * t68;
t39 = t55 * t48 - t54 * t63;
t64 = qJD(4) - t39;
t51 = qJD(1) + qJD(3);
t50 = t51 * qJD(4);
t35 = t50 + t66;
t60 = t39 * t51 - t66;
t57 = t49 * t65 - t54 * t62;
t56 = t54 * t49 + t55 * t68;
t40 = t55 * t63 + t67;
t42 = t56 * qJD(3);
t41 = qJD(4) + t57;
t37 = t51 * qJ(4) + t40;
t36 = -t51 * pkin(3) + t64;
t1 = [(-t57 * t51 - t66) * MDP(7) + (t41 * t51 + t35) * MDP(9) + (t35 * (qJ(4) + t56) + t37 * t41 + t38 * (-t55 * t49 + t54 * t68 - pkin(3)) + t36 * t42) * MDP(10) + t69 * (-t42 * t51 - t38); 0; t60 * MDP(7) + (0.2e1 * t50 - t60) * MDP(9) + (-t38 * pkin(3) + t35 * qJ(4) - t36 * t40 + t64 * t37) * MDP(10) + t69 * (t40 * t51 - t38); -t51 ^ 2 * MDP(9) + (-t37 * t51 + t38) * MDP(10);];
tauc  = t1;
