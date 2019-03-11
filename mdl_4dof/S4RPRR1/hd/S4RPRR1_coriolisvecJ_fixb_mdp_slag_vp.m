% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRR1
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
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:50
% EndTime: 2019-03-08 18:31:51
% DurationCPUTime: 0.22s
% Computational Cost: add. (181->48), mult. (440->71), div. (0->0), fcn. (234->6), ass. (0->29)
t48 = cos(pkin(7)) * pkin(1) + pkin(2);
t47 = t48 * qJD(1);
t54 = sin(qJ(3));
t56 = cos(qJ(3));
t66 = pkin(1) * sin(pkin(7));
t63 = qJD(1) * t66;
t39 = t56 * t47 - t54 * t63;
t50 = qJD(1) + qJD(3);
t49 = qJD(4) + t50;
t67 = qJD(4) - t49;
t40 = t54 * t47 + t56 * t63;
t53 = sin(qJ(4));
t65 = t53 * t40;
t55 = cos(qJ(4));
t64 = t55 * t40;
t35 = t50 * pkin(3) + t39;
t62 = -pkin(3) * t49 - t35;
t37 = t39 * qJD(3);
t38 = t40 * qJD(3);
t61 = -t53 * t37 - t55 * t38;
t59 = -t53 * t35 - t64;
t44 = t54 * t48 + t56 * t66;
t58 = t56 * t48 - t54 * t66;
t36 = qJD(4) * t65;
t57 = -t55 * t37 + t53 * t38 + t36;
t43 = pkin(3) + t58;
t42 = t44 * qJD(3);
t41 = t58 * qJD(3);
t1 = [(-t42 * t50 - t38) * MDP(6) + (-t41 * t50 - t37) * MDP(7) + ((-t53 * t41 - t55 * t42) * t49 + t61) * MDP(9) + (-(t55 * t41 - t53 * t42) * t49 + t57) * MDP(10) + (((-t43 * t53 - t44 * t55) * t49 + t59) * MDP(9) + (-(t43 * t55 - t44 * t53) * t49 - t55 * t35) * MDP(10)) * qJD(4); 0; (t40 * t50 - t38) * MDP(6) + (t39 * t50 - t37) * MDP(7) + (-(-t53 * t39 - t64) * t49 + t61) * MDP(9) + ((t55 * t39 - t65) * t49 + t57) * MDP(10) + (t62 * MDP(9) * t53 + (t62 * MDP(10) - t40 * MDP(9)) * t55) * qJD(4); (t67 * t59 + t61) * MDP(9) + (t36 + (-t40 * t49 + t38) * t53 + (-t67 * t35 - t37) * t55) * MDP(10);];
tauc  = t1;
