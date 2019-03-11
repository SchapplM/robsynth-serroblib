% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRRR1
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
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:20
% EndTime: 2019-03-08 18:25:21
% DurationCPUTime: 0.16s
% Computational Cost: add. (85->31), mult. (198->40), div. (0->0), fcn. (90->4), ass. (0->22)
t35 = sin(qJ(4));
t58 = MDP(10) * t35;
t47 = qJD(3) + qJD(4);
t37 = cos(qJ(4));
t53 = t37 * MDP(10);
t57 = t35 * MDP(9) + t53;
t56 = t57 * qJD(4);
t55 = pkin(2) * qJD(2);
t33 = qJD(2) + t47;
t52 = -qJD(2) - t33;
t34 = qJD(2) + qJD(3);
t51 = -qJD(2) - t34;
t49 = -qJD(3) + t34;
t46 = t33 * t58;
t36 = sin(qJ(3));
t43 = t47 * t36 * t55 * t58;
t42 = t52 * MDP(9);
t40 = t37 * t42 + t46;
t39 = -t46 + (t33 - t47) * MDP(9) * t37;
t38 = cos(qJ(3));
t29 = t34 * pkin(3) + t38 * t55;
t1 = [0; t43 + (t40 * t36 * qJD(4) + ((t51 * MDP(6) + t40) * t36 + (t51 * MDP(7) + t35 * t42 + t52 * t53) * t38) * qJD(3)) * pkin(2) + (-(t38 * pkin(2) + pkin(3)) * t33 - t29) * t56; t43 + ((t49 * MDP(7) + t57 * (-qJD(3) + t33)) * t38 + (t49 * MDP(6) + t39) * t36) * t55 + (-pkin(3) * t33 - t29) * t56; t43 + (-t57 * t38 * qJD(3) + t39 * t36) * t55 + t57 * t29 * (-qJD(4) + t33);];
tauc  = t1;
