% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(2,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:25
% EndTime: 2019-07-18 13:27:25
% DurationCPUTime: 0.18s
% Computational Cost: add. (85->31), mult. (198->40), div. (0->0), fcn. (90->4), ass. (0->22)
t35 = sin(qJ(4));
t58 = MDP(10) * t35;
t47 = qJD(3) + qJD(4);
t37 = cos(qJ(4));
t53 = t37 * MDP(10);
t57 = t35 * MDP(9) + t53;
t56 = t57 * qJD(4);
t55 = pkin(1) * qJD(2);
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
t29 = t34 * pkin(2) + t38 * t55;
t1 = [0; t43 + (t40 * t36 * qJD(4) + ((t51 * MDP(6) + t40) * t36 + (t51 * MDP(7) + t35 * t42 + t52 * t53) * t38) * qJD(3)) * pkin(1) + (-(t38 * pkin(1) + pkin(2)) * t33 - t29) * t56; t43 + ((t49 * MDP(7) + t57 * (-qJD(3) + t33)) * t38 + (t49 * MDP(6) + t39) * t36) * t55 + (-pkin(2) * t33 - t29) * t56; t43 + (-qJD(3) * t38 * t57 + t39 * t36) * t55 + t57 * t29 * (-qJD(4) + t33);];
tauc  = t1;
