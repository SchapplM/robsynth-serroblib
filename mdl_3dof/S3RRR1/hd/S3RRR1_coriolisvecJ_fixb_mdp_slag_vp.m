% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [3x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S3RRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:08
% EndTime: 2019-03-08 18:08:09
% DurationCPUTime: 0.16s
% Computational Cost: add. (85->31), mult. (198->40), div. (0->0), fcn. (90->4), ass. (0->22)
t35 = sin(qJ(3));
t58 = MDP(9) * t35;
t47 = qJD(2) + qJD(3);
t37 = cos(qJ(3));
t53 = t37 * MDP(9);
t57 = t35 * MDP(8) + t53;
t56 = t57 * qJD(3);
t55 = pkin(1) * qJD(1);
t33 = qJD(1) + t47;
t52 = -qJD(1) - t33;
t34 = qJD(1) + qJD(2);
t51 = -qJD(1) - t34;
t49 = -qJD(2) + t34;
t46 = t33 * t58;
t36 = sin(qJ(2));
t43 = t47 * t36 * t55 * t58;
t42 = t52 * MDP(8);
t40 = t37 * t42 + t46;
t39 = -t46 + (t33 - t47) * MDP(8) * t37;
t38 = cos(qJ(2));
t29 = t34 * pkin(2) + t38 * t55;
t1 = [t43 + (t40 * t36 * qJD(3) + ((t51 * MDP(5) + t40) * t36 + (t51 * MDP(6) + t35 * t42 + t52 * t53) * t38) * qJD(2)) * pkin(1) + (-(t38 * pkin(1) + pkin(2)) * t33 - t29) * t56; t43 + ((t49 * MDP(6) + t57 * (-qJD(2) + t33)) * t38 + (t49 * MDP(5) + t39) * t36) * t55 + (-pkin(2) * t33 - t29) * t56; t43 + (-qJD(2) * t38 * t57 + t39 * t36) * t55 + t57 * t29 * (-qJD(3) + t33);];
tauc  = t1;
