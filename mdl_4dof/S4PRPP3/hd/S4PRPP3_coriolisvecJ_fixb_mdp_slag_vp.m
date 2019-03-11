% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRPP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:19:54
% EndTime: 2019-03-08 18:19:54
% DurationCPUTime: 0.06s
% Computational Cost: add. (52->27), mult. (112->32), div. (0->0), fcn. (40->2), ass. (0->16)
t39 = MDP(7) + MDP(10);
t34 = sin(qJ(2));
t42 = t34 * qJD(1);
t33 = qJD(2) * qJ(3) + t42;
t35 = cos(qJ(2));
t46 = t33 * t35;
t43 = qJD(2) * pkin(2);
t41 = t35 * qJD(1);
t40 = MDP(6) + MDP(9);
t38 = qJD(2) * (-pkin(2) - pkin(3));
t37 = qJD(3) - t41;
t36 = qJD(2) ^ 2;
t32 = t37 - t43;
t31 = (qJD(3) + t41) * qJD(2);
t28 = t38 + t37;
t1 = [((t32 - t41) * MDP(7) + (t28 - t41) * MDP(10)) * t34 * qJD(2) + ((-MDP(4) + t40) * t35 + (-MDP(3) - MDP(5) - MDP(8)) * t34) * t36 + t39 * (qJD(2) * t46 + t31 * t34); 0.2e1 * t40 * qJD(3) * qJD(2) + (-t39 * t46 + ((-t32 - t43) * MDP(7) + (t38 - t28) * MDP(10)) * t34) * qJD(1) + t39 * (t31 * qJ(3) + t33 * qJD(3)); -t40 * t36 + t39 * (-t33 + t42) * qJD(2); 0;];
tauc  = t1;
