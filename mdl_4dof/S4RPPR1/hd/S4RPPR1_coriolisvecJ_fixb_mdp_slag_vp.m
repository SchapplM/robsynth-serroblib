% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:31
% EndTime: 2019-03-08 18:27:31
% DurationCPUTime: 0.08s
% Computational Cost: add. (69->25), mult. (120->42), div. (0->0), fcn. (45->4), ass. (0->13)
t20 = sin(pkin(6)) * pkin(1) + qJ(3);
t18 = qJD(1) * t20;
t31 = t18 * MDP(7);
t30 = qJD(3) * qJD(1);
t29 = qJD(1) - qJD(4);
t19 = -cos(pkin(6)) * pkin(1) - pkin(2) - pkin(3);
t28 = qJD(3) * (qJD(1) + t29);
t17 = t19 * qJD(1) + qJD(3);
t24 = sin(qJ(4));
t25 = cos(qJ(4));
t27 = t25 * t17 - t24 * t18;
t26 = -t24 * t17 - t25 * t18;
t1 = [0.2e1 * MDP(6) * t30 + 0.2e1 * qJD(3) * t31 + (t24 * t28 + (-(-t19 * t24 - t20 * t25) * t29 - t26) * qJD(4)) * MDP(9) + (t25 * t28 + ((t19 * t25 - t20 * t24) * t29 + t27) * qJD(4)) * MDP(10); 0; (-MDP(6) * qJD(1) - t31) * qJD(1) - (t25 * MDP(10) + t24 * MDP(9)) * t29 ^ 2; (-t24 * t30 + t26 * t29) * MDP(9) + (-t25 * t30 - t27 * t29) * MDP(10) + (-t27 * MDP(10) + t26 * MDP(9)) * qJD(4);];
tauc  = t1;
