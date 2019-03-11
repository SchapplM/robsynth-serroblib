% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:12:24
% EndTime: 2019-03-08 18:12:25
% DurationCPUTime: 0.04s
% Computational Cost: add. (24->17), mult. (58->25), div. (0->0), fcn. (21->2), ass. (0->11)
t22 = qJD(3) * pkin(3);
t16 = sin(qJ(3));
t21 = t16 * qJD(2);
t17 = cos(qJ(3));
t20 = t17 * qJD(2);
t19 = MDP(8) * qJD(3);
t18 = qJD(3) ^ 2;
t15 = qJD(3) * qJ(4) + t21;
t14 = qJD(4) - t20 - t22;
t13 = (qJD(4) + t20) * qJD(3);
t1 = [0; (t15 * t19 + (-MDP(5) + MDP(7)) * t18) * t17 + ((t13 + (t14 - t20) * qJD(3)) * MDP(8) + (-MDP(4) - MDP(6)) * t18) * t16; 0.2e1 * qJD(3) * qJD(4) * MDP(7) + (t13 * qJ(4) + t15 * qJD(4) + (-t15 * t17 + (-t14 - t22) * t16) * qJD(2)) * MDP(8); -t18 * MDP(7) + (-t15 + t21) * t19;];
tauc  = t1;
