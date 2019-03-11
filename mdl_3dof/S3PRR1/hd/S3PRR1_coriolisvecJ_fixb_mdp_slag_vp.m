% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% MDP [7x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3PRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [3x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S3PRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [7 1]), ...
  'S3PRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [7x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:04:00
% EndTime: 2019-03-08 18:04:01
% DurationCPUTime: 0.09s
% Computational Cost: add. (47->19), mult. (110->33), div. (0->0), fcn. (68->4), ass. (0->15)
t25 = qJD(2) + qJD(3);
t26 = sin(qJ(3));
t27 = sin(qJ(2));
t42 = t26 * t27;
t45 = t25 * MDP(7) * t42;
t28 = cos(qJ(3));
t43 = t26 * MDP(6) + t28 * MDP(7);
t41 = t27 * t28;
t29 = cos(qJ(2));
t38 = qJD(2) * t29;
t34 = qJD(1) * t45;
t33 = -t26 * t29 - t41;
t32 = t28 * t29 - t42;
t24 = qJD(2) * pkin(2) + t29 * qJD(1);
t1 = [(-t27 * MDP(3) - t29 * MDP(4)) * qJD(2) ^ 2 + t25 ^ 2 * (t33 * MDP(6) - t32 * MDP(7)); t34 + ((-t26 * t38 + (-t33 - t41) * t25) * MDP(6) + (t32 * t25 - t28 * t38) * MDP(7)) * qJD(1) + t43 * qJD(3) * (-pkin(2) * t25 - t24); t34 + (-t43 * t38 - t45) * qJD(1) + t43 * t24 * (-qJD(3) + t25);];
tauc  = t1;
