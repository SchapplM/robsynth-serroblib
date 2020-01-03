% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:27
% EndTime: 2019-12-31 16:17:27
% DurationCPUTime: 0.10s
% Computational Cost: add. (42->16), mult. (131->38), div. (0->0), fcn. (58->4), ass. (0->14)
t49 = -2 * qJD(3);
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t48 = -t31 * MDP(11) + t29 * MDP(12);
t47 = (t29 ^ 2 - t31 ^ 2) * MDP(7);
t44 = qJD(3) * pkin(3);
t43 = t31 * MDP(6);
t42 = MDP(12) * t31;
t37 = t44 * qJD(3);
t36 = MDP(11) * t29 + t42;
t35 = -2 * t44;
t34 = qJD(3) ^ 2;
t33 = qJD(4) ^ 2;
t1 = [t36 * t33; (qJD(4) * t36 * t49 - (t34 * MDP(5))) * cos(qJ(3)) + (-(t34 * MDP(4)) + t48 * (t33 + t34)) * sin(qJ(3)); (t31 * MDP(8) - t29 * MDP(9) + t48 * pkin(5)) * t33 + (t47 * t49 + t35 * t42 + (MDP(11) * t35 + 0.2e1 * qJD(3) * t43) * t29) * qJD(4); t34 * t47 + t37 * t42 + (MDP(11) * t37 - t34 * t43) * t29;];
tauc = t1;
