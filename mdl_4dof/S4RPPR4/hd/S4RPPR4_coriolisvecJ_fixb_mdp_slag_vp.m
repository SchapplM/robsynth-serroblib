% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPPR4
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
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:54
% EndTime: 2019-12-31 16:38:54
% DurationCPUTime: 0.08s
% Computational Cost: add. (58->22), mult. (141->37), div. (0->0), fcn. (49->4), ass. (0->14)
t32 = cos(qJ(4));
t37 = t32 * MDP(14);
t31 = sin(qJ(4));
t38 = t31 * MDP(13);
t46 = t37 + t38;
t45 = -MDP(13) * t32 + MDP(14) * t31;
t44 = (t31 ^ 2 - t32 ^ 2) * MDP(9);
t26 = sin(pkin(6)) * pkin(1) + qJ(3);
t24 = qJD(1) * t26;
t39 = t24 * qJD(1);
t36 = t32 * t31 * MDP(8);
t34 = qJD(1) ^ 2;
t33 = qJD(4) ^ 2;
t1 = [(qJD(3) * MDP(7) - qJD(4) * t45) * t24 + (-t31 * MDP(10) - t32 * MDP(11) - t46 * (-cos(pkin(6)) * pkin(1) - pkin(2) - pkin(5))) * t33 + ((t26 * MDP(7) + (2 * MDP(6)) + 0.2e1 * t37 + 0.2e1 * t38) * qJD(3) + (-t26 * t45 - 0.2e1 * t36 + 0.2e1 * t44) * qJD(4)) * qJD(1); t45 * t33; -t34 * MDP(6) - MDP(7) * t39 + t46 * (-t33 - t34); t45 * t39 + (t36 - t44) * t34;];
tauc = t1;
