% Calculate joint inertia matrix for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPPR6_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:44
% EndTime: 2019-12-31 16:40:45
% DurationCPUTime: 0.08s
% Computational Cost: add. (60->36), mult. (114->51), div. (0->0), fcn. (82->4), ass. (0->18)
t47 = (MDP(7) + MDP(11));
t35 = cos(pkin(6));
t34 = sin(pkin(6));
t39 = t34 * qJ(3) + pkin(1);
t46 = 0.2e1 * (pkin(2) + pkin(3)) * t35 + 0.2e1 * t39;
t45 = -0.2e1 * t34;
t44 = pkin(1) * MDP(7);
t43 = -pkin(5) + qJ(2);
t36 = sin(qJ(4));
t37 = cos(qJ(4));
t24 = t34 * t36 + t35 * t37;
t41 = t24 * MDP(17);
t26 = -t35 * pkin(2) - t39;
t40 = t26 * MDP(11);
t29 = t43 * t35;
t28 = t43 * t34;
t25 = t34 * t37 - t35 * t36;
t1 = [MDP(1) + t41 * t46 + (MDP(10) * t45 - 0.2e1 * t35 * MDP(8) + t40) * t26 + (0.2e1 * t35 * MDP(4) + MDP(5) * t45 + t44) * pkin(1) + (MDP(12) * t25 - 0.2e1 * t24 * MDP(13) + MDP(18) * t46) * t25 + (t47 * qJ(2) + 2 * MDP(6) + 2 * MDP(9)) * (t34 ^ 2 + t35 ^ 2) * qJ(2); t40 - t41 - t25 * MDP(18) - t44 + (-MDP(4) - MDP(8)) * t35 + (-MDP(10) + MDP(5)) * t34; t47; (MDP(11) * qJ(2) + MDP(9)) * t34; 0; MDP(11); t25 * MDP(14) - t24 * MDP(15) + (t37 * t28 - t36 * t29) * MDP(17) + (-t36 * t28 - t37 * t29) * MDP(18); 0; t37 * MDP(17) - t36 * MDP(18); MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
