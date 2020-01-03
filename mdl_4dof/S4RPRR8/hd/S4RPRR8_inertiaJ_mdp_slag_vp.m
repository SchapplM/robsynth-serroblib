% Calculate joint inertia matrix for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S4RPRR8_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:15
% EndTime: 2019-12-31 16:55:15
% DurationCPUTime: 0.05s
% Computational Cost: add. (69->34), mult. (112->50), div. (0->0), fcn. (93->4), ass. (0->16)
t30 = sin(qJ(3));
t39 = 0.2e1 * t30 * pkin(3) + (2 * qJ(2));
t33 = -pkin(1) - pkin(5);
t38 = -pkin(6) + t33;
t37 = (MDP(6) * pkin(1));
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t32 = cos(qJ(3));
t24 = t29 * t32 + t31 * t30;
t25 = -t29 * t30 + t31 * t32;
t36 = t25 * MDP(19) - t24 * MDP(20);
t26 = t38 * t30;
t27 = t38 * t32;
t35 = t25 * MDP(16) - t24 * MDP(17) + (-t29 * t26 + t31 * t27) * MDP(19) + (-t31 * t26 - t29 * t27) * MDP(20);
t34 = (MDP(19) * t31 - MDP(20) * t29) * pkin(3);
t1 = [t24 * MDP(19) * t39 + MDP(1) + (MDP(7) * t32 - 0.2e1 * t30 * MDP(8)) * t32 + ((-2 * MDP(4) + t37) * pkin(1)) + (MDP(14) * t25 - 0.2e1 * t24 * MDP(15) + MDP(20) * t39) * t25 + (0.2e1 * t30 * MDP(12) + 0.2e1 * t32 * MDP(13) + (MDP(6) * qJ(2)) + (2 * MDP(5))) * qJ(2); MDP(4) - t37; MDP(6); (MDP(12) * t33 + MDP(9)) * t32 + (-MDP(13) * t33 - MDP(10)) * t30 + t35; t32 * MDP(12) - t30 * MDP(13) + t36; MDP(11) + MDP(18) + 0.2e1 * t34; t35; t36; MDP(18) + t34; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
