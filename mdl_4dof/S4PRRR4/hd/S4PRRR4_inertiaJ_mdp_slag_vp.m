% Calculate joint inertia matrix for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4PRRR4_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:39
% EndTime: 2019-12-31 16:32:40
% DurationCPUTime: 0.05s
% Computational Cost: add. (54->28), mult. (104->44), div. (0->0), fcn. (93->4), ass. (0->16)
t33 = cos(qJ(3));
t39 = -0.2e1 * t33 * pkin(3) - (2 * pkin(2));
t38 = pkin(5) + pkin(6);
t30 = sin(qJ(4));
t31 = sin(qJ(3));
t32 = cos(qJ(4));
t24 = t30 * t31 - t32 * t33;
t20 = t24 * MDP(17);
t25 = t30 * t33 + t32 * t31;
t37 = -t25 * MDP(18) - t20;
t36 = t33 * MDP(10);
t26 = t38 * t31;
t27 = t38 * t33;
t35 = t25 * MDP(14) - t24 * MDP(15) + (-t32 * t26 - t30 * t27) * MDP(17) + (t30 * t26 - t32 * t27) * MDP(18);
t34 = (MDP(17) * t32 - MDP(18) * t30) * pkin(3);
t1 = [MDP(1); 0; 0.2e1 * pkin(2) * t36 + t20 * t39 + MDP(2) + (-(2 * MDP(11) * pkin(2)) + MDP(5) * t31 + 0.2e1 * MDP(6) * t33) * t31 + (MDP(12) * t25 - 0.2e1 * MDP(13) * t24 + MDP(18) * t39) * t25; -t31 * MDP(11) + t36 + t37; t31 * MDP(7) + t33 * MDP(8) + (-MDP(10) * t31 - MDP(11) * t33) * pkin(5) + t35; MDP(16) + MDP(9) + 0.2e1 * t34; t37; t35; MDP(16) + t34; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
