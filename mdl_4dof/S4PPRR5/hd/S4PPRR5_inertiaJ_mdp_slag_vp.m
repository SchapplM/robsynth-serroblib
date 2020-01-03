% Calculate joint inertia matrix for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PPRR5_inertiaJ_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:52
% EndTime: 2019-12-31 16:19:52
% DurationCPUTime: 0.02s
% Computational Cost: add. (16->12), mult. (37->21), div. (0->0), fcn. (29->4), ass. (0->8)
t12 = cos(qJ(4));
t16 = MDP(11) * t12;
t10 = sin(qJ(4));
t15 = -MDP(11) * t10 - MDP(12) * t12;
t14 = -MDP(12) * t10 + MDP(4) + t16;
t13 = cos(qJ(3));
t11 = sin(qJ(3));
t1 = [MDP(1) + MDP(2); 0; MDP(2); -t13 * MDP(5) - t14 * t11; -t11 * MDP(5) + t14 * t13; 0.2e1 * pkin(3) * t16 + MDP(3) + (-0.2e1 * MDP(12) * pkin(3) + MDP(6) * t10 + 0.2e1 * MDP(7) * t12) * t10; t15 * t13; t15 * t11; t10 * MDP(8) + t12 * MDP(9) + t15 * pkin(5); MDP(10);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
