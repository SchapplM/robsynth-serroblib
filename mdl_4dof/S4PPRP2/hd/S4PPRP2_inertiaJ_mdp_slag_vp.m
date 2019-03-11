% Calculate joint inertia matrix for
% S4PPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRP2_inertiaJ_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:13:20
% EndTime: 2019-03-08 18:13:20
% DurationCPUTime: 0.03s
% Computational Cost: add. (24->16), mult. (40->22), div. (0->0), fcn. (38->4), ass. (0->8)
t15 = cos(qJ(3));
t14 = -pkin(3) * MDP(8) - MDP(6);
t13 = sin(qJ(3));
t12 = cos(pkin(5));
t11 = sin(pkin(5));
t9 = t15 * t11 + t13 * t12;
t8 = t13 * t11 - t15 * t12;
t1 = [MDP(1) + (t11 ^ 2 + t12 ^ 2) * MDP(2) + (t8 ^ 2 + t9 ^ 2) * MDP(8); 0; MDP(2) + MDP(8); (MDP(8) * qJ(4) - MDP(5) + MDP(7)) * t9 + (-MDP(4) + t14) * t8; 0; MDP(3) + (2 * pkin(3) * MDP(6)) + 0.2e1 * qJ(4) * MDP(7) + ((pkin(3) ^ 2) + qJ(4) ^ 2) * MDP(8); t8 * MDP(8); 0; t14; MDP(8);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
