% Calculate joint inertia matrix for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPPR6_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:52
% EndTime: 2019-12-31 17:47:53
% DurationCPUTime: 0.27s
% Computational Cost: add. (199->76), mult. (380->103), div. (0->0), fcn. (329->6), ass. (0->37)
t101 = pkin(3) + qJ(2);
t80 = cos(pkin(7));
t76 = t80 ^ 2;
t78 = sin(pkin(7));
t100 = t78 ^ 2 + t76;
t79 = cos(pkin(8));
t98 = t79 * t80;
t77 = sin(pkin(8));
t97 = t77 * t80;
t89 = -t78 * qJ(3) - pkin(1);
t62 = (-pkin(2) - qJ(4)) * t80 + t89;
t66 = t101 * t78;
t57 = t79 * t62 + t77 * t66;
t96 = t100 * qJ(2) ^ 2;
t67 = t101 * t80;
t56 = -t77 * t62 + t79 * t66;
t95 = t56 * MDP(15);
t94 = t57 * MDP(15);
t93 = t77 * MDP(13);
t73 = t77 ^ 2;
t75 = t79 ^ 2;
t92 = MDP(11) + (t75 + t73) * MDP(15);
t81 = sin(qJ(5));
t91 = t81 * t97;
t82 = cos(qJ(5));
t90 = t82 * t97;
t60 = t78 * t82 + t91;
t61 = t78 * t81 - t90;
t88 = t61 * MDP(18) + t60 * MDP(19);
t55 = t78 * pkin(6) + t57;
t58 = (pkin(4) * t79 + pkin(6) * t77) * t80 + t67;
t87 = (-t81 * t55 + t82 * t58) * MDP(21) - (t82 * t55 + t81 * t58) * MDP(22);
t86 = -t60 * MDP(21) + t61 * MDP(22);
t85 = t82 * MDP(21) - t81 * MDP(22);
t84 = -t81 * MDP(21) - t82 * MDP(22);
t64 = -t80 * pkin(2) + t89;
t1 = [MDP(1) + (pkin(1) ^ 2 + t96) * MDP(7) + (t64 ^ 2 + t96) * MDP(11) + (t56 ^ 2 + t57 ^ 2 + t67 ^ 2) * MDP(15) + t75 * t76 * MDP(20) + (MDP(16) * t61 + 0.2e1 * t60 * MDP(17)) * t61 + 0.2e1 * t86 * (-t78 * pkin(4) - t56) + 0.2e1 * (MDP(6) + MDP(8)) * t100 * qJ(2) + 0.2e1 * t88 * t98 + 0.2e1 * (-t64 * MDP(10) + t56 * MDP(12) - t57 * MDP(13) - pkin(1) * MDP(5)) * t78 + 0.2e1 * (pkin(1) * MDP(4) + t64 * MDP(9) + (-t67 * MDP(13) + t56 * MDP(14)) * t77 + (t67 * MDP(12) - t57 * MDP(14) + t87) * t79) * t80; t79 * t94 + t64 * MDP(11) - pkin(1) * MDP(7) + (t86 - t95) * t77 + (-t73 * MDP(14) - MDP(4) + MDP(9) + (-MDP(14) + t84) * t75) * t80 + (-t77 * MDP(12) - t79 * MDP(13) - MDP(10) + MDP(5)) * t78; MDP(7) + t92; t77 * t94 + (qJ(2) * MDP(11) + MDP(8) - t93) * t78 + (t78 * MDP(12) + t95 + (t60 - t91) * MDP(21) + (-t61 - t90) * MDP(22)) * t79; 0; t92; t67 * MDP(15) + (-t93 + (MDP(12) + t85) * t79) * t80; 0; 0; MDP(15); MDP(20) * t98 + t87 + t88; t84 * t79; t84 * t77; t85; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
