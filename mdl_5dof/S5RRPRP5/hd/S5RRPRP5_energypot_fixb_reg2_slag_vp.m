% Calculate inertial parameters regressor of potential energy for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:56
% EndTime: 2019-12-31 19:54:56
% DurationCPUTime: 0.07s
% Computational Cost: add. (119->43), mult. (106->50), div. (0->0), fcn. (93->8), ass. (0->25)
t83 = g(3) * pkin(5);
t73 = sin(qJ(2));
t82 = t73 * pkin(2) + pkin(5);
t75 = cos(qJ(2));
t63 = t75 * pkin(2) + pkin(1);
t72 = -pkin(6) - qJ(3);
t71 = qJ(2) + pkin(8);
t65 = cos(t71);
t57 = pkin(3) * t65 + t63;
t70 = -pkin(7) + t72;
t74 = sin(qJ(1));
t76 = cos(qJ(1));
t81 = t74 * t57 + t76 * t70;
t64 = sin(t71);
t80 = pkin(3) * t64 + t82;
t79 = t76 * t57 - t74 * t70;
t78 = g(1) * t76 + g(2) * t74;
t66 = qJ(4) + t71;
t61 = sin(t66);
t62 = cos(t66);
t77 = pkin(4) * t62 + qJ(5) * t61;
t58 = g(1) * t74 - g(2) * t76;
t54 = -g(3) * t61 - t78 * t62;
t53 = -g(3) * t62 + t78 * t61;
t1 = [0, 0, 0, 0, 0, 0, -t78, t58, -g(3), -t83, 0, 0, 0, 0, 0, 0, -g(3) * t73 - t78 * t75, -g(3) * t75 + t78 * t73, -t58, -g(1) * (t76 * pkin(1) + t74 * pkin(6)) - g(2) * (t74 * pkin(1) - t76 * pkin(6)) - t83, 0, 0, 0, 0, 0, 0, -g(3) * t64 - t78 * t65, -g(3) * t65 + t78 * t64, -t58, -g(1) * (t76 * t63 - t74 * t72) - g(2) * (t74 * t63 + t76 * t72) - g(3) * t82, 0, 0, 0, 0, 0, 0, t54, t53, -t58, -g(1) * t79 - g(2) * t81 - g(3) * t80, 0, 0, 0, 0, 0, 0, t54, -t58, -t53, -g(1) * (t77 * t76 + t79) - g(2) * (t77 * t74 + t81) - g(3) * (t61 * pkin(4) - t62 * qJ(5) + t80);];
U_reg = t1;
