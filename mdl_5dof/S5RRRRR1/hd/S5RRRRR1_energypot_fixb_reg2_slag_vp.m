% Calculate inertial parameters regressor of potential energy for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 14:52:09
% EndTime: 2018-11-16 14:52:09
% DurationCPUTime: 0.09s
% Computational Cost: add. (112->33), mult. (111->46), div. (0->0), fcn. (102->10), ass. (0->26)
t59 = sin(qJ(1));
t62 = cos(qJ(1));
t76 = -g(1) * t62 - g(2) * t59;
t75 = g(3) * pkin(5);
t56 = qJ(2) + qJ(3);
t54 = qJ(4) + t56;
t49 = sin(t54);
t72 = g(3) * t49;
t61 = cos(qJ(2));
t71 = t61 * pkin(2) + pkin(1);
t57 = sin(qJ(5));
t70 = t59 * t57;
t60 = cos(qJ(5));
t69 = t59 * t60;
t68 = t62 * t57;
t67 = t62 * t60;
t58 = sin(qJ(2));
t66 = -t58 * pkin(2) + pkin(5);
t52 = sin(t56);
t64 = -pkin(3) * t52 + t66;
t53 = cos(t56);
t50 = cos(t54);
t48 = g(1) * t59 - g(2) * t62;
t47 = pkin(3) * t53 + t71;
t46 = g(3) * t50 - t76 * t49;
t1 = [0, 0, 0, 0, 0, 0, t76, t48, -g(3), -t75, 0, 0, 0, 0, 0, 0, g(3) * t58 + t61 * t76, g(3) * t61 - t58 * t76, t48, pkin(1) * t76 - t75, 0, 0, 0, 0, 0, 0, g(3) * t52 + t53 * t76, g(3) * t53 - t52 * t76, t48, -g(3) * t66 + t71 * t76, 0, 0, 0, 0, 0, 0, t76 * t50 + t72, t46, t48, -g(3) * t64 + t47 * t76, 0, 0, 0, 0, 0, 0, -g(1) * (t50 * t67 - t70) - g(2) * (t50 * t69 + t68) + t60 * t72, -g(1) * (-t50 * t68 - t69) - g(2) * (-t50 * t70 + t67) - t57 * t72, -t46, -g(3) * (-t49 * pkin(4) + t50 * pkin(6) + t64) + t76 * (pkin(4) * t50 + pkin(6) * t49 + t47);];
U_reg  = t1;
