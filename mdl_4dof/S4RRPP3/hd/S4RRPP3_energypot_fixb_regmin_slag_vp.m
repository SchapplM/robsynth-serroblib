% Calculate minimal parameter regressor of potential energy for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:36:03
% EndTime: 2021-01-15 10:36:03
% DurationCPUTime: 0.04s
% Computational Cost: add. (58->29), mult. (70->34), div. (0->0), fcn. (64->8), ass. (0->20)
t53 = sin(qJ(2));
t58 = t53 * pkin(2) + pkin(4);
t54 = sin(qJ(1));
t56 = cos(qJ(1));
t57 = g(1) * t56 + g(2) * t54;
t55 = cos(qJ(2));
t52 = pkin(5) + qJ(3);
t51 = cos(pkin(6));
t50 = sin(pkin(6));
t49 = qJ(2) + pkin(6);
t46 = cos(t49);
t45 = sin(t49);
t44 = t55 * pkin(2) + pkin(1);
t43 = t52 * t54;
t42 = t56 * t52;
t41 = g(1) * t54 - g(2) * t56;
t40 = -g(3) * t45 - t57 * t46;
t39 = -g(3) * t46 + t57 * t45;
t38 = (pkin(3) * t51 + qJ(4) * t50 + pkin(2)) * t55 + (-t50 * pkin(3) + qJ(4) * t51) * t53 + pkin(1);
t1 = [0, -t57, t41, 0, 0, 0, 0, 0, -g(3) * t53 - t57 * t55, -g(3) * t55 + t57 * t53, t40, t39, -t41, -g(1) * (t56 * t44 + t43) - g(2) * (t54 * t44 - t42) - g(3) * t58, t40, -t41, -t39, -g(1) * (t38 * t56 + t43) - g(2) * (t38 * t54 - t42) - g(3) * (t45 * pkin(3) - t46 * qJ(4) + t58);];
U_reg = t1;
