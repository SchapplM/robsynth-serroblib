% Calculate inertial parameters regressor of potential energy for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:53
% EndTime: 2019-12-31 17:04:53
% DurationCPUTime: 0.05s
% Computational Cost: add. (67->34), mult. (71->41), div. (0->0), fcn. (61->8), ass. (0->19)
t62 = g(3) * pkin(4);
t56 = sin(qJ(2));
t61 = t56 * pkin(2) + pkin(4);
t58 = cos(qJ(2));
t46 = t58 * pkin(2) + pkin(1);
t55 = -pkin(5) - qJ(3);
t54 = qJ(2) + pkin(7);
t57 = sin(qJ(1));
t59 = cos(qJ(1));
t60 = g(1) * t59 + g(2) * t57;
t53 = -pkin(6) + t55;
t49 = qJ(4) + t54;
t48 = cos(t54);
t47 = sin(t54);
t45 = cos(t49);
t44 = sin(t49);
t43 = g(1) * t57 - g(2) * t59;
t42 = pkin(3) * t48 + t46;
t1 = [0, 0, 0, 0, 0, 0, -t60, t43, -g(3), -t62, 0, 0, 0, 0, 0, 0, -g(3) * t56 - t60 * t58, -g(3) * t58 + t60 * t56, -t43, -g(1) * (t59 * pkin(1) + t57 * pkin(5)) - g(2) * (t57 * pkin(1) - t59 * pkin(5)) - t62, 0, 0, 0, 0, 0, 0, -g(3) * t47 - t60 * t48, -g(3) * t48 + t60 * t47, -t43, -g(1) * (t59 * t46 - t57 * t55) - g(2) * (t57 * t46 + t59 * t55) - g(3) * t61, 0, 0, 0, 0, 0, 0, -g(3) * t44 - t60 * t45, -g(3) * t45 + t60 * t44, -t43, -g(1) * (t59 * t42 - t57 * t53) - g(2) * (t57 * t42 + t59 * t53) - g(3) * (pkin(3) * t47 + t61);];
U_reg = t1;
