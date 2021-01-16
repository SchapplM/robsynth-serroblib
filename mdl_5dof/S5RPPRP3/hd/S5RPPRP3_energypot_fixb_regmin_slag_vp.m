% Calculate minimal parameter regressor of potential energy for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:04:40
% EndTime: 2021-01-15 17:04:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (64->26), mult. (60->31), div. (0->0), fcn. (51->6), ass. (0->19)
t45 = sin(qJ(4));
t54 = pkin(4) * t45;
t44 = qJ(2) + pkin(5);
t53 = g(3) * t44;
t42 = qJ(1) + pkin(7);
t38 = sin(t42);
t46 = sin(qJ(1));
t52 = t46 * pkin(1) + t38 * pkin(2);
t39 = cos(t42);
t48 = cos(qJ(1));
t51 = t48 * pkin(1) + t39 * pkin(2) + t38 * qJ(3);
t50 = -g(1) * t38 + g(2) * t39;
t49 = -g(1) * t48 - g(2) * t46;
t47 = cos(qJ(4));
t43 = -qJ(5) - pkin(6);
t34 = g(1) * t39 + g(2) * t38;
t33 = g(3) * t45 + t50 * t47;
t32 = -g(3) * t47 + t50 * t45;
t1 = [0, t49, g(1) * t46 - g(2) * t48, t49 * pkin(1) - t53, t34, t50, -g(1) * t51 - g(2) * (-t39 * qJ(3) + t52) - t53, 0, 0, 0, 0, 0, t32, t33, t32, t33, -t34, -g(1) * (t38 * t54 - t39 * t43 + t51) - g(2) * (-t38 * t43 + (-qJ(3) - t54) * t39 + t52) - g(3) * (t47 * pkin(4) + pkin(3) + t44);];
U_reg = t1;
