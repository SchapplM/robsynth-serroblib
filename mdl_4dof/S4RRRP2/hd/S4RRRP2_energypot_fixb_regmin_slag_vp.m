% Calculate minimal parameter regressor of potential energy for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:04
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:04:28
% EndTime: 2021-01-15 11:04:28
% DurationCPUTime: 0.03s
% Computational Cost: add. (43->19), mult. (42->23), div. (0->0), fcn. (39->6), ass. (0->14)
t29 = qJ(1) + qJ(2);
t27 = sin(t29);
t28 = cos(t29);
t35 = g(1) * t28 + g(2) * t27;
t34 = cos(qJ(1));
t33 = cos(qJ(3));
t32 = sin(qJ(1));
t31 = sin(qJ(3));
t30 = -qJ(4) - pkin(6);
t26 = t33 * pkin(3) + pkin(2);
t25 = g(1) * t27 - g(2) * t28;
t24 = -g(3) * t31 - t35 * t33;
t23 = -g(3) * t33 + t35 * t31;
t1 = [0, -g(1) * t34 - g(2) * t32, g(1) * t32 - g(2) * t34, 0, -t35, t25, 0, 0, 0, 0, 0, t24, t23, t24, t23, -t25, -g(1) * (t34 * pkin(1) + t28 * t26 - t27 * t30) - g(2) * (t32 * pkin(1) + t27 * t26 + t28 * t30) - g(3) * (t31 * pkin(3) + pkin(4) + pkin(5));];
U_reg = t1;
