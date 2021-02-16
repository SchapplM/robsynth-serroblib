% Calculate minimal parameter regressor of potential energy for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:27:39
% EndTime: 2021-01-15 10:27:39
% DurationCPUTime: 0.03s
% Computational Cost: add. (29->19), mult. (47->24), div. (0->0), fcn. (41->4), ass. (0->11)
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t19 = g(1) * t24 - g(2) * t26;
t25 = cos(qJ(3));
t23 = sin(qJ(3));
t22 = pkin(1) + pkin(5) + qJ(4);
t21 = t23 * pkin(3) + qJ(2);
t20 = g(1) * t26 + g(2) * t24;
t18 = g(3) * t23 - t19 * t25;
t17 = -g(3) * t25 - t19 * t23;
t1 = [0, -t20, t19, t20, -t19, -g(1) * (t26 * pkin(1) + t24 * qJ(2)) - g(2) * (t24 * pkin(1) - t26 * qJ(2)) - g(3) * pkin(4), 0, 0, 0, 0, 0, t17, t18, t17, t18, -t20, -g(1) * (t21 * t24 + t22 * t26) - g(2) * (-t21 * t26 + t22 * t24) - g(3) * (t25 * pkin(3) + pkin(2) + pkin(4));];
U_reg = t1;
