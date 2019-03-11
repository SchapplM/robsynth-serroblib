% Calculate minimal parameter regressor of potential energy for
% S4PRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
% 
% Output:
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_energypot_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:19:53
% EndTime: 2019-03-08 18:19:53
% DurationCPUTime: 0.02s
% Computational Cost: add. (22->16), mult. (29->17), div. (0->0), fcn. (22->2), ass. (0->7)
t39 = sin(qJ(2));
t40 = cos(qJ(2));
t42 = t40 * pkin(2) + t39 * qJ(3) + pkin(1);
t41 = t39 * pkin(2) - t40 * qJ(3) + qJ(1);
t34 = -g(1) * t40 - g(2) * t39;
t33 = g(1) * t39 - g(2) * t40;
t1 = [-g(2) * qJ(1), 0, t34, t33, t34, -t33, -g(3) * pkin(4) - g(1) * t42 - g(2) * t41, t34, -t33, -g(1) * (t40 * pkin(3) + t42) - g(2) * (t39 * pkin(3) + t41) - g(3) * (-qJ(4) + pkin(4));];
U_reg  = t1;
