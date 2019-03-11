% Calculate minimal parameter regressor of potential energy for
% S4PPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta1]';
% 
% Output:
% U_reg [1x6]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PPPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:09:04
% EndTime: 2019-03-08 18:09:04
% DurationCPUTime: 0.04s
% Computational Cost: add. (18->14), mult. (29->20), div. (0->0), fcn. (26->4), ass. (0->10)
t43 = g(3) * qJ(1);
t37 = sin(pkin(5));
t38 = cos(pkin(5));
t42 = t38 * pkin(1) + t37 * qJ(2);
t41 = t37 * pkin(1) - t38 * qJ(2);
t40 = cos(qJ(4));
t39 = sin(qJ(4));
t33 = t37 * t40 + t38 * t39;
t32 = -t37 * t39 + t38 * t40;
t1 = [-t43, -g(1) * t42 - g(2) * t41 - t43, -g(1) * (t38 * qJ(3) + t42) - g(2) * (t37 * qJ(3) + t41) - g(3) * (pkin(2) + qJ(1)) 0, -g(1) * t33 + g(2) * t32, -g(1) * t32 - g(2) * t33;];
U_reg  = t1;
