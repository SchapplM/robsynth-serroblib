% Calculate minimal parameter regressor of potential energy for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% U_reg [1x12]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PPRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:26
% EndTime: 2019-12-31 16:17:26
% DurationCPUTime: 0.05s
% Computational Cost: add. (19->11), mult. (38->19), div. (0->0), fcn. (42->6), ass. (0->11)
t66 = cos(qJ(3));
t65 = sin(qJ(3));
t64 = g(3) * qJ(1);
t59 = sin(pkin(6));
t60 = cos(pkin(6));
t55 = -t59 * t65 - t60 * t66;
t56 = t59 * t66 - t60 * t65;
t63 = g(1) * t55 - g(2) * t56;
t62 = cos(qJ(4));
t61 = sin(qJ(4));
t1 = [-t64, -g(1) * (t60 * pkin(1) + t59 * qJ(2)) - g(2) * (t59 * pkin(1) - t60 * qJ(2)) - t64, 0, t63, -g(1) * t56 - g(2) * t55, 0, 0, 0, 0, 0, g(3) * t61 + t63 * t62, g(3) * t62 - t63 * t61;];
U_reg = t1;
