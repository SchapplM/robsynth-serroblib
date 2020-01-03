% Calculate minimal parameter regressor of potential energy for
% S4PPRR5
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PPRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:51
% EndTime: 2019-12-31 16:19:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (17->16), mult. (38->29), div. (0->0), fcn. (38->6), ass. (0->12)
t65 = cos(qJ(3));
t70 = g(3) * t65;
t69 = g(3) * qJ(1);
t62 = sin(qJ(4));
t63 = sin(qJ(3));
t68 = t62 * t63;
t64 = cos(qJ(4));
t67 = t63 * t64;
t60 = sin(pkin(6));
t61 = cos(pkin(6));
t66 = -g(1) * t60 + g(2) * t61;
t1 = [-t69, -g(1) * (t61 * pkin(1) + t60 * qJ(2)) - g(2) * (t60 * pkin(1) - t61 * qJ(2)) - t69, 0, t63 * t66 - t70, g(3) * t63 + t65 * t66, 0, 0, 0, 0, 0, -g(1) * (t60 * t67 + t61 * t62) - g(2) * (t60 * t62 - t61 * t67) - t64 * t70, -g(1) * (-t60 * t68 + t61 * t64) - g(2) * (t60 * t64 + t61 * t68) + t62 * t70;];
U_reg = t1;
