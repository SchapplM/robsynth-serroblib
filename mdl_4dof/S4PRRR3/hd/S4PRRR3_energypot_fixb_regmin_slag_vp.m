% Calculate minimal parameter regressor of potential energy for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x14]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:38
% EndTime: 2019-12-31 16:31:38
% DurationCPUTime: 0.05s
% Computational Cost: add. (29->10), mult. (19->13), div. (0->0), fcn. (18->6), ass. (0->10)
t62 = pkin(7) + qJ(2);
t61 = qJ(3) + t62;
t57 = sin(t61);
t58 = cos(t61);
t65 = g(1) * t58 + g(2) * t57;
t64 = cos(qJ(4));
t63 = sin(qJ(4));
t60 = cos(t62);
t59 = sin(t62);
t1 = [-g(3) * qJ(1), 0, -g(1) * t60 - g(2) * t59, g(1) * t59 - g(2) * t60, 0, -t65, g(1) * t57 - g(2) * t58, 0, 0, 0, 0, 0, -g(3) * t63 - t65 * t64, -g(3) * t64 + t65 * t63;];
U_reg = t1;
