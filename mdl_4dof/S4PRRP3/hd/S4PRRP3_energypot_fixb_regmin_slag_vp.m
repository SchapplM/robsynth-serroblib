% Calculate minimal parameter regressor of potential energy for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% U_reg [1x13]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:55
% EndTime: 2019-12-31 16:26:55
% DurationCPUTime: 0.03s
% Computational Cost: add. (34->18), mult. (29->20), div. (0->0), fcn. (25->6), ass. (0->10)
t62 = pkin(6) + qJ(2);
t60 = sin(t62);
t61 = cos(t62);
t66 = g(1) * t61 + g(2) * t60;
t65 = cos(qJ(3));
t64 = sin(qJ(3));
t63 = -pkin(5) - qJ(4);
t59 = t65 * pkin(3) + pkin(2);
t58 = g(1) * t60 - g(2) * t61;
t1 = [-g(3) * qJ(1), 0, -t66, t58, 0, 0, 0, 0, 0, -g(3) * t64 - t66 * t65, -g(3) * t65 + t66 * t64, -t58, -g(1) * (t61 * t59 - t60 * t63 + cos(pkin(6)) * pkin(1)) - g(2) * (t60 * t59 + t61 * t63 + sin(pkin(6)) * pkin(1)) - g(3) * (t64 * pkin(3) + pkin(4) + qJ(1));];
U_reg = t1;
