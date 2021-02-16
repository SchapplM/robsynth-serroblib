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
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:27
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
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
% StartTime: 2021-01-14 22:27:04
% EndTime: 2021-01-14 22:27:04
% DurationCPUTime: 0.03s
% Computational Cost: add. (42->18), mult. (39->20), div. (0->0), fcn. (35->6), ass. (0->12)
t64 = pkin(6) + qJ(2);
t62 = sin(t64);
t63 = cos(t64);
t68 = g(1) * t63 + g(2) * t62;
t67 = cos(qJ(3));
t66 = sin(qJ(3));
t65 = -pkin(5) - qJ(4);
t61 = t67 * pkin(3) + pkin(2);
t60 = g(1) * t62 - g(2) * t63;
t59 = -g(3) * t66 - t68 * t67;
t58 = -g(3) * t67 + t68 * t66;
t1 = [-g(3) * qJ(1), 0, -t68, t60, 0, 0, 0, 0, 0, t59, t58, t59, t58, -t60, -g(1) * (t63 * t61 - t62 * t65 + cos(pkin(6)) * pkin(1)) - g(2) * (t62 * t61 + t63 * t65 + sin(pkin(6)) * pkin(1)) - g(3) * (t66 * pkin(3) + pkin(4) + qJ(1));];
U_reg = t1;
