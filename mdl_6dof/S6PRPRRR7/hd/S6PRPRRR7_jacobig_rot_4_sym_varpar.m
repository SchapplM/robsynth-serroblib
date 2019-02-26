% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR7_jacobig_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobig_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobig_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:30
% EndTime: 2019-02-26 19:57:30
% DurationCPUTime: 0.08s
% Computational Cost: add. (19->17), mult. (54->39), div. (0->0), fcn. (78->12), ass. (0->20)
t116 = sin(pkin(13));
t119 = sin(pkin(6));
t131 = t116 * t119;
t118 = sin(pkin(7));
t130 = t118 * t119;
t121 = cos(pkin(13));
t129 = t121 * t119;
t124 = cos(pkin(6));
t125 = sin(qJ(2));
t128 = t124 * t125;
t126 = cos(qJ(2));
t127 = t124 * t126;
t123 = cos(pkin(7));
t122 = cos(pkin(8));
t120 = cos(pkin(14));
t117 = sin(pkin(8));
t115 = sin(pkin(14));
t114 = -t116 * t127 - t121 * t125;
t113 = -t116 * t125 + t121 * t127;
t1 = [0, t131, 0 -(-(-t116 * t128 + t121 * t126) * t115 + (t114 * t123 + t116 * t130) * t120) * t117 + (-t114 * t118 + t123 * t131) * t122, 0, 0; 0, -t129, 0 -(-(t116 * t126 + t121 * t128) * t115 + (t113 * t123 - t118 * t129) * t120) * t117 + (-t113 * t118 - t123 * t129) * t122, 0, 0; 0, t124, 0 -(t124 * t118 * t120 + (t120 * t123 * t126 - t115 * t125) * t119) * t117 + (t124 * t123 - t126 * t130) * t122, 0, 0;];
Jg_rot  = t1;
