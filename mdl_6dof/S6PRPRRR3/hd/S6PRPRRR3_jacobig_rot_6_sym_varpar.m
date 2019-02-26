% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR3_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:48
% EndTime: 2019-02-26 19:54:49
% DurationCPUTime: 0.08s
% Computational Cost: add. (24->11), mult. (31->20), div. (0->0), fcn. (52->8), ass. (0->17)
t118 = sin(pkin(11));
t119 = sin(pkin(6));
t128 = t118 * t119;
t123 = cos(qJ(2));
t127 = t119 * t123;
t120 = cos(pkin(11));
t126 = t120 * t119;
t121 = cos(pkin(6));
t122 = sin(qJ(2));
t125 = t121 * t122;
t124 = t121 * t123;
t117 = pkin(12) + qJ(4) + qJ(5);
t116 = cos(t117);
t115 = sin(t117);
t114 = t118 * t124 + t120 * t122;
t113 = t118 * t122 - t120 * t124;
t1 = [0, t128, 0, t114, t114 (-t118 * t125 + t120 * t123) * t115 - t116 * t128; 0, -t126, 0, t113, t113 (t118 * t123 + t120 * t125) * t115 + t116 * t126; 0, t121, 0, -t127, -t127, t119 * t122 * t115 - t121 * t116;];
Jg_rot  = t1;
