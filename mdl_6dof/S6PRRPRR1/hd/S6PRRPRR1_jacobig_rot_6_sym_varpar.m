% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:06
% EndTime: 2019-02-26 20:04:06
% DurationCPUTime: 0.03s
% Computational Cost: add. (24->11), mult. (31->20), div. (0->0), fcn. (52->8), ass. (0->17)
t121 = sin(pkin(11));
t122 = sin(pkin(6));
t131 = t121 * t122;
t126 = cos(qJ(2));
t130 = t122 * t126;
t123 = cos(pkin(11));
t129 = t123 * t122;
t124 = cos(pkin(6));
t125 = sin(qJ(2));
t128 = t124 * t125;
t127 = t124 * t126;
t120 = qJ(3) + pkin(12) + qJ(5);
t119 = cos(t120);
t118 = sin(t120);
t117 = t121 * t127 + t123 * t125;
t116 = t121 * t125 - t123 * t127;
t1 = [0, t131, t117, 0, t117 (-t121 * t128 + t123 * t126) * t118 - t119 * t131; 0, -t129, t116, 0, t116 (t121 * t126 + t123 * t128) * t118 + t119 * t129; 0, t124, -t130, 0, -t130, t122 * t125 * t118 - t124 * t119;];
Jg_rot  = t1;
