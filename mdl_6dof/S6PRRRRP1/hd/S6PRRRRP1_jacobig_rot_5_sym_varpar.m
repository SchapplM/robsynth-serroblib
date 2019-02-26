% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRP1_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobig_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:13
% EndTime: 2019-02-26 20:15:13
% DurationCPUTime: 0.09s
% Computational Cost: add. (18->11), mult. (31->20), div. (0->0), fcn. (52->8), ass. (0->17)
t117 = sin(pkin(11));
t118 = sin(pkin(6));
t127 = t117 * t118;
t122 = cos(qJ(2));
t126 = t118 * t122;
t119 = cos(pkin(11));
t125 = t119 * t118;
t120 = cos(pkin(6));
t121 = sin(qJ(2));
t124 = t120 * t121;
t123 = t120 * t122;
t116 = qJ(3) + qJ(4);
t115 = cos(t116);
t114 = sin(t116);
t113 = t117 * t123 + t119 * t121;
t112 = t117 * t121 - t119 * t123;
t1 = [0, t127, t113, t113 (-t117 * t124 + t119 * t122) * t114 - t115 * t127, 0; 0, -t125, t112, t112 (t117 * t122 + t119 * t124) * t114 + t115 * t125, 0; 0, t120, -t126, -t126, t118 * t121 * t114 - t120 * t115, 0;];
Jg_rot  = t1;
