% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR14_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:40
% EndTime: 2019-02-26 22:23:40
% DurationCPUTime: 0.03s
% Computational Cost: add. (14->9), mult. (39->18), div. (0->0), fcn. (63->8), ass. (0->18)
t117 = sin(pkin(6));
t121 = sin(qJ(1));
t130 = t121 * t117;
t120 = sin(qJ(2));
t129 = t121 * t120;
t123 = cos(qJ(2));
t128 = t121 * t123;
t124 = cos(qJ(1));
t127 = t124 * t117;
t126 = t124 * t120;
t125 = t124 * t123;
t122 = cos(qJ(3));
t119 = sin(qJ(3));
t118 = cos(pkin(6));
t116 = t117 * t120 * t122 + t118 * t119;
t115 = (-t118 * t129 + t125) * t122 + t119 * t130;
t114 = (t118 * t126 + t128) * t122 - t119 * t127;
t1 = [0, t130, t118 * t128 + t126, 0, t115, t115; 0, -t127, -t118 * t125 + t129, 0, t114, t114; 1, t118, -t117 * t123, 0, t116, t116;];
Jg_rot  = t1;
