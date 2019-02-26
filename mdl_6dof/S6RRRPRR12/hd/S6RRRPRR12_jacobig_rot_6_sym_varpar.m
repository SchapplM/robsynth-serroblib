% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR12_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:22:35
% EndTime: 2019-02-26 22:22:35
% DurationCPUTime: 0.03s
% Computational Cost: add. (14->9), mult. (39->18), div. (0->0), fcn. (63->8), ass. (0->18)
t119 = sin(pkin(6));
t123 = sin(qJ(1));
t132 = t123 * t119;
t122 = sin(qJ(2));
t131 = t123 * t122;
t125 = cos(qJ(2));
t130 = t123 * t125;
t126 = cos(qJ(1));
t129 = t126 * t119;
t128 = t126 * t122;
t127 = t126 * t125;
t124 = cos(qJ(3));
t121 = sin(qJ(3));
t120 = cos(pkin(6));
t118 = t119 * t122 * t121 - t120 * t124;
t117 = (-t120 * t131 + t127) * t121 - t124 * t132;
t116 = (t120 * t128 + t130) * t121 + t124 * t129;
t1 = [0, t132, t120 * t130 + t128, 0, t117, t117; 0, -t129, -t120 * t127 + t131, 0, t116, t116; 1, t120, -t119 * t125, 0, t118, t118;];
Jg_rot  = t1;
