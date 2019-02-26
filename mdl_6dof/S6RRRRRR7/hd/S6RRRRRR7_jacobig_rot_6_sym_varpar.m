% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR7_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:42
% EndTime: 2019-02-26 22:50:42
% DurationCPUTime: 0.03s
% Computational Cost: add. (19->9), mult. (54->18), div. (0->0), fcn. (86->8), ass. (0->18)
t122 = sin(pkin(6));
t126 = sin(qJ(1));
t135 = t126 * t122;
t125 = sin(qJ(2));
t134 = t126 * t125;
t128 = cos(qJ(2));
t133 = t126 * t128;
t129 = cos(qJ(1));
t132 = t129 * t122;
t131 = t129 * t125;
t130 = t129 * t128;
t127 = cos(qJ(3));
t124 = sin(qJ(3));
t123 = cos(pkin(6));
t121 = t122 * t125 * t124 - t123 * t127;
t120 = (-t123 * t134 + t130) * t124 - t127 * t135;
t119 = (t123 * t131 + t133) * t124 + t127 * t132;
t1 = [0, t135, t123 * t133 + t131, t120, t120, t120; 0, -t132, -t123 * t130 + t134, t119, t119, t119; 1, t123, -t122 * t128, t121, t121, t121;];
Jg_rot  = t1;
