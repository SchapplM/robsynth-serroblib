% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRP9_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:44:26
% EndTime: 2019-02-26 22:44:26
% DurationCPUTime: 0.03s
% Computational Cost: add. (14->9), mult. (39->18), div. (0->0), fcn. (63->8), ass. (0->18)
t121 = sin(pkin(6));
t125 = sin(qJ(1));
t134 = t125 * t121;
t124 = sin(qJ(2));
t133 = t125 * t124;
t127 = cos(qJ(2));
t132 = t125 * t127;
t128 = cos(qJ(1));
t131 = t128 * t121;
t130 = t128 * t124;
t129 = t128 * t127;
t126 = cos(qJ(3));
t123 = sin(qJ(3));
t122 = cos(pkin(6));
t120 = t121 * t124 * t123 - t122 * t126;
t119 = (-t122 * t133 + t129) * t123 - t126 * t134;
t118 = (t122 * t130 + t132) * t123 + t126 * t131;
t1 = [0, t134, t122 * t132 + t130, t119, t119, 0; 0, -t131, -t122 * t129 + t133, t118, t118, 0; 1, t122, -t121 * t127, t120, t120, 0;];
Jg_rot  = t1;
