% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRR1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:37
% EndTime: 2019-02-26 20:18:37
% DurationCPUTime: 0.03s
% Computational Cost: add. (27->12), mult. (38->20), div. (0->0), fcn. (64->8), ass. (0->17)
t124 = sin(pkin(12));
t125 = sin(pkin(6));
t134 = t124 * t125;
t129 = cos(qJ(2));
t133 = t125 * t129;
t126 = cos(pkin(12));
t132 = t126 * t125;
t127 = cos(pkin(6));
t128 = sin(qJ(2));
t131 = t127 * t128;
t130 = t127 * t129;
t123 = qJ(3) + qJ(4) + qJ(5);
t122 = cos(t123);
t121 = sin(t123);
t120 = t124 * t130 + t126 * t128;
t119 = t124 * t128 - t126 * t130;
t1 = [0, t134, t120, t120, t120 (-t124 * t131 + t126 * t129) * t121 - t122 * t134; 0, -t132, t119, t119, t119 (t124 * t129 + t126 * t131) * t121 + t122 * t132; 0, t127, -t133, -t133, -t133, t125 * t128 * t121 - t127 * t122;];
Jg_rot  = t1;
