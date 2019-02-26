% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR8_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:07
% EndTime: 2019-02-26 20:08:07
% DurationCPUTime: 0.04s
% Computational Cost: add. (16->14), mult. (48->33), div. (0->0), fcn. (72->10), ass. (0->18)
t120 = sin(pkin(12));
t122 = sin(pkin(6));
t134 = t120 * t122;
t121 = sin(pkin(7));
t133 = t121 * t122;
t123 = cos(pkin(12));
t132 = t123 * t122;
t125 = cos(pkin(6));
t127 = sin(qJ(2));
t131 = t125 * t127;
t129 = cos(qJ(2));
t130 = t125 * t129;
t128 = cos(qJ(3));
t126 = sin(qJ(3));
t124 = cos(pkin(7));
t119 = -t120 * t130 - t123 * t127;
t118 = -t120 * t127 + t123 * t130;
t1 = [0, t134, -t119 * t121 + t124 * t134, 0 (-t120 * t131 + t123 * t129) * t128 + (t119 * t124 + t120 * t133) * t126, 0; 0, -t132, -t118 * t121 - t124 * t132, 0 (t120 * t129 + t123 * t131) * t128 + (t118 * t124 - t121 * t132) * t126, 0; 0, t125, t125 * t124 - t129 * t133, 0, t125 * t121 * t126 + (t124 * t126 * t129 + t127 * t128) * t122, 0;];
Jg_rot  = t1;
