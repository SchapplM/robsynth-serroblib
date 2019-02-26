% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR13_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:01:09
% EndTime: 2019-02-26 22:01:09
% DurationCPUTime: 0.03s
% Computational Cost: add. (13->8), mult. (39->18), div. (0->0), fcn. (63->8), ass. (0->18)
t118 = sin(pkin(6));
t122 = sin(qJ(1));
t131 = t122 * t118;
t121 = sin(qJ(2));
t130 = t122 * t121;
t124 = cos(qJ(2));
t129 = t122 * t124;
t125 = cos(qJ(1));
t128 = t125 * t118;
t127 = t125 * t121;
t126 = t125 * t124;
t123 = cos(qJ(4));
t120 = sin(qJ(4));
t119 = cos(pkin(6));
t117 = t118 * t124 * t123 + t119 * t120;
t116 = -t120 * t128 - (-t119 * t126 + t130) * t123;
t115 = t120 * t131 - (t119 * t129 + t127) * t123;
t1 = [0, t131, 0, -t119 * t130 + t126, t115, t115; 0, -t128, 0, t119 * t127 + t129, t116, t116; 1, t119, 0, t118 * t121, t117, t117;];
Jg_rot  = t1;
