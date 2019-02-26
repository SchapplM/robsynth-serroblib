% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRP4_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:01
% EndTime: 2019-02-26 19:52:01
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->10), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->14)
t120 = sin(pkin(10));
t121 = sin(pkin(6));
t129 = t120 * t121;
t122 = cos(pkin(10));
t128 = t122 * t121;
t123 = cos(pkin(6));
t124 = sin(qJ(2));
t127 = t123 * t124;
t125 = cos(qJ(2));
t126 = t123 * t125;
t119 = pkin(11) + qJ(4);
t118 = cos(t119);
t117 = sin(t119);
t1 = [0, t129, 0, t120 * t126 + t122 * t124 (-t120 * t127 + t122 * t125) * t117 - t118 * t129, 0; 0, -t128, 0, t120 * t124 - t122 * t126 (t120 * t125 + t122 * t127) * t117 + t118 * t128, 0; 0, t123, 0, -t121 * t125, t121 * t124 * t117 - t123 * t118, 0;];
Jg_rot  = t1;
