% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRP2_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:48
% EndTime: 2019-02-26 20:01:48
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->10), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->14)
t123 = sin(pkin(10));
t124 = sin(pkin(6));
t132 = t123 * t124;
t125 = cos(pkin(10));
t131 = t125 * t124;
t126 = cos(pkin(6));
t127 = sin(qJ(2));
t130 = t126 * t127;
t128 = cos(qJ(2));
t129 = t126 * t128;
t122 = qJ(3) + pkin(11);
t121 = cos(t122);
t120 = sin(t122);
t1 = [0, t132, t123 * t129 + t125 * t127, 0 (-t123 * t130 + t125 * t128) * t120 - t121 * t132, 0; 0, -t131, t123 * t127 - t125 * t129, 0 (t123 * t128 + t125 * t130) * t120 + t121 * t131, 0; 0, t126, -t124 * t128, 0, t124 * t127 * t120 - t126 * t121, 0;];
Jg_rot  = t1;
