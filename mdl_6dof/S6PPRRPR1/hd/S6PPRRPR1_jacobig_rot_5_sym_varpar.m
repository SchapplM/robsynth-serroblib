% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRPR1_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:29
% EndTime: 2019-02-26 19:40:29
% DurationCPUTime: 0.04s
% Computational Cost: add. (15->13), mult. (46->32), div. (0->0), fcn. (67->10), ass. (0->17)
t123 = sin(pkin(11));
t129 = cos(pkin(6));
t135 = t123 * t129;
t124 = sin(pkin(7));
t125 = sin(pkin(6));
t134 = t124 * t125;
t128 = cos(pkin(7));
t133 = t125 * t128;
t127 = cos(pkin(11));
t132 = t127 * t129;
t131 = cos(qJ(3));
t130 = sin(qJ(3));
t126 = cos(pkin(12));
t122 = sin(pkin(12));
t121 = -t127 * t122 - t126 * t135;
t120 = -t123 * t122 + t126 * t132;
t1 = [0, 0, -t121 * t124 + t123 * t133 (-t122 * t135 + t127 * t126) * t130 + (-t121 * t128 - t123 * t134) * t131, 0, 0; 0, 0, -t120 * t124 - t127 * t133 (t122 * t132 + t123 * t126) * t130 + (-t120 * t128 + t127 * t134) * t131, 0, 0; 0, 0, -t126 * t134 + t129 * t128, -t129 * t124 * t131 + (-t126 * t128 * t131 + t122 * t130) * t125, 0, 0;];
Jg_rot  = t1;
