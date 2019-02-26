% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR10_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:17
% EndTime: 2019-02-26 21:59:18
% DurationCPUTime: 0.03s
% Computational Cost: add. (26->10), mult. (39->18), div. (0->0), fcn. (63->8), ass. (0->19)
t127 = sin(pkin(6));
t130 = sin(qJ(1));
t138 = t130 * t127;
t129 = sin(qJ(2));
t137 = t130 * t129;
t131 = cos(qJ(2));
t136 = t130 * t131;
t132 = cos(qJ(1));
t135 = t132 * t127;
t134 = t132 * t129;
t133 = t132 * t131;
t128 = cos(pkin(6));
t126 = pkin(12) + qJ(4);
t125 = cos(t126);
t124 = sin(t126);
t123 = t127 * t129 * t124 - t128 * t125;
t122 = (-t128 * t137 + t133) * t124 - t125 * t138;
t121 = (t128 * t134 + t136) * t124 + t125 * t135;
t1 = [0, t138, 0, t128 * t136 + t134, t122, t122; 0, -t135, 0, -t128 * t133 + t137, t121, t121; 1, t128, 0, -t127 * t131, t123, t123;];
Jg_rot  = t1;
