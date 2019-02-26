% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRP1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:13
% EndTime: 2019-02-26 20:15:13
% DurationCPUTime: 0.03s
% Computational Cost: add. (18->11), mult. (31->20), div. (0->0), fcn. (52->8), ass. (0->17)
t126 = sin(pkin(11));
t127 = sin(pkin(6));
t136 = t126 * t127;
t131 = cos(qJ(2));
t135 = t127 * t131;
t128 = cos(pkin(11));
t134 = t128 * t127;
t129 = cos(pkin(6));
t130 = sin(qJ(2));
t133 = t129 * t130;
t132 = t129 * t131;
t125 = qJ(3) + qJ(4);
t124 = cos(t125);
t123 = sin(t125);
t122 = t126 * t132 + t128 * t130;
t121 = t126 * t130 - t128 * t132;
t1 = [0, t136, t122, t122 (-t126 * t133 + t128 * t131) * t123 - t124 * t136, 0; 0, -t134, t121, t121 (t126 * t131 + t128 * t133) * t123 + t124 * t134, 0; 0, t129, -t135, -t135, t127 * t130 * t123 - t129 * t124, 0;];
Jg_rot  = t1;
