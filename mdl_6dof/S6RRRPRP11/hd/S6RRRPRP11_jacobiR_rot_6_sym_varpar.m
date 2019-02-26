% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP11_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:59
% EndTime: 2019-02-26 22:14:59
% DurationCPUTime: 0.09s
% Computational Cost: add. (71->31), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->36)
t133 = cos(pkin(6));
t136 = sin(qJ(2));
t141 = cos(qJ(1));
t144 = t141 * t136;
t137 = sin(qJ(1));
t140 = cos(qJ(2));
t146 = t137 * t140;
t128 = t133 * t144 + t146;
t135 = sin(qJ(3));
t139 = cos(qJ(3));
t132 = sin(pkin(6));
t151 = t132 * t141;
t120 = t128 * t135 + t139 * t151;
t143 = t141 * t140;
t147 = t137 * t136;
t127 = -t133 * t143 + t147;
t134 = sin(qJ(5));
t138 = cos(qJ(5));
t157 = -t120 * t134 - t127 * t138;
t156 = t120 * t138 - t127 * t134;
t153 = t132 * t135;
t152 = t132 * t139;
t150 = t134 * t135;
t149 = t134 * t140;
t148 = t135 * t138;
t145 = t138 * t140;
t142 = -t128 * t139 + t135 * t151;
t130 = -t133 * t147 + t143;
t129 = t133 * t146 + t144;
t126 = t133 * t135 + t136 * t152;
t125 = -t133 * t139 + t136 * t153;
t124 = t130 * t139 + t137 * t153;
t123 = t130 * t135 - t137 * t152;
t119 = t123 * t134 + t129 * t138;
t118 = t123 * t138 - t129 * t134;
t1 = [t157, -t129 * t150 + t130 * t138, t124 * t134, 0, t118, 0; t119, -t127 * t150 + t128 * t138, -t142 * t134, 0, t156, 0; 0 (t135 * t149 + t136 * t138) * t132, t126 * t134, 0, t125 * t138 + t132 * t149, 0; -t156, -t129 * t148 - t130 * t134, t124 * t138, 0, -t119, 0; t118, -t127 * t148 - t128 * t134, -t142 * t138, 0, t157, 0; 0 (-t134 * t136 + t135 * t145) * t132, t126 * t138, 0, -t125 * t134 + t132 * t145, 0; t142, -t129 * t139, -t123, 0, 0, 0; t124, -t127 * t139, -t120, 0, 0, 0; 0, t140 * t152, -t125, 0, 0, 0;];
JR_rot  = t1;
