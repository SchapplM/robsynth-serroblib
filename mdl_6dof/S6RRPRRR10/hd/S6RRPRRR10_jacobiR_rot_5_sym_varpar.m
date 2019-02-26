% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
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
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR10_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:18
% EndTime: 2019-02-26 21:59:18
% DurationCPUTime: 0.10s
% Computational Cost: add. (125->31), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
t142 = cos(pkin(6));
t144 = sin(qJ(2));
t148 = cos(qJ(1));
t150 = t148 * t144;
t145 = sin(qJ(1));
t147 = cos(qJ(2));
t152 = t145 * t147;
t133 = t142 * t150 + t152;
t140 = pkin(12) + qJ(4);
t138 = sin(t140);
t139 = cos(t140);
t141 = sin(pkin(6));
t155 = t141 * t148;
t127 = -t133 * t139 + t138 * t155;
t149 = t148 * t147;
t153 = t145 * t144;
t132 = -t142 * t149 + t153;
t143 = sin(qJ(5));
t146 = cos(qJ(5));
t163 = t127 * t143 + t132 * t146;
t162 = t127 * t146 - t132 * t143;
t159 = t139 * t143;
t158 = t139 * t146;
t157 = t141 * t144;
t156 = t141 * t145;
t154 = t143 * t147;
t151 = t146 * t147;
t125 = -t133 * t138 - t139 * t155;
t135 = -t142 * t153 + t149;
t134 = t142 * t152 + t150;
t131 = t142 * t138 + t139 * t157;
t130 = -t138 * t157 + t142 * t139;
t129 = t135 * t139 + t138 * t156;
t128 = t135 * t138 - t139 * t156;
t124 = t129 * t146 + t134 * t143;
t123 = -t129 * t143 + t134 * t146;
t1 = [t162, -t134 * t158 + t135 * t143, 0, -t128 * t146, t123, 0; t124, -t132 * t158 + t133 * t143, 0, t125 * t146, t163, 0; 0 (t139 * t151 + t143 * t144) * t141, 0, t130 * t146, -t131 * t143 - t141 * t151, 0; -t163, t134 * t159 + t135 * t146, 0, t128 * t143, -t124, 0; t123, t132 * t159 + t133 * t146, 0, -t125 * t143, t162, 0; 0 (-t139 * t154 + t144 * t146) * t141, 0, -t130 * t143, -t131 * t146 + t141 * t154, 0; t125, -t134 * t138, 0, t129, 0, 0; t128, -t132 * t138, 0, -t127, 0, 0; 0, t141 * t147 * t138, 0, t131, 0, 0;];
JR_rot  = t1;
