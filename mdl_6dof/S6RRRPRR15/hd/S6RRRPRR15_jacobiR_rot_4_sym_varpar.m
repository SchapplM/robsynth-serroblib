% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR15_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiR_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:20
% EndTime: 2019-02-26 22:24:20
% DurationCPUTime: 0.08s
% Computational Cost: add. (58->24), mult. (178->57), div. (0->0), fcn. (255->10), ass. (0->32)
t143 = sin(qJ(2));
t144 = sin(qJ(1));
t146 = cos(qJ(2));
t147 = cos(qJ(1));
t163 = cos(pkin(6));
t152 = t147 * t163;
t134 = t143 * t152 + t144 * t146;
t142 = sin(qJ(3));
t145 = cos(qJ(3));
t133 = t144 * t143 - t146 * t152;
t139 = sin(pkin(7));
t141 = cos(pkin(7));
t140 = sin(pkin(6));
t160 = t140 * t147;
t150 = t133 * t141 + t139 * t160;
t164 = t134 * t142 + t150 * t145;
t161 = t140 * t144;
t159 = t141 * t142;
t158 = t141 * t145;
t157 = t142 * t143;
t156 = t142 * t146;
t155 = t143 * t145;
t154 = t145 * t146;
t153 = t144 * t163;
t151 = t163 * t139;
t135 = -t147 * t143 - t146 * t153;
t149 = t135 * t141 + t139 * t161;
t148 = t134 * t145 - t150 * t142;
t136 = -t143 * t153 + t147 * t146;
t132 = t136 * t145 + t149 * t142;
t131 = t136 * t142 - t149 * t145;
t1 = [-t133 * t139 + t141 * t160, t136 * t139, 0, 0, 0, 0; -t135 * t139 + t141 * t161, t134 * t139, 0, 0, 0, 0; 0, t140 * t143 * t139, 0, 0, 0, 0; t148, -t135 * t145 + t136 * t159, t131, 0, 0, 0; -t132, t133 * t145 + t134 * t159, t164, 0, 0, 0; 0 (t141 * t157 - t154) * t140, -t145 * t151 + (-t141 * t154 + t157) * t140, 0, 0, 0; -t164, t135 * t142 + t136 * t158, t132, 0, 0, 0; t131, -t133 * t142 + t134 * t158, t148, 0, 0, 0; 0 (t141 * t155 + t156) * t140, t142 * t151 + (t141 * t156 + t155) * t140, 0, 0, 0;];
JR_rot  = t1;
