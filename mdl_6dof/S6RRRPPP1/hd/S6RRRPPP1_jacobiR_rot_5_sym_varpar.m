% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPP1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:48
% EndTime: 2019-02-26 22:02:48
% DurationCPUTime: 0.11s
% Computational Cost: add. (53->27), mult. (183->71), div. (0->0), fcn. (258->10), ass. (0->35)
t145 = sin(pkin(10));
t148 = cos(pkin(6));
t174 = t145 * t148;
t146 = sin(pkin(6));
t150 = sin(qJ(2));
t173 = t146 * t150;
t147 = cos(pkin(10));
t172 = t147 * t148;
t152 = cos(qJ(3));
t171 = t148 * t152;
t149 = sin(qJ(3));
t153 = cos(qJ(2));
t170 = t149 * t153;
t169 = t150 * t148;
t151 = sin(qJ(1));
t168 = t150 * t151;
t167 = t150 * t152;
t154 = cos(qJ(1));
t166 = t150 * t154;
t165 = t151 * t149;
t164 = t152 * t153;
t163 = t154 * t149;
t162 = t154 * t152;
t141 = t153 * t165 + t162;
t161 = -t141 * t148 + t146 * t168;
t143 = t151 * t152 - t153 * t163;
t160 = -t143 * t148 - t146 * t166;
t159 = t148 * t170 - t173;
t158 = -t146 * t153 - t149 * t169;
t157 = t148 * t153 - t149 * t173;
t156 = -t145 * t167 + t158 * t147;
t155 = t158 * t145 + t147 * t167;
t144 = t153 * t162 + t165;
t142 = -t151 * t164 + t163;
t1 = [-t141 * t146 - t148 * t168, t157 * t154, t144 * t146, 0, 0, 0; -t143 * t146 + t148 * t166, t157 * t151, -t142 * t146, 0, 0, 0; 0, t146 * t170 + t169, t146 * t167, 0, 0, 0; -t142 * t147 + t161 * t145, t155 * t154, -t143 * t147 + t144 * t174, 0, 0, 0; -t144 * t147 + t160 * t145, t155 * t151, t141 * t147 - t142 * t174, 0, 0, 0; 0, t159 * t145 - t147 * t164 (t145 * t171 + t147 * t149) * t150, 0, 0, 0; t142 * t145 + t161 * t147, t156 * t154, t143 * t145 + t144 * t172, 0, 0, 0; t144 * t145 + t160 * t147, t156 * t151, -t141 * t145 - t142 * t172, 0, 0, 0; 0, t145 * t164 + t159 * t147 (-t145 * t149 + t147 * t171) * t150, 0, 0, 0;];
JR_rot  = t1;
