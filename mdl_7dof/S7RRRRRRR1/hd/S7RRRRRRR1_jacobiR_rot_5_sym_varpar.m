% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JR_rot [9x7]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S7RRRRRRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiR_rot_5_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiR_rot_5_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:53
% EndTime: 2018-11-26 21:20:53
% DurationCPUTime: 0.17s
% Computational Cost: add. (90->37), mult. (280->82), div. (0->0), fcn. (405->10), ass. (0->43)
t166 = cos(qJ(3));
t161 = sin(qJ(3));
t168 = cos(qJ(1));
t172 = t168 * t161;
t163 = sin(qJ(1));
t167 = cos(qJ(2));
t176 = t163 * t167;
t151 = t166 * t176 + t172;
t165 = cos(qJ(4));
t160 = sin(qJ(4));
t162 = sin(qJ(2));
t180 = t162 * t160;
t142 = t151 * t165 + t163 * t180;
t171 = t168 * t166;
t150 = t161 * t176 - t171;
t159 = sin(qJ(5));
t164 = cos(qJ(5));
t187 = t142 * t159 + t150 * t164;
t186 = -t142 * t164 + t150 * t159;
t183 = t159 * t165;
t182 = t161 * t162;
t181 = t161 * t167;
t179 = t162 * t165;
t178 = t162 * t166;
t177 = t162 * t168;
t175 = t164 * t165;
t174 = t167 * t160;
t173 = t167 * t165;
t170 = t159 * t182;
t169 = t164 * t182;
t141 = -t151 * t160 + t163 * t179;
t149 = t165 * t178 - t174;
t148 = -t160 * t178 - t173;
t154 = -t163 * t161 + t167 * t171;
t153 = -t163 * t166 - t167 * t172;
t152 = t166 * t173 + t180;
t147 = t149 * t168;
t146 = t149 * t163;
t145 = t154 * t165 + t160 * t177;
t144 = t154 * t160 - t165 * t177;
t140 = t145 * t164 + t153 * t159;
t139 = -t145 * t159 + t153 * t164;
t1 = [t186, -t147 * t164 + t168 * t170, t153 * t175 - t154 * t159, -t144 * t164, t139, 0, 0; t140, -t146 * t164 + t163 * t170, -t150 * t175 - t151 * t159, t141 * t164, -t187, 0, 0; 0, t152 * t164 - t159 * t181 (-t159 * t166 - t161 * t175) * t162, t148 * t164, -t149 * t159 - t169, 0, 0; t187, t147 * t159 + t168 * t169, -t153 * t183 - t154 * t164, t144 * t159, -t140, 0, 0; t139, t146 * t159 + t163 * t169, t150 * t183 - t151 * t164, -t141 * t159, t186, 0, 0; 0, -t152 * t159 - t164 * t181 (t161 * t183 - t164 * t166) * t162, -t148 * t159, -t149 * t164 + t170, 0, 0; t141, t148 * t168, t153 * t160, t145, 0, 0, 0; t144, t148 * t163, -t150 * t160, t142, 0, 0, 0; 0, t166 * t174 - t179, -t161 * t180, t149, 0, 0, 0;];
JR_rot  = t1;
