% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR14V3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiR_rot_6_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:09
% EndTime: 2019-04-12 15:12:09
% DurationCPUTime: 0.15s
% Computational Cost: add. (90->38), mult. (280->82), div. (0->0), fcn. (405->10), ass. (0->43)
t158 = cos(qJ(4));
t153 = sin(qJ(4));
t160 = cos(qJ(1));
t164 = t160 * t153;
t155 = sin(qJ(1));
t159 = cos(qJ(2));
t168 = t155 * t159;
t143 = t158 * t168 - t164;
t157 = cos(qJ(5));
t152 = sin(qJ(5));
t154 = sin(qJ(2));
t172 = t154 * t152;
t134 = t143 * t157 + t155 * t172;
t163 = t160 * t158;
t142 = t153 * t168 + t163;
t151 = sin(qJ(6));
t156 = cos(qJ(6));
t179 = t134 * t151 - t142 * t156;
t178 = -t134 * t156 - t142 * t151;
t175 = t151 * t157;
t174 = t153 * t154;
t173 = t153 * t159;
t171 = t154 * t157;
t170 = t154 * t158;
t169 = t154 * t160;
t167 = t156 * t157;
t166 = t159 * t152;
t165 = t159 * t157;
t162 = t151 * t174;
t161 = t156 * t174;
t133 = -t143 * t152 + t155 * t171;
t141 = t157 * t170 - t166;
t140 = -t152 * t170 - t165;
t146 = t155 * t153 + t159 * t163;
t145 = -t155 * t158 + t159 * t164;
t144 = t158 * t165 + t172;
t139 = t141 * t160;
t138 = t141 * t155;
t137 = t146 * t157 + t152 * t169;
t136 = t146 * t152 - t157 * t169;
t132 = t137 * t156 + t145 * t151;
t131 = -t137 * t151 + t145 * t156;
t1 = [t178, -t139 * t156 - t160 * t162, 0, -t145 * t167 + t146 * t151, -t136 * t156, t131; t132, -t138 * t156 - t155 * t162, 0, -t142 * t167 + t143 * t151, t133 * t156, -t179; 0, t144 * t156 + t151 * t173, 0 (t151 * t158 - t153 * t167) * t154, t140 * t156, -t141 * t151 + t161; t179, t139 * t151 - t160 * t161, 0, t145 * t175 + t146 * t156, t136 * t151, -t132; t131, t138 * t151 - t155 * t161, 0, t142 * t175 + t143 * t156, -t133 * t151, t178; 0, -t144 * t151 + t156 * t173, 0 (t153 * t175 + t156 * t158) * t154, -t140 * t151, -t141 * t156 - t162; t133, t140 * t160, 0, -t145 * t152, t137, 0; t136, t140 * t155, 0, -t142 * t152, t134, 0; 0, t158 * t166 - t171, 0, -t153 * t172, t141, 0;];
JR_rot  = t1;
