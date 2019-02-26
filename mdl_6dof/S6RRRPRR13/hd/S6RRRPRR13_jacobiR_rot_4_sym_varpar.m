% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR13_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobiR_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:04
% EndTime: 2019-02-26 22:23:04
% DurationCPUTime: 0.11s
% Computational Cost: add. (100->37), mult. (304->83), div. (0->0), fcn. (425->12), ass. (0->43)
t168 = sin(qJ(2));
t169 = sin(qJ(1));
t171 = cos(qJ(2));
t172 = cos(qJ(1));
t190 = cos(pkin(6));
t175 = t172 * t190;
t156 = t168 * t175 + t169 * t171;
t167 = sin(qJ(3));
t170 = cos(qJ(3));
t155 = t169 * t168 - t171 * t175;
t163 = sin(pkin(7));
t166 = cos(pkin(7));
t164 = sin(pkin(6));
t185 = t164 * t172;
t173 = t155 * t166 + t163 * t185;
t145 = -t156 * t170 + t173 * t167;
t162 = sin(pkin(13));
t188 = t162 * t163;
t165 = cos(pkin(13));
t187 = t163 * t165;
t186 = t164 * t169;
t184 = t166 * t167;
t183 = t166 * t170;
t182 = t167 * t168;
t181 = t167 * t171;
t180 = t168 * t170;
t179 = t170 * t171;
t178 = t163 * t164 * t168;
t177 = t163 * t186;
t176 = t169 * t190;
t174 = t190 * t163;
t144 = -t156 * t167 - t173 * t170;
t158 = -t168 * t176 + t172 * t171;
t157 = -t172 * t168 - t171 * t176;
t154 = (-t166 * t182 + t179) * t164;
t152 = -t157 * t163 + t166 * t186;
t151 = -t155 * t163 + t166 * t185;
t150 = t170 * t174 + (t166 * t179 - t182) * t164;
t149 = t157 * t170 - t158 * t184;
t148 = -t155 * t170 - t156 * t184;
t147 = t158 * t170 + (t157 * t166 + t177) * t167;
t146 = -t157 * t183 + t158 * t167 - t170 * t177;
t1 = [t145 * t165 + t151 * t162, t149 * t165 + t158 * t188, -t146 * t165, 0, 0, 0; t147 * t165 + t152 * t162, t148 * t165 + t156 * t188, t144 * t165, 0, 0, 0; 0, t154 * t165 + t162 * t178, t150 * t165, 0, 0, 0; -t145 * t162 + t151 * t165, -t149 * t162 + t158 * t187, t146 * t162, 0, 0, 0; -t147 * t162 + t152 * t165, -t148 * t162 + t156 * t187, -t144 * t162, 0, 0, 0; 0, -t154 * t162 + t165 * t178, -t150 * t162, 0, 0, 0; t144, t157 * t167 + t158 * t183, t147, 0, 0, 0; t146, -t155 * t167 + t156 * t183, -t145, 0, 0, 0; 0 (t166 * t180 + t181) * t164, t167 * t174 + (t166 * t181 + t180) * t164, 0, 0, 0;];
JR_rot  = t1;
