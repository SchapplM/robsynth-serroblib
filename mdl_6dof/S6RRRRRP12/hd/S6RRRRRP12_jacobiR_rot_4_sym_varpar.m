% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:32
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP12_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobiR_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:32:47
% EndTime: 2019-02-22 12:32:47
% DurationCPUTime: 0.17s
% Computational Cost: add. (136->42), mult. (408->90), div. (0->0), fcn. (571->12), ass. (0->50)
t184 = sin(qJ(2));
t185 = sin(qJ(1));
t188 = cos(qJ(2));
t189 = cos(qJ(1));
t208 = cos(pkin(6));
t192 = t189 * t208;
t173 = t184 * t192 + t185 * t188;
t183 = sin(qJ(3));
t187 = cos(qJ(3));
t172 = t185 * t184 - t188 * t192;
t179 = sin(pkin(7));
t181 = cos(pkin(7));
t180 = sin(pkin(6));
t202 = t180 * t189;
t190 = t172 * t181 + t179 * t202;
t159 = -t173 * t187 + t190 * t183;
t166 = -t172 * t179 + t181 * t202;
t182 = sin(qJ(4));
t186 = cos(qJ(4));
t212 = t159 * t182 - t166 * t186;
t211 = t159 * t186 + t166 * t182;
t206 = t179 * t180;
t205 = t179 * t182;
t204 = t179 * t186;
t203 = t180 * t185;
t201 = t181 * t183;
t200 = t181 * t187;
t199 = t183 * t184;
t198 = t183 * t188;
t197 = t184 * t187;
t196 = t187 * t188;
t195 = t184 * t206;
t194 = t179 * t203;
t193 = t185 * t208;
t191 = t208 * t179;
t157 = -t173 * t183 - t190 * t187;
t175 = -t184 * t193 + t189 * t188;
t174 = -t189 * t184 - t188 * t193;
t171 = t208 * t181 - t188 * t206;
t170 = (-t181 * t199 + t196) * t180;
t168 = -t174 * t179 + t181 * t203;
t165 = t183 * t191 + (t181 * t198 + t197) * t180;
t164 = t187 * t191 + (t181 * t196 - t199) * t180;
t163 = t174 * t187 - t175 * t201;
t162 = -t172 * t187 - t173 * t201;
t161 = t175 * t187 + (t174 * t181 + t194) * t183;
t160 = -t174 * t200 + t175 * t183 - t187 * t194;
t156 = t161 * t186 + t168 * t182;
t155 = -t161 * t182 + t168 * t186;
t1 = [t211, t163 * t186 + t175 * t205, -t160 * t186, t155, 0, 0; t156, t162 * t186 + t173 * t205, t157 * t186, t212, 0, 0; 0, t170 * t186 + t182 * t195, t164 * t186, -t165 * t182 + t171 * t186, 0, 0; -t212, -t163 * t182 + t175 * t204, t160 * t182, -t156, 0, 0; t155, -t162 * t182 + t173 * t204, -t157 * t182, t211, 0, 0; 0, -t170 * t182 + t186 * t195, -t164 * t182, -t165 * t186 - t171 * t182, 0, 0; t157, t174 * t183 + t175 * t200, t161, 0, 0, 0; t160, -t172 * t183 + t173 * t200, -t159, 0, 0, 0; 0 (t181 * t197 + t198) * t180, t165, 0, 0, 0;];
JR_rot  = t1;
