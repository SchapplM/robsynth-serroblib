% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR9_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:24
% EndTime: 2019-02-26 20:53:24
% DurationCPUTime: 0.15s
% Computational Cost: add. (168->34), mult. (475->63), div. (0->0), fcn. (661->14), ass. (0->44)
t180 = sin(pkin(7));
t178 = sin(pkin(13));
t182 = cos(pkin(13));
t187 = sin(qJ(3));
t190 = cos(qJ(3));
t192 = t190 * t178 + t187 * t182;
t166 = t192 * t180;
t184 = cos(pkin(7));
t168 = t192 * t184;
t185 = cos(pkin(6));
t183 = cos(pkin(12));
t191 = cos(qJ(1));
t193 = t191 * t183;
t179 = sin(pkin(12));
t188 = sin(qJ(1));
t196 = t188 * t179;
t170 = -t185 * t193 + t196;
t194 = t191 * t179;
t195 = t188 * t183;
t171 = t185 * t194 + t195;
t174 = t187 * t178 - t190 * t182;
t181 = sin(pkin(6));
t197 = t181 * t191;
t157 = t166 * t197 + t170 * t168 + t171 * t174;
t162 = -t170 * t180 + t184 * t197;
t186 = sin(qJ(5));
t189 = cos(qJ(5));
t200 = t157 * t189 + t162 * t186;
t199 = t157 * t186 - t162 * t189;
t198 = t181 * t188;
t165 = t174 * t180;
t167 = t174 * t184;
t155 = t165 * t197 + t170 * t167 - t171 * t192;
t172 = -t185 * t195 - t194;
t173 = -t185 * t196 + t193;
t159 = t166 * t198 + t172 * t168 - t173 * t174;
t161 = t185 * t166 + (t168 * t183 - t174 * t179) * t181;
t169 = -t181 * t183 * t180 + t185 * t184;
t164 = -t172 * t180 + t184 * t198;
t160 = -t185 * t165 + (-t167 * t183 - t179 * t192) * t181;
t158 = -t165 * t198 - t172 * t167 - t173 * t192;
t154 = t159 * t189 + t164 * t186;
t153 = -t159 * t186 + t164 * t189;
t1 = [t200, 0, t158 * t189, 0, t153, 0; t154, 0, t155 * t189, 0, t199, 0; 0, 0, t160 * t189, 0, -t161 * t186 + t169 * t189, 0; -t199, 0, -t158 * t186, 0, -t154, 0; t153, 0, -t155 * t186, 0, t200, 0; 0, 0, -t160 * t186, 0, -t161 * t189 - t169 * t186, 0; t155, 0, t159, 0, 0, 0; -t158, 0, -t157, 0, 0, 0; 0, 0, t161, 0, 0, 0;];
JR_rot  = t1;
