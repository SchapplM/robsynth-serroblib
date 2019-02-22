% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:24
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRPR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:24:12
% EndTime: 2019-02-22 09:24:12
% DurationCPUTime: 0.15s
% Computational Cost: add. (172->43), mult. (516->91), div. (0->0), fcn. (712->14), ass. (0->46)
t209 = cos(qJ(3));
t183 = sin(pkin(11));
t189 = cos(pkin(6));
t208 = t183 * t189;
t184 = sin(pkin(7));
t207 = t184 * t189;
t185 = sin(pkin(6));
t206 = t185 * t184;
t188 = cos(pkin(7));
t205 = t185 * t188;
t186 = cos(pkin(12));
t204 = t186 * t188;
t187 = cos(pkin(11));
t203 = t187 * t189;
t190 = sin(qJ(6));
t191 = sin(qJ(4));
t202 = t190 * t191;
t193 = cos(qJ(6));
t201 = t191 * t193;
t200 = t185 * t209;
t199 = t184 * t200;
t182 = sin(pkin(12));
t198 = -t183 * t182 + t186 * t203;
t197 = t187 * t182 + t186 * t208;
t196 = t198 * t188;
t195 = t197 * t188;
t194 = cos(qJ(4));
t192 = sin(qJ(3));
t178 = -t182 * t208 + t187 * t186;
t177 = t182 * t203 + t183 * t186;
t176 = -t186 * t206 + t189 * t188;
t173 = t183 * t205 + t197 * t184;
t172 = -t198 * t184 - t187 * t205;
t171 = t192 * t207 + (t209 * t182 + t192 * t204) * t185;
t170 = t185 * t182 * t192 - t200 * t204 - t209 * t207;
t169 = t171 * t194 + t176 * t191;
t168 = t171 * t191 - t176 * t194;
t167 = t178 * t209 + (t183 * t206 - t195) * t192;
t166 = t178 * t192 - t183 * t199 + t209 * t195;
t165 = t177 * t209 + (-t187 * t206 + t196) * t192;
t164 = t177 * t192 + t187 * t199 - t209 * t196;
t163 = t167 * t194 + t173 * t191;
t162 = t167 * t191 - t173 * t194;
t161 = t165 * t194 + t172 * t191;
t160 = t165 * t191 - t172 * t194;
t1 = [0, 0, -t166 * t202 + t167 * t193, t163 * t190, 0, t162 * t193 - t166 * t190; 0, 0, -t164 * t202 + t165 * t193, t161 * t190, 0, t160 * t193 - t164 * t190; 0, 0, -t170 * t202 + t171 * t193, t169 * t190, 0, t168 * t193 - t170 * t190; 0, 0, -t166 * t201 - t167 * t190, t163 * t193, 0, -t162 * t190 - t166 * t193; 0, 0, -t164 * t201 - t165 * t190, t161 * t193, 0, -t160 * t190 - t164 * t193; 0, 0, -t170 * t201 - t171 * t190, t169 * t193, 0, -t168 * t190 - t170 * t193; 0, 0, -t166 * t194, -t162, 0, 0; 0, 0, -t164 * t194, -t160, 0, 0; 0, 0, -t170 * t194, -t168, 0, 0;];
JR_rot  = t1;
