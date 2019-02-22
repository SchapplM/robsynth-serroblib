% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:24
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR13_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:24:35
% EndTime: 2019-02-22 12:24:35
% DurationCPUTime: 0.24s
% Computational Cost: add. (202->50), mult. (569->90), div. (0->0), fcn. (804->12), ass. (0->54)
t183 = cos(pkin(6));
t187 = sin(qJ(2));
t193 = cos(qJ(1));
t202 = t193 * t187;
t188 = sin(qJ(1));
t192 = cos(qJ(2));
t205 = t188 * t192;
t177 = t183 * t202 + t205;
t186 = sin(qJ(3));
t191 = cos(qJ(3));
t182 = sin(pkin(6));
t208 = t182 * t193;
t166 = t177 * t191 - t186 * t208;
t201 = t193 * t192;
t206 = t188 * t187;
t176 = -t183 * t201 + t206;
t185 = sin(qJ(4));
t190 = cos(qJ(4));
t153 = t166 * t185 - t176 * t190;
t154 = t166 * t190 + t176 * t185;
t184 = sin(qJ(6));
t189 = cos(qJ(6));
t200 = t153 * t189 - t154 * t184;
t199 = t153 * t184 + t154 * t189;
t211 = t182 * t186;
t210 = t182 * t191;
t209 = t182 * t192;
t207 = t185 * t191;
t204 = t190 * t191;
t203 = t191 * t192;
t179 = -t183 * t206 + t201;
t169 = t179 * t191 + t188 * t211;
t178 = t183 * t205 + t202;
t157 = t169 * t185 - t178 * t190;
t158 = t169 * t190 + t178 * t185;
t151 = t157 * t189 - t158 * t184;
t152 = t157 * t184 + t158 * t189;
t175 = t183 * t186 + t187 * t210;
t163 = t175 * t185 + t190 * t209;
t164 = t175 * t190 - t185 * t209;
t198 = t163 * t189 - t164 * t184;
t197 = t163 * t184 + t164 * t189;
t196 = -t184 * t190 + t185 * t189;
t195 = t184 * t185 + t189 * t190;
t194 = t177 * t186 + t191 * t208;
t174 = t183 * t191 - t187 * t211;
t171 = (t185 * t187 + t190 * t203) * t182;
t170 = (t185 * t203 - t187 * t190) * t182;
t168 = -t179 * t186 + t188 * t210;
t162 = -t178 * t204 + t179 * t185;
t161 = -t178 * t207 - t179 * t190;
t160 = -t176 * t204 + t177 * t185;
t159 = -t176 * t207 - t177 * t190;
t1 = [-t199, t161 * t184 + t162 * t189, t195 * t168, -t151, 0, t151; t152, t159 * t184 + t160 * t189, -t195 * t194, -t200, 0, t200; 0, t170 * t184 + t171 * t189, t195 * t174, -t198, 0, t198; -t200, t161 * t189 - t162 * t184, t196 * t168, t152, 0, -t152; t151, t159 * t189 - t160 * t184, -t196 * t194, t199, 0, -t199; 0, t170 * t189 - t171 * t184, t196 * t174, t197, 0, -t197; t194, t178 * t186, -t169, 0, 0, 0; t168, t176 * t186, -t166, 0, 0, 0; 0, -t186 * t209, -t175, 0, 0, 0;];
JR_rot  = t1;
