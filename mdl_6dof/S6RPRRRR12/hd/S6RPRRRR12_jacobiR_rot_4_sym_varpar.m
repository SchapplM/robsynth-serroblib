% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:05
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR12_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobiR_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:05:34
% EndTime: 2019-02-22 11:05:34
% DurationCPUTime: 0.25s
% Computational Cost: add. (183->38), mult. (540->80), div. (0->0), fcn. (741->14), ass. (0->46)
t184 = cos(pkin(6));
t181 = cos(pkin(14));
t190 = cos(qJ(1));
t198 = t190 * t181;
t177 = sin(pkin(14));
t187 = sin(qJ(1));
t201 = t187 * t177;
t171 = -t184 * t198 + t201;
t199 = t190 * t177;
t200 = t187 * t181;
t172 = t184 * t199 + t200;
t186 = sin(qJ(3));
t189 = cos(qJ(3));
t179 = sin(pkin(7));
t180 = sin(pkin(6));
t205 = t180 * t190;
t197 = t179 * t205;
t183 = cos(pkin(7));
t202 = t183 * t186;
t162 = t171 * t202 - t172 * t189 + t186 * t197;
t185 = sin(qJ(4));
t188 = cos(qJ(4));
t161 = (t171 * t183 + t197) * t189 + t172 * t186;
t167 = -t171 * t179 + t183 * t205;
t178 = sin(pkin(8));
t182 = cos(pkin(8));
t195 = t161 * t182 + t167 * t178;
t213 = t162 * t188 + t195 * t185;
t212 = -t162 * t185 + t195 * t188;
t207 = t179 * t184;
t206 = t180 * t187;
t204 = t182 * t185;
t203 = t182 * t188;
t174 = -t184 * t201 + t198;
t173 = -t184 * t200 - t199;
t191 = t173 * t183 + t179 * t206;
t163 = -t174 * t186 + t191 * t189;
t169 = -t173 * t179 + t183 * t206;
t194 = t163 * t182 + t169 * t178;
t165 = t189 * t207 + (t181 * t183 * t189 - t177 * t186) * t180;
t193 = t165 * t182 + (-t180 * t181 * t179 + t184 * t183) * t178;
t166 = t186 * t207 + (t177 * t189 + t181 * t202) * t180;
t164 = t174 * t189 + t191 * t186;
t158 = t164 * t188 + t194 * t185;
t157 = -t164 * t185 + t194 * t188;
t1 = [t213, 0, t163 * t188 - t164 * t204, t157, 0, 0; t158, 0, -t161 * t188 + t162 * t204, -t212, 0, 0; 0, 0, t165 * t188 - t166 * t204, -t166 * t185 + t193 * t188, 0, 0; t212, 0, -t163 * t185 - t164 * t203, -t158, 0, 0; t157, 0, t161 * t185 + t162 * t203, t213, 0, 0; 0, 0, -t165 * t185 - t166 * t203, -t166 * t188 - t193 * t185, 0, 0; -t161 * t178 + t167 * t182, 0, t164 * t178, 0, 0, 0; -t163 * t178 + t169 * t182, 0, -t162 * t178, 0, 0, 0; 0, 0, t166 * t178, 0, 0, 0;];
JR_rot  = t1;
