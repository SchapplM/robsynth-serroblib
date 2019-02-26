% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR10_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:50
% EndTime: 2019-02-26 21:19:50
% DurationCPUTime: 0.15s
% Computational Cost: add. (179->32), mult. (411->59), div. (0->0), fcn. (576->12), ass. (0->43)
t186 = cos(pkin(6));
t181 = sin(pkin(13));
t190 = cos(qJ(1));
t194 = t190 * t181;
t184 = cos(pkin(13));
t188 = sin(qJ(1));
t195 = t188 * t184;
t173 = t186 * t194 + t195;
t187 = sin(qJ(3));
t189 = cos(qJ(3));
t193 = t190 * t184;
t196 = t188 * t181;
t172 = -t186 * t193 + t196;
t182 = sin(pkin(7));
t185 = cos(pkin(7));
t183 = sin(pkin(6));
t198 = t183 * t190;
t191 = t172 * t185 + t182 * t198;
t162 = -t173 * t189 + t187 * t191;
t167 = -t172 * t182 + t185 * t198;
t180 = qJ(4) + qJ(5);
t178 = sin(t180);
t179 = cos(t180);
t154 = t162 * t178 - t167 * t179;
t155 = t162 * t179 + t167 * t178;
t200 = t182 * t186;
t199 = t183 * t188;
t197 = t185 * t189;
t192 = t182 * t199;
t160 = -t173 * t187 - t189 * t191;
t175 = -t186 * t196 + t193;
t174 = -t186 * t195 - t194;
t171 = -t182 * t183 * t184 + t185 * t186;
t169 = -t174 * t182 + t185 * t199;
t166 = t187 * t200 + (t184 * t185 * t187 + t181 * t189) * t183;
t165 = t189 * t200 + (-t181 * t187 + t184 * t197) * t183;
t164 = t175 * t189 + (t174 * t185 + t192) * t187;
t163 = -t174 * t197 + t175 * t187 - t189 * t192;
t159 = -t166 * t179 - t171 * t178;
t158 = -t166 * t178 + t171 * t179;
t157 = t164 * t179 + t169 * t178;
t156 = -t164 * t178 + t169 * t179;
t1 = [t155, 0, -t163 * t179, t156, t156, 0; t157, 0, t160 * t179, t154, t154, 0; 0, 0, t165 * t179, t158, t158, 0; -t154, 0, t163 * t178, -t157, -t157, 0; t156, 0, -t160 * t178, t155, t155, 0; 0, 0, -t165 * t178, t159, t159, 0; t160, 0, t164, 0, 0, 0; t163, 0, -t162, 0, 0, 0; 0, 0, t166, 0, 0, 0;];
JR_rot  = t1;
