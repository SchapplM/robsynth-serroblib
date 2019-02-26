% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPPRRR1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobig_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:44
% EndTime: 2019-02-26 19:38:44
% DurationCPUTime: 0.08s
% Computational Cost: add. (101->32), mult. (294->70), div. (0->0), fcn. (404->16), ass. (0->44)
t186 = sin(pkin(12));
t195 = cos(pkin(6));
t210 = t186 * t195;
t188 = sin(pkin(7));
t189 = sin(pkin(6));
t209 = t188 * t189;
t208 = t188 * t195;
t194 = cos(pkin(7));
t207 = t189 * t194;
t191 = cos(pkin(13));
t206 = t191 * t194;
t192 = cos(pkin(12));
t205 = t192 * t195;
t185 = sin(pkin(13));
t181 = t185 * t205 + t186 * t191;
t184 = sin(pkin(14));
t190 = cos(pkin(14));
t180 = -t186 * t185 + t191 * t205;
t201 = t180 * t194 - t192 * t209;
t170 = -t181 * t184 + t201 * t190;
t177 = -t180 * t188 - t192 * t207;
t187 = sin(pkin(8));
t193 = cos(pkin(8));
t204 = t170 * t193 + t177 * t187;
t183 = -t185 * t210 + t192 * t191;
t182 = -t192 * t185 - t191 * t210;
t200 = t182 * t194 + t186 * t209;
t172 = -t183 * t184 + t200 * t190;
t178 = -t182 * t188 + t186 * t207;
t203 = t172 * t193 + t178 * t187;
t175 = t190 * t208 + (-t184 * t185 + t190 * t206) * t189;
t179 = -t191 * t209 + t195 * t194;
t202 = t175 * t193 + t179 * t187;
t199 = cos(qJ(4));
t198 = cos(qJ(5));
t197 = sin(qJ(4));
t196 = sin(qJ(5));
t176 = t189 * t185 * t190 + (t189 * t206 + t208) * t184;
t174 = -t175 * t187 + t179 * t193;
t173 = t183 * t190 + t200 * t184;
t171 = t181 * t190 + t201 * t184;
t169 = -t172 * t187 + t178 * t193;
t168 = -t170 * t187 + t177 * t193;
t1 = [0, 0, 0, t169, t173 * t197 - t203 * t199 (t173 * t199 + t203 * t197) * t196 - t169 * t198; 0, 0, 0, t168, t171 * t197 - t204 * t199 (t171 * t199 + t204 * t197) * t196 - t168 * t198; 0, 0, 0, t174, t176 * t197 - t202 * t199 (t176 * t199 + t202 * t197) * t196 - t174 * t198;];
Jg_rot  = t1;
