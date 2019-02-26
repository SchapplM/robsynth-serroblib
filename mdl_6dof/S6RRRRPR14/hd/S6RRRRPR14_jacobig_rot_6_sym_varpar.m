% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPR14_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:12
% EndTime: 2019-02-26 22:38:12
% DurationCPUTime: 0.08s
% Computational Cost: add. (34->21), mult. (100->45), div. (0->0), fcn. (145->12), ass. (0->30)
t189 = sin(pkin(7));
t192 = cos(pkin(6));
t210 = t189 * t192;
t191 = cos(pkin(7));
t199 = cos(qJ(2));
t209 = t191 * t199;
t190 = sin(pkin(6));
t196 = sin(qJ(1));
t208 = t196 * t190;
t195 = sin(qJ(2));
t207 = t196 * t195;
t206 = t196 * t199;
t200 = cos(qJ(1));
t205 = t200 * t190;
t204 = t200 * t195;
t203 = t200 * t199;
t185 = t192 * t203 - t207;
t202 = -t185 * t191 + t189 * t205;
t187 = -t192 * t206 - t204;
t201 = t187 * t191 + t189 * t208;
t198 = cos(qJ(3));
t197 = cos(qJ(4));
t194 = sin(qJ(3));
t193 = sin(qJ(4));
t188 = -t192 * t207 + t203;
t186 = t192 * t204 + t206;
t184 = -t190 * t199 * t189 + t192 * t191;
t183 = -t187 * t189 + t191 * t208;
t182 = -t185 * t189 - t191 * t205;
t1 = [0, t208, t183, t188 * t194 - t201 * t198, 0 (t188 * t198 + t201 * t194) * t193 - t183 * t197; 0, -t205, t182, t186 * t194 + t202 * t198, 0 (t186 * t198 - t202 * t194) * t193 - t182 * t197; 1, t192, t184, -t198 * t210 + (t194 * t195 - t198 * t209) * t190, 0 (t194 * t210 + (t194 * t209 + t195 * t198) * t190) * t193 - t184 * t197;];
Jg_rot  = t1;
