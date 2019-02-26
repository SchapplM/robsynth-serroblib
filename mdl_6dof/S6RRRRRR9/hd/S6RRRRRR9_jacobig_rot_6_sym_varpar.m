% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR9_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:07
% EndTime: 2019-02-26 22:52:08
% DurationCPUTime: 0.09s
% Computational Cost: add. (52->21), mult. (152->45), div. (0->0), fcn. (218->12), ass. (0->33)
t205 = sin(pkin(7));
t208 = cos(pkin(6));
t226 = t205 * t208;
t207 = cos(pkin(7));
t215 = cos(qJ(2));
t225 = t207 * t215;
t206 = sin(pkin(6));
t212 = sin(qJ(1));
t224 = t212 * t206;
t211 = sin(qJ(2));
t223 = t212 * t211;
t222 = t212 * t215;
t216 = cos(qJ(1));
t221 = t216 * t206;
t220 = t216 * t211;
t219 = t216 * t215;
t201 = t208 * t219 - t223;
t218 = -t201 * t207 + t205 * t221;
t203 = -t208 * t222 - t220;
t217 = t203 * t207 + t205 * t224;
t214 = cos(qJ(3));
t213 = cos(qJ(4));
t210 = sin(qJ(3));
t209 = sin(qJ(4));
t204 = -t208 * t223 + t219;
t202 = t208 * t220 + t222;
t200 = -t206 * t215 * t205 + t208 * t207;
t199 = -t203 * t205 + t207 * t224;
t198 = -t201 * t205 - t207 * t221;
t197 = (t210 * t226 + (t210 * t225 + t211 * t214) * t206) * t209 - t200 * t213;
t196 = (t204 * t214 + t210 * t217) * t209 - t199 * t213;
t195 = (t202 * t214 - t210 * t218) * t209 - t198 * t213;
t1 = [0, t224, t199, t204 * t210 - t214 * t217, t196, t196; 0, -t221, t198, t202 * t210 + t214 * t218, t195, t195; 1, t208, t200, -t214 * t226 + (t210 * t211 - t214 * t225) * t206, t197, t197;];
Jg_rot  = t1;
