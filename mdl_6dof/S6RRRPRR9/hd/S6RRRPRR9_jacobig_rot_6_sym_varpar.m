% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR9_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:38
% EndTime: 2019-02-26 22:20:38
% DurationCPUTime: 0.09s
% Computational Cost: add. (52->25), mult. (148->51), div. (0->0), fcn. (211->14), ass. (0->34)
t204 = sin(pkin(6));
t211 = sin(qJ(1));
t222 = t211 * t204;
t210 = sin(qJ(2));
t221 = t211 * t210;
t214 = cos(qJ(2));
t220 = t211 * t214;
t215 = cos(qJ(1));
t219 = t215 * t204;
t218 = t215 * t210;
t217 = t215 * t214;
t202 = sin(pkin(13));
t205 = cos(pkin(13));
t209 = sin(qJ(3));
t213 = cos(qJ(3));
t216 = t213 * t202 + t209 * t205;
t201 = -t209 * t202 + t213 * t205;
t212 = cos(qJ(5));
t208 = sin(qJ(5));
t207 = cos(pkin(6));
t206 = cos(pkin(7));
t203 = sin(pkin(7));
t199 = -t207 * t221 + t217;
t198 = -t207 * t220 - t218;
t197 = t207 * t218 + t220;
t196 = t207 * t217 - t221;
t195 = -t204 * t214 * t203 + t207 * t206;
t194 = t216 * t206;
t193 = t201 * t206;
t192 = t216 * t203;
t191 = t201 * t203;
t190 = -t198 * t203 + t206 * t222;
t189 = -t196 * t203 - t206 * t219;
t1 = [0, t222, t190, 0, -t191 * t222 - t198 * t193 + t199 * t216 (t192 * t222 + t198 * t194 + t199 * t201) * t208 - t190 * t212; 0, -t219, t189, 0, t191 * t219 - t196 * t193 + t197 * t216 (-t192 * t219 + t196 * t194 + t197 * t201) * t208 - t189 * t212; 1, t207, t195, 0, -t207 * t191 + (-t193 * t214 + t210 * t216) * t204 (t207 * t192 + (t194 * t214 + t201 * t210) * t204) * t208 - t195 * t212;];
Jg_rot  = t1;
