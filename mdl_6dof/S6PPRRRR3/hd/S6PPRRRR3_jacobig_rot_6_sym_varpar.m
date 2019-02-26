% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRRR3_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobig_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:56
% EndTime: 2019-02-26 19:43:57
% DurationCPUTime: 0.10s
% Computational Cost: add. (106->32), mult. (309->70), div. (0->0), fcn. (427->16), ass. (0->44)
t206 = sin(pkin(13));
t214 = cos(pkin(6));
t231 = t206 * t214;
t208 = sin(pkin(7));
t209 = sin(pkin(6));
t230 = t208 * t209;
t229 = t208 * t214;
t213 = cos(pkin(7));
t228 = t209 * t213;
t210 = cos(pkin(14));
t227 = t210 * t213;
t211 = cos(pkin(13));
t226 = t211 * t214;
t205 = sin(pkin(14));
t202 = t205 * t226 + t206 * t210;
t217 = sin(qJ(3));
t220 = cos(qJ(3));
t201 = -t206 * t205 + t210 * t226;
t222 = t201 * t213 - t211 * t230;
t191 = -t202 * t217 + t222 * t220;
t198 = -t201 * t208 - t211 * t228;
t207 = sin(pkin(8));
t212 = cos(pkin(8));
t225 = t191 * t212 + t198 * t207;
t204 = -t205 * t231 + t211 * t210;
t203 = -t211 * t205 - t210 * t231;
t221 = t203 * t213 + t206 * t230;
t193 = -t204 * t217 + t221 * t220;
t199 = -t203 * t208 + t206 * t228;
t224 = t193 * t212 + t199 * t207;
t196 = t220 * t229 + (-t205 * t217 + t220 * t227) * t209;
t200 = -t210 * t230 + t214 * t213;
t223 = t196 * t212 + t200 * t207;
t219 = cos(qJ(4));
t218 = cos(qJ(5));
t216 = sin(qJ(4));
t215 = sin(qJ(5));
t197 = t217 * t229 + (t205 * t220 + t217 * t227) * t209;
t195 = -t196 * t207 + t200 * t212;
t194 = t204 * t220 + t221 * t217;
t192 = t202 * t220 + t222 * t217;
t190 = -t193 * t207 + t199 * t212;
t189 = -t191 * t207 + t198 * t212;
t1 = [0, 0, t199, t190, t194 * t216 - t224 * t219 (t194 * t219 + t224 * t216) * t215 - t190 * t218; 0, 0, t198, t189, t192 * t216 - t225 * t219 (t192 * t219 + t225 * t216) * t215 - t189 * t218; 0, 0, t200, t195, t197 * t216 - t223 * t219 (t197 * t219 + t223 * t216) * t215 - t195 * t218;];
Jg_rot  = t1;
