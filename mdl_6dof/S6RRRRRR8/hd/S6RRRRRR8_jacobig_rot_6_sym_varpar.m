% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR8
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
% Datum: 2019-02-26 22:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR8_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:51:24
% EndTime: 2019-02-26 22:51:24
% DurationCPUTime: 0.07s
% Computational Cost: add. (50->22), mult. (131->45), div. (0->0), fcn. (189->12), ass. (0->34)
t220 = sin(pkin(7));
t223 = cos(pkin(6));
t239 = t220 * t223;
t222 = cos(pkin(7));
t228 = cos(qJ(2));
t238 = t222 * t228;
t221 = sin(pkin(6));
t226 = sin(qJ(1));
t237 = t226 * t221;
t225 = sin(qJ(2));
t236 = t226 * t225;
t235 = t226 * t228;
t229 = cos(qJ(1));
t234 = t229 * t221;
t233 = t229 * t225;
t232 = t229 * t228;
t213 = t223 * t232 - t236;
t231 = -t213 * t222 + t220 * t234;
t215 = -t223 * t235 - t233;
t230 = t215 * t222 + t220 * t237;
t227 = cos(qJ(3));
t224 = sin(qJ(3));
t219 = qJ(4) + qJ(5);
t218 = cos(t219);
t217 = sin(t219);
t216 = -t223 * t236 + t232;
t214 = t223 * t233 + t235;
t212 = -t221 * t228 * t220 + t223 * t222;
t211 = -t215 * t220 + t222 * t237;
t210 = -t213 * t220 - t222 * t234;
t209 = -t227 * t239 + (t224 * t225 - t227 * t238) * t221;
t208 = t216 * t224 - t230 * t227;
t207 = t214 * t224 + t231 * t227;
t1 = [0, t237, t211, t208, t208 (t216 * t227 + t230 * t224) * t217 - t211 * t218; 0, -t234, t210, t207, t207 (t214 * t227 - t231 * t224) * t217 - t210 * t218; 1, t223, t212, t209, t209 (t224 * t239 + (t224 * t238 + t225 * t227) * t221) * t217 - t212 * t218;];
Jg_rot  = t1;
