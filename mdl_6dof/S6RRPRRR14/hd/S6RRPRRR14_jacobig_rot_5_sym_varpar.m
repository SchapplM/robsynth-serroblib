% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR14_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobig_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:46
% EndTime: 2019-02-26 22:55:46
% DurationCPUTime: 0.09s
% Computational Cost: add. (50->27), mult. (146->57), div. (0->0), fcn. (204->14), ass. (0->35)
t222 = sin(pkin(7));
t227 = cos(pkin(6));
t243 = t222 * t227;
t226 = cos(pkin(7));
t232 = cos(qJ(2));
t242 = t226 * t232;
t223 = sin(pkin(6));
t230 = sin(qJ(1));
t241 = t230 * t223;
t229 = sin(qJ(2));
t240 = t230 * t229;
t239 = t230 * t232;
t233 = cos(qJ(1));
t238 = t233 * t223;
t237 = t233 * t229;
t236 = t233 * t232;
t216 = t227 * t236 - t240;
t235 = t216 * t226 - t222 * t238;
t218 = -t227 * t239 - t237;
t234 = t218 * t226 + t222 * t241;
t231 = cos(qJ(4));
t228 = sin(qJ(4));
t225 = cos(pkin(8));
t224 = cos(pkin(14));
t221 = sin(pkin(8));
t220 = sin(pkin(14));
t219 = -t227 * t240 + t236;
t217 = t227 * t237 + t239;
t215 = -t223 * t232 * t222 + t227 * t226;
t214 = -t218 * t222 + t226 * t241;
t213 = -t216 * t222 - t226 * t238;
t212 = t224 * t243 + (-t220 * t229 + t224 * t242) * t223;
t211 = -t219 * t220 + t234 * t224;
t210 = -t217 * t220 + t235 * t224;
t1 = [0, t241, 0, -t211 * t221 + t214 * t225 (t219 * t224 + t234 * t220) * t228 + (-t211 * t225 - t214 * t221) * t231, 0; 0, -t238, 0, -t210 * t221 + t213 * t225 (t217 * t224 + t235 * t220) * t228 + (-t210 * t225 - t213 * t221) * t231, 0; 1, t227, 0, -t212 * t221 + t215 * t225 (t223 * t229 * t224 + (t223 * t242 + t243) * t220) * t228 + (-t212 * t225 - t215 * t221) * t231, 0;];
Jg_rot  = t1;
