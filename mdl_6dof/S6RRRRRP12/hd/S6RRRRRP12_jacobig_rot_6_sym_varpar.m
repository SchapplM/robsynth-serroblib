% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRP12_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:46:23
% EndTime: 2019-02-26 22:46:23
% DurationCPUTime: 0.06s
% Computational Cost: add. (34->21), mult. (100->45), div. (0->0), fcn. (145->12), ass. (0->30)
t226 = sin(pkin(7));
t229 = cos(pkin(6));
t247 = t226 * t229;
t228 = cos(pkin(7));
t236 = cos(qJ(2));
t246 = t228 * t236;
t227 = sin(pkin(6));
t233 = sin(qJ(1));
t245 = t233 * t227;
t232 = sin(qJ(2));
t244 = t233 * t232;
t243 = t233 * t236;
t237 = cos(qJ(1));
t242 = t237 * t227;
t241 = t237 * t232;
t240 = t237 * t236;
t222 = t229 * t240 - t244;
t239 = -t222 * t228 + t226 * t242;
t224 = -t229 * t243 - t241;
t238 = t224 * t228 + t226 * t245;
t235 = cos(qJ(3));
t234 = cos(qJ(4));
t231 = sin(qJ(3));
t230 = sin(qJ(4));
t225 = -t229 * t244 + t240;
t223 = t229 * t241 + t243;
t221 = -t227 * t236 * t226 + t229 * t228;
t220 = -t224 * t226 + t228 * t245;
t219 = -t222 * t226 - t228 * t242;
t1 = [0, t245, t220, t225 * t231 - t238 * t235 (t225 * t235 + t238 * t231) * t230 - t220 * t234, 0; 0, -t242, t219, t223 * t231 + t239 * t235 (t223 * t235 - t239 * t231) * t230 - t219 * t234, 0; 1, t229, t221, -t235 * t247 + (t231 * t232 - t235 * t246) * t227 (t231 * t247 + (t231 * t246 + t232 * t235) * t227) * t230 - t221 * t234, 0;];
Jg_rot  = t1;
