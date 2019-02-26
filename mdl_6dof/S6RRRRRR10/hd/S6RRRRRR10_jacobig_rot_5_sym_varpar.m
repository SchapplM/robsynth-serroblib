% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR10_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobig_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:48
% EndTime: 2019-02-26 22:52:48
% DurationCPUTime: 0.09s
% Computational Cost: add. (55->27), mult. (161->57), div. (0->0), fcn. (227->14), ass. (0->35)
t241 = sin(pkin(7));
t245 = cos(pkin(6));
t263 = t241 * t245;
t244 = cos(pkin(7));
t252 = cos(qJ(2));
t262 = t244 * t252;
t242 = sin(pkin(6));
t249 = sin(qJ(1));
t261 = t249 * t242;
t248 = sin(qJ(2));
t260 = t249 * t248;
t259 = t249 * t252;
t253 = cos(qJ(1));
t258 = t253 * t242;
t257 = t253 * t248;
t256 = t253 * t252;
t236 = t245 * t256 - t260;
t255 = t236 * t244 - t241 * t258;
t238 = -t245 * t259 - t257;
t254 = t238 * t244 + t241 * t261;
t251 = cos(qJ(3));
t250 = cos(qJ(4));
t247 = sin(qJ(3));
t246 = sin(qJ(4));
t243 = cos(pkin(8));
t240 = sin(pkin(8));
t239 = -t245 * t260 + t256;
t237 = t245 * t257 + t259;
t235 = -t242 * t252 * t241 + t245 * t244;
t234 = -t238 * t241 + t244 * t261;
t233 = -t236 * t241 - t244 * t258;
t232 = t251 * t263 + (-t247 * t248 + t251 * t262) * t242;
t231 = -t239 * t247 + t254 * t251;
t230 = -t237 * t247 + t255 * t251;
t1 = [0, t261, t234, -t231 * t240 + t234 * t243 (t239 * t251 + t254 * t247) * t246 + (-t231 * t243 - t234 * t240) * t250, 0; 0, -t258, t233, -t230 * t240 + t233 * t243 (t237 * t251 + t255 * t247) * t246 + (-t230 * t243 - t233 * t240) * t250, 0; 1, t245, t235, -t232 * t240 + t235 * t243 (t247 * t263 + (t247 * t262 + t248 * t251) * t242) * t246 + (-t232 * t243 - t235 * t240) * t250, 0;];
Jg_rot  = t1;
