% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR7_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobig_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:31
% EndTime: 2019-02-26 19:57:31
% DurationCPUTime: 0.12s
% Computational Cost: add. (102->33), mult. (296->71), div. (0->0), fcn. (409->16), ass. (0->45)
t231 = sin(pkin(13));
t234 = sin(pkin(6));
t257 = t231 * t234;
t233 = sin(pkin(7));
t256 = t233 * t234;
t239 = cos(pkin(6));
t255 = t233 * t239;
t236 = cos(pkin(13));
t254 = t236 * t234;
t238 = cos(pkin(7));
t245 = cos(qJ(2));
t253 = t238 * t245;
t242 = sin(qJ(2));
t252 = t239 * t242;
t251 = t239 * t245;
t227 = t231 * t245 + t236 * t252;
t230 = sin(pkin(14));
t235 = cos(pkin(14));
t226 = -t231 * t242 + t236 * t251;
t247 = t226 * t238 - t233 * t254;
t217 = -t227 * t230 + t235 * t247;
t223 = -t226 * t233 - t238 * t254;
t232 = sin(pkin(8));
t237 = cos(pkin(8));
t250 = t217 * t237 + t223 * t232;
t229 = -t231 * t252 + t236 * t245;
t228 = -t231 * t251 - t236 * t242;
t246 = t228 * t238 + t231 * t256;
t219 = -t229 * t230 + t235 * t246;
t224 = -t228 * t233 + t238 * t257;
t249 = t219 * t237 + t224 * t232;
t221 = t235 * t255 + (-t230 * t242 + t235 * t253) * t234;
t225 = t239 * t238 - t245 * t256;
t248 = t221 * t237 + t225 * t232;
t244 = cos(qJ(4));
t243 = cos(qJ(5));
t241 = sin(qJ(4));
t240 = sin(qJ(5));
t222 = t234 * t242 * t235 + (t234 * t253 + t255) * t230;
t220 = t229 * t235 + t230 * t246;
t218 = t227 * t235 + t230 * t247;
t216 = -t221 * t232 + t225 * t237;
t215 = -t219 * t232 + t224 * t237;
t214 = -t217 * t232 + t223 * t237;
t1 = [0, t257, 0, t215, t220 * t241 - t244 * t249 (t220 * t244 + t241 * t249) * t240 - t215 * t243; 0, -t254, 0, t214, t218 * t241 - t244 * t250 (t218 * t244 + t241 * t250) * t240 - t214 * t243; 0, t239, 0, t216, t222 * t241 - t244 * t248 (t222 * t244 + t241 * t248) * t240 - t216 * t243;];
Jg_rot  = t1;
