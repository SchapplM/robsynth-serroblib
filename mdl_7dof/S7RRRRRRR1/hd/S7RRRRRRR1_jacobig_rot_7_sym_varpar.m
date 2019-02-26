% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 7 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% Jg_rot [3x7]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S7RRRRRRR1_jacobig_rot_7_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobig_rot_7_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobig_rot_7_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_7_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:31
% EndTime: 2019-02-26 22:54:31
% DurationCPUTime: 0.07s
% Computational Cost: add. (39->23), mult. (101->43), div. (0->0), fcn. (153->12), ass. (0->30)
t254 = sin(qJ(3));
t255 = sin(qJ(2));
t269 = t255 * t254;
t260 = cos(qJ(3));
t268 = t255 * t260;
t256 = sin(qJ(1));
t267 = t256 * t255;
t261 = cos(qJ(2));
t266 = t256 * t261;
t262 = cos(qJ(1));
t265 = t262 * t254;
t264 = t262 * t255;
t263 = t262 * t260;
t259 = cos(qJ(4));
t258 = cos(qJ(5));
t257 = cos(qJ(6));
t253 = sin(qJ(4));
t252 = sin(qJ(5));
t251 = sin(qJ(6));
t250 = -t256 * t254 + t261 * t263;
t249 = -t256 * t260 - t261 * t265;
t248 = t260 * t266 + t265;
t247 = -t254 * t266 + t263;
t246 = -t261 * t253 + t259 * t268;
t245 = t253 * t268 + t261 * t259;
t244 = t250 * t259 + t253 * t264;
t243 = t250 * t253 - t259 * t264;
t242 = t248 * t259 + t253 * t267;
t241 = t248 * t253 - t259 * t267;
t1 = [0, t256, -t264, t249, t243, t244 * t252 - t249 * t258 -(t244 * t258 + t249 * t252) * t251 + t243 * t257; 0, -t262, -t267, t247, t241, t242 * t252 - t247 * t258 -(t242 * t258 + t247 * t252) * t251 + t241 * t257; 1, 0, t261, -t269, t245, t246 * t252 + t258 * t269 -(t246 * t258 - t252 * t269) * t251 + t245 * t257;];
Jg_rot  = t1;
